const DEFAULTS = Dict(
 "input"=>Dict("path"=>"", "delimiter"=>"auto", "header"=>"auto", "time_column"=>1,"value_column"=>2,"time_scale"=>1.0,"dt"=>0.01,"skip_lines"=>0,"quantity"=>"Signal","unit"=>"a.u."),
 "preprocess"=>Dict("tolerance"=>1e-4,"resample"=>false,"remove_mean"=>true,"detrend"=>false,"filter"=>"none","low_cut"=>1.0,"high_cut"=>20.0,"integrate"=>0),
 "fft"=>Dict("enabled"=>true,"window"=>"hann","padding"=>"none"),
 "welch"=>Dict("enabled"=>true,"segment_length"=>1024,"overlap"=>.5,"window"=>"hann"),
 "peaks"=>Dict("enabled"=>true,"threshold"=>.02,"distance"=>.5,"count"=>8),
 "stft"=>Dict("enabled"=>true,"segment_length"=>512,"overlap"=>.75,"window"=>"hann","format"=>"binary"),
 "frs"=>Dict("enabled"=>false,"minimum"=>.1,"maximum"=>50.0,"count"=>200,"spacing"=>"log","damping_percent"=>[2.0,5.0],"steps_per_period"=>100,"residual_cycles"=>5.0,"unit"=>"","parallel"=>true,"checkpoint_directory"=>""),
 "output"=>Dict("directory"=>"results","overwrite"=>false,"history"=>true,"report"=>true)
)
function complete_config(raw)
    config=deepcopy(DEFAULTS)
    for (section,values) in raw
        require(haskey(config,section),"Unknown config section: $section")
        require(values isa AbstractDict,"Config section $section must be a table")
        for (k,v) in values
            require(haskey(config[section],k),"Unknown setting: $section.$k")
            config[section][k]=v
        end
    end
    for (section,values) in config, (k,v) in values
        default=DEFAULTS[section][k]
        if default isa Bool
            require(v isa Bool,"$section.$k must be true or false")
        elseif default isa Integer
            require(v isa Integer && !(v isa Bool),"$section.$k must be an integer")
        elseif default isa AbstractFloat
            require(v isa Real && !(v isa Bool) && isfinite(v),"$section.$k must be finite numeric data")
        elseif default isa String
            require(v isa String,"$section.$k must be a string")
        end
    end
    config
end
csvcell(v) = begin
    s=string(v)
    v isa AbstractString && occursin(r"^[=+\-@]",s) && (s="'"*s)
    occursin(r"[,\"\r\n]",s) ? "\""*replace(s,"\""=>"\"\"")*"\"" : s
end
csvrow(io,values)=println(io,join(csvcell.(values),','))
function atomic_file(f,path)
    tmp=path*".partial"; open(f,tmp,"w"); mv(tmp,path;force=true)
end
function write_columns(path,headers,columns)
    require(all(c->length(c)==length(columns[1]),columns),"Output column lengths differ")
    atomic_file(path) do io
        csvrow(io,headers)
        for i in eachindex(columns[1]); csvrow(io,(c[i] for c in columns)); end
    end
end
asdict(x::NamedTuple)=Dict(string(k)=>v for (k,v) in pairs(x))
function export_stft(t,x,fs,c,out)
    s=segment_settings(length(x),c["segment_length"],c["overlap"])
    require(c["format"] in ("binary","csv"),"STFT format must be binary or csv")
    bins=s.length÷2+1; nx=min(s.segments,160); ny=min(bins,120)
    preview=fill(-300.0,ny,nx); frame=0
    format=c["format"]; target=joinpath(out,format=="binary" ? "stft.f64" : "stft.csv")
    @info "STFT output" frames=s.segments bins binary_bytes=big(8)*bins*s.segments format
    atomic_file(joinpath(out,"stft-times.csv")) do times_io
        csvrow(times_io,["frame_1based","time_s"])
        atomic_file(target) do io
            format=="csv" && csvrow(io,["time_s","frequency_Hz","PSD"])
            stft_each(t,x,fs;segment_length=c["segment_length"],overlap=c["overlap"],window=c["window"]) do time,f,p
                frame+=1; csvrow(times_io,[frame,time])
                if format=="binary"
                    write(io,p) # native-endian Float64, frame-major (Julia column-major bins × frames)
                else
                    for k in eachindex(p); csvrow(io,[time,f[k],p[k]]); end
                end
                ix=min(nx,(frame-1)*nx÷s.segments+1)
                for k in eachindex(p)
                    iy=min(ny,(k-1)*ny÷bins+1)
                    preview[iy,ix]=max(preview[iy,ix],10log10(max(1e-30,p[k])))
                end
            end
        end
    end
    frequency=collect(0:bins-1).*(fs/s.length)
    write_columns(joinpath(out,"stft-frequencies.csv"),["frequency_Hz"],[frequency])
    meta=Dict("format"=>format,"dtype"=>"Float64","endian"=>(Base.ENDIAN_BOM==0x04030201 ? "little" : "big"),"layout"=>"frequency bins contiguous within each frame; bins x frames column-major","bins"=>bins,"frames"=>s.segments,"segment_length"=>s.length,"hop_samples"=>s.step,"window"=>c["window"],"scale"=>"PSD; input_unit^2/Hz","preview"=>"maximum dB within each time-frequency display block; full output is not decimated")
    atomic_toml(joinpath(out,"stft-metadata.toml"),meta)
    (;preview,frequency,frames=s.segments,start=(t[1]+t[s.length])/2,stop=(t[(s.segments-1)*s.step+1]+t[(s.segments-1)*s.step+s.length])/2)
end
function run_analysis(config_path::AbstractString)
    config=complete_config(TOML.parsefile(config_path)); base=dirname(abspath(config_path))
    for (section,key) in [("input","path"),("output","directory"),("frs","checkpoint_directory")]
        val=config[section][key]; isempty(val) || (config[section][key]=abspath(base,val))
    end
    run_analysis(config)
end
"Execute numerical processing from a config dictionary; paths in dictionaries are relative to pwd()."
function run_analysis(raw::AbstractDict)
    config=complete_config(raw); inp=config["input"]; pre=config["preprocess"]; fc=config["frs"]; oc=config["output"]
    require(!isempty(inp["path"]),"input.path is required")
    require(pre["tolerance"]>=0 && pre["integrate"] in (0,1,2),"Invalid preprocessing settings")
    require(!(fc["enabled"] && pre["integrate"]!=0),"FRS requires acceleration; set preprocess.integrate=0")
    require(!config["peaks"]["enabled"] || (config["fft"]["enabled"] && config["welch"]["enabled"]),"Peak detection requires FFT and Welch")
    input_path=abspath(inp["path"]); out=abspath(oc["directory"])
    require(input_path!=out && !startswith(input_path,out*(Sys.iswindows() ? "\\" : "/")),"Input must be outside output directory")
    require(!isdir(out) || isempty(readdir(out)) || oc["overwrite"],"Output directory is not empty. Choose another directory or explicitly set output.overwrite=true (also used when resuming).")
    mkpath(out)
    runlock=joinpath(out,"analysis.lock"); require(!ispath(runlock),"Another run or stale lock exists: $runlock"); mkdir(runlock)
    metadata=Dict{String,Any}("status"=>"running","config"=>config,"julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"output_files"=>String[])
    start=time()
    try
        metadata["input_sha256"]=open(sha256,input_path) |> bytes2hex
        atomic_toml(joinpath(out,"run.toml"),metadata)
        @info "Reading numerical columns" input=input_path
        data=read_signal(input_path;delimiter=inp["delimiter"],header=inp["header"],time_column=inp["time_column"],value_column=inp["value_column"],time_scale=inp["time_scale"],dt=inp["dt"],skip_lines=inp["skip_lines"])
        t=data.t; x=data.x; initial=merge(inspect_sampling(t),signal_stats(x)); warnings=String[]
        tolerance=fc["enabled"] ? min(pre["tolerance"],1e-4) : pre["tolerance"]
        if initial.jitter>tolerance
            require(pre["resample"],"Nonuniform sampling: enable preprocess.resample explicitly")
            t,x=resample_signal(t,x); push!(warnings,"Linear resampling applied; no reconstruction of missing data or anti-alias filtering.")
        end
        sampling=inspect_sampling(t); fs=sampling.fs
        before=preprocess(x;remove_mean=pre["remove_mean"],detrend=pre["detrend"])
        processed=filter_signal(before,fs;kind=pre["filter"],low_cut=pre["low_cut"],high_cut=pre["high_cut"])
        if pre["filter"]!="none"; push!(warnings,"Centered <=257-tap FIR with reflected edges; check edge transients and transition band."); end
        pre["integrate"]!=0 && (processed=integrate_signal(t,processed,pre["integrate"]); push!(warnings,"Zero-initial-value trapezoidal integration; drift is not seismic baseline correction."))
        processed_stats=signal_stats(processed); metadata["input_summary"]=asdict(initial); metadata["processed_summary"]=merge(asdict(sampling),asdict(processed_stats)); metadata["headers"]=data.headers
        outputs=metadata["output_files"]
        if oc["history"]
            write_columns(joinpath(out,"history.csv"),["time_s","input_resampled_if_requested","processed"],[t,x,processed]); push!(outputs,"history.csv")
        end
        fft=nothing; psd=nothing; oldfft=nothing; oldpsd=nothing; peaks=NamedTuple[]
        compare=pre["filter"]!="none" || pre["integrate"]!=0
        c=config["fft"]
        if c["enabled"]
            @info "FFT" samples=length(processed)
            fft=amplitude_spectrum(processed,fs;window=c["window"],padding=c["padding"])
            write_columns(joinpath(out,"fft.csv"),["frequency_Hz","amplitude"],[fft.frequency,fft.values]); push!(outputs,"fft.csv")
            if compare
                oldfft=amplitude_spectrum(before,fs;window=c["window"],padding=c["padding"])
                write_columns(joinpath(out,"fft-before.csv"),["frequency_Hz","amplitude_before_filter_integration"],[oldfft.frequency,oldfft.values]); push!(outputs,"fft-before.csv")
            end
        end
        c=config["welch"]
        if c["enabled"]
            @info "Welch PSD"
            psd=welch_psd(processed,fs;segment_length=c["segment_length"],overlap=c["overlap"],window=c["window"])
            write_columns(joinpath(out,"welch.csv"),["frequency_Hz","PSD"],[psd.frequency,psd.values]); push!(outputs,"welch.csv")
            if compare
                oldpsd=welch_psd(before,fs;segment_length=c["segment_length"],overlap=c["overlap"],window=c["window"])
                write_columns(joinpath(out,"welch-before.csv"),["frequency_Hz","PSD_before_filter_integration"],[oldpsd.frequency,oldpsd.values]); push!(outputs,"welch-before.csv")
            end
        end
        if config["peaks"]["enabled"]
            c=config["peaks"]; peaks=find_peaks(fft,psd;threshold=c["threshold"],distance=c["distance"],count=c["count"])
            atomic_file(joinpath(out,"peaks.csv")) do io
                csvrow(io,["frequency_Hz","period_s","amplitude","Welch_PSD_interpolated"])
                for p in peaks; csvrow(io,values(p)); end
            end
            push!(outputs,"peaks.csv")
        end
        stft=nothing
        if config["stft"]["enabled"]
            stft=export_stft(t,processed,fs,config["stft"],out)
            append!(outputs,[config["stft"]["format"]=="binary" ? "stft.f64" : "stft.csv","stft-times.csv","stft-frequencies.csv","stft-metadata.toml"])
        end
        frs=nothing
        if fc["enabled"]
            require(fc["damping_percent"] isa AbstractVector && all(z->z isa Real && isfinite(z),fc["damping_percent"]),"damping_percent must be a numeric list")
            o=FRSOptions(minimum=fc["minimum"],maximum=fc["maximum"],count=fc["count"],spacing=fc["spacing"],damping=Float64.(fc["damping_percent"])./100,steps_per_period=fc["steps_per_period"],residual_cycles=fc["residual_cycles"],unit=fc["unit"])
            estimated=estimate_frs(length(processed),sampling.dt,o)
            @info "FRS workload (no 100-million-step ceiling)" integration_steps=string(estimated.work) jobs=o.count*length(o.damping) threads=Threads.nthreads()
            checkpoint=isempty(fc["checkpoint_directory"]) ? joinpath(out,"frs-checkpoints") : fc["checkpoint_directory"]
            last_report=Ref(time()); progress=(done,total)->begin
                if done==total || done==0 || time()-last_report[]>2
                    @info "FRS progress" completed=done total; last_report[]=time()
                end
            end
            frs=response_spectrum(processed,sampling.dt;options=o,parallel=fc["parallel"],checkpoint_directory=checkpoint,progress)
            append!(warnings,frs.warnings)
            atomic_file(joinpath(out,"frs.csv")) do io
                csvrow(io,["frequency_Hz","period_s","damping_ratio","Sa_absolute_m_per_s2","PSa_pseudo_m_per_s2","Sd_relative_m","Sv_relative_m_per_s","PSv_pseudo_m_per_s","substeps"])
                for r in frs.rows; csvrow(io,(getfield(r,k) for k in fieldnames(FRSRow))); end
            end
            push!(outputs,"frs.csv"); metadata["frs"]=Dict("pfa_m_per_s2"=>frs.pfa,"integration_steps"=>string(frs.work),"resumed_jobs"=>frs.resumed,"checkpoint_directory"=>abspath(checkpoint),"fingerprint"=>frs.fingerprint,"method"=>"Newmark beta=1/4 gamma=1/2; zero initial relative displacement and velocity; piecewise-linear input; zero floor acceleration after record")
        end
        metadata["warnings"]=warnings
        if oc["report"]
            write_report(joinpath(out,"report.html"),t,x,processed,fft,psd,oldfft,oldpsd,peaks,stft,frs,metadata)
            push!(outputs,"report.html")
        end
        metadata["elapsed_seconds"]=time()-start; metadata["status"]="complete"
        atomic_toml(joinpath(out,"run.toml"),metadata)
        @info "Analysis complete" output=out elapsed=metadata["elapsed_seconds"]
        (;output=out,metadata,fft,psd,peaks,frs)
    catch e
        metadata["status"]="incomplete"; metadata["error"]=sprint(showerror,e); metadata["elapsed_seconds"]=time()-start
        atomic_toml(joinpath(out,"run.toml"),metadata)
        rethrow()
    finally
        rm(runlock)
    end
end
