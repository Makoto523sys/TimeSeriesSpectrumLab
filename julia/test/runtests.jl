using Test, TimeSeriesSpectrumLab, TOML
const Lab=TimeSeriesSpectrumLab
@testset "Parser and sampling" begin
    mktempdir() do d
        p=joinpath(d,"data.csv")
        write(p,"# metadata\n\ufeff\"time\",\"a,x\",other\n"*join(["$(i*10),$(i+1),text" for i in 0:15],"\n"))
        r=read_signal(p;time_scale=.001)
        @test r.headers==["time","a,x","other"]
        @test inspect_sampling(r.t).dt≈.01
        @test r.x==collect(1.0:16)
        write(p,"meta\nNS EW UD\n"*join(["$i $(2i) $(3i)" for i in 1:16],"\n"))
        r=read_signal(p;skip_lines=1,time_column=0,value_column=2,dt=.01)
        @test r.t≈collect(0:.01:.15)
        @test r.x==collect(2.0:2:32)
        write(p,"t,x\n"*join(["$i,$(i==5 ? "" : string(i))" for i in 1:10],"\n"))
        @test_throws ArgumentError read_signal(p)
        @test_throws ArgumentError Lab.split_fields("\"unclosed,x","comma")
    end
    @test_throws ArgumentError inspect_sampling([0,1,2,3,3,5,6,7])
    @test_throws ArgumentError signal_stats([1,NaN])
    @test signal_stats([1,2,3]).std≈sqrt(2/3)
    # Browser size limit is removed, even when no actual spectral calculation is requested.
    @test inspect_sampling(range(0,step=.001,length=600_000)).n==600_000
end
@testset "FFT, Welch, STFT, peaks" begin
    for n in (128,129,257), window in ("rectangular","hann","hamming","blackman")
        fs=128.0; t=collect(0:n-1)./fs; x=@. 2sin(2pi*5*128/n*t)
        a=amplitude_spectrum(x,fs;window)
        @test a.values[6]≈2 atol=2e-12
        p=welch_psd(x,fs;segment_length=n,overlap=0.,window)
        w=Lab.window_values(n,window)
        @test sum(p.values)*p.bin_width≈sum(abs2,x.*w)/sum(abs2,w) rtol=1e-12
        frames=[]
        stft_each(t,x,fs;segment_length=n,overlap=0.,window) do tc,f,v
            push!(frames,(tc,copy(f),copy(v)))
        end
        @test frames[1][3]≈p.values
        @test frames[1][1]≈(t[1]+t[end])/2
        padded=amplitude_spectrum(x,fs;window,padding="double")
        @test padded.nfft==2nextpow(2,n)
    end
    x=[(-1.0)^i for i in 0:127]
    @test amplitude_spectrum(x,128.;window="rectangular").values[end]≈1
    @test amplitude_spectrum(fill(3.,128),128.;window="rectangular").values[1]≈3
    x=[sin(2pi*2i/1000)+.5sin(2pi*10i/1000)+.2sin(2pi*25i/1000) for i in 0:3999]
    a=amplitude_spectrum(x,1000.); p=welch_psd(x,1000.)
    @test [r.frequency for r in find_peaks(a,p;distance=1.)][1:3]==[2.,10.,25.]
    @test Lab.segment_settings(20,9,.5).step==5
    @test_throws ArgumentError welch_psd(x,1000.;overlap=1.)
end
@testset "Preprocess, resampling, filters, integration" begin
    t=collect(0:.01:2); x=@. 2+3t
    @test maximum(abs,preprocess(x;detrend=true))<1e-13
    @test integrate_signal(t,fill(2.,length(t)),1)≈2t
    @test integrate_signal(t,fill(2.,length(t)),2)≈t.^2
    irregular=@. t+.0002sin(3t); values=@. 2+3irregular
    rt,rx=resample_signal(irregular,values)
    @test rx≈(2 .+ 3rt)
    n=8192; fs=1000.; mixed=[sin(2pi*40i/fs)+sin(2pi*200i/fs) for i in 0:n-1]
    for (kind,lo,hi,pass,stop) in [("lowpass",0.,100.,40.,200.),("highpass",100.,0.,200.,40.),("bandpass",150.,250.,200.,40.)]
        y=filter_signal(mixed,fs;kind,low_cut=lo,high_cut=hi)
        project(f)=2abs(sum(y[i]*cis(2pi*f*(i-1)/fs) for i in 513:n-512))/(n-1024)
        @test abs(project(pass)-1)<.01
        @test project(stop)<.01
    end
end
@testset "FRS analytic, units, deterministic threading and checkpoints" begin
    t=collect(0:.0005:2); w=4pi; x=sin.(w.*t)
    o=FRSOptions(minimum=2.,maximum=3.,count=2,damping=[0.],residual_cycles=0.)
    r=response_spectrum(x,.0005;options=o,parallel=false)
    sd=maximum(abs,@. t*cos(w*t)/(2w)-sin(w*t)/(2w*w))
    @test r.rows[1].sd≈sd rtol=2e-4
    @test r.rows[1].sa≈w*w*sd rtol=2e-4
    @test r.rows[1].sa≈r.rows[1].psa
    @test response_spectrum(x,.0005;options=o,parallel=true).rows==r.rows
    for (unit,scale) in [("g",9.80665),("gal",.01),("mm/s2",.001)]
        ou=FRSOptions(minimum=2.,maximum=3.,count=2,damping=[0.],residual_cycles=0.,unit=unit)
        ru=response_spectrum(-3x./scale,.0005;options=ou)
        @test ru.rows[1].sa≈3r.rows[1].sa rtol=1e-11
    end
    for z in (.05,.2,.9)
        os=FRSOptions(minimum=2.,maximum=3.,count=2,damping=[z],residual_cycles=0.)
        rs=response_spectrum(ones(20001),.0005;options=os)
        @test rs.rows[1].psa≈1+exp(-pi*z/sqrt(1-z*z)) rtol=2e-4
    end
    tail=FRSOptions(minimum=1.,maximum=2.,count=2,damping=[0.],residual_cycles=2.)
    @test response_spectrum(ones(51),.001;options=tail).rows[1].sa≈2sin(pi*.05) rtol=2e-4
    mktempdir() do d
        a=response_spectrum(x,.0005;options=o,checkpoint_directory=d)
        b=response_spectrum(x,.0005;options=o,checkpoint_directory=d)
        @test b.resumed==2
        @test a.rows==b.rows
        rm(joinpath(d,"2.toml"))
        @test response_spectrum(x,.0005;options=o,checkpoint_directory=d).resumed==1
        @test_throws ArgumentError response_spectrum(2x,.0005;options=o,checkpoint_directory=d)
        @test !ispath(joinpath(d,"run.lock"))
        broken=TOML.parsefile(joinpath(d,"1.toml")); broken["sa"]=NaN
        open(joinpath(d,"1.toml"),"w") do io; TOML.print(io,broken); end
        @test_throws ArgumentError response_spectrum(x,.0005;options=o,checkpoint_directory=d)
    end
    @test estimate_frs(1_000_001,.001,FRSOptions(minimum=1.,maximum=20.,count=64,damping=[.02,.05],residual_cycles=0.)).work>100_000_000
    @test_throws ArgumentError estimate_frs(50,.01,FRSOptions(maximum=51.))
    @test_throws ArgumentError estimate_frs(50,.01,FRSOptions(damping=[1.]))
    @test_throws ArgumentError response_spectrum([0.,NaN],.001)
end
@testset "End-to-end artifacts, restart and explicit validation" begin
    mktempdir() do d
        input=joinpath(d,"data.csv"); out=joinpath(d,"out")
        t=collect(0:127).*.001; x=sin.(2pi*20 .* t)
        open(input,"w") do io
            println(io,"time,acc"); for i in eachindex(t); println(io,t[i],',',x[i]); end
        end
        cfg=Dict("input"=>Dict("path"=>input,"unit"=>"m/s2"),"output"=>Dict{String,Any}("directory"=>out),"frs"=>Dict("enabled"=>true,"unit"=>"m/s2","minimum"=>1.,"maximum"=>40.,"count"=>6,"residual_cycles"=>1.),"stft"=>Dict{String,Any}("segment_length"=>32,"overlap"=>.5))
        r=run_analysis(cfg)
        @test r.metadata["status"]=="complete"
        @test isfile(joinpath(out,"report.html"))
        @test occursin("Sa absolute",read(joinpath(out,"report.html"),String))
        @test isfile(joinpath(out,"frs-sa.svg"))
        meta=TOML.parsefile(joinpath(out,"stft-metadata.toml"))
        @test filesize(joinpath(out,"stft.f64"))==8meta["bins"]*meta["frames"]
        binary=open(joinpath(out,"stft.f64")) do io; read!(io,Vector{Float64}(undef,meta["bins"]*meta["frames"])); end
        expected=[]
        stft_each(t,preprocess(x),1000.;segment_length=32,overlap=.5) do _,_,p; append!(expected,p); end
        @test binary≈expected
        @test_throws ArgumentError run_analysis(cfg)
        cfg["output"]["overwrite"]=true
        @test run_analysis(cfg).frs.resumed==12
        cfg["stft"]["format"]="csv"
        @test run_analysis(cfg).metadata["status"]=="complete"
        @test countlines(joinpath(out,"stft.csv"))==1+meta["bins"]*meta["frames"]
        cfg["preprocess"]=Dict("integrate"=>1)
        @test_throws ArgumentError run_analysis(cfg)
        cfg["frs"]["enabled"]=false
        @test run_analysis(cfg).metadata["status"]=="complete"
        @test isfile(joinpath(out,"fft-before.csv"))
        @test_throws ArgumentError Lab.complete_config(Dict("frs"=>Dict("typo"=>true)))
    end
end
