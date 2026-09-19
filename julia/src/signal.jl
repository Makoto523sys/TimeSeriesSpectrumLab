# Parse one physical line at a time; only selected numerical columns are retained.
# Like the HTML parser, multiline quoted fields are explicitly unsupported.
function split_fields(line, delimiter)
    delimiter == "space" && return String.(split(line))
    sep = delimiter == "comma" ? ',' : '\t'
    fields=String[]; field=IOBuffer(); quoted=false; chars=collect(line); i=1
    while i <= length(chars)
        c=chars[i]
        if c == '"'
            if quoted && i < length(chars) && chars[i+1] == '"'
                write(field,'"'); i+=1
            else
                quoted=!quoted
            end
        elseif c == sep && !quoted
            push!(fields,strip(String(take!(field))))
        else
            write(field,c)
        end
        i+=1
    end
    require(!quoted,"Unclosed quote (multiline CSV is not supported)")
    push!(fields,strip(String(take!(field)))); fields
end
function read_signal(path; delimiter="auto", header="auto", time_column=1,
                     value_column=2, time_scale=1.0, dt=nothing, skip_lines=0)
    require(delimiter in ("auto","comma","tab","space"),"Invalid delimiter")
    require(header in ("auto","yes","no"),"Invalid header setting")
    require(time_column isa Integer && time_column >= 0 && value_column isa Integer && value_column > 0,"Columns use positive 1-based indices; time_column=0 synthesizes time")
    require(time_column != value_column,"Time and value columns must differ")
    require(finitepositive(time_scale),"time_scale must be positive")
    require(skip_lines isa Integer && skip_lines >= 0,"skip_lines must be nonnegative")
    time_column == 0 && require(dt !== nothing && finitepositive(dt),"dt in seconds is required when time_column=0")
    t=Float64[]; x=Float64[]; names=String[]; width=0; firstrow=true
    for (lineno, raw) in enumerate(eachline(path))
        lineno <= skip_lines && continue
        line=strip(replace(raw,'\ufeff'=>""))
        (isempty(line) || startswith(line,"#")) && continue
        if firstrow && delimiter == "auto"
            delimiter=occursin(',',line) ? "comma" : occursin('\t',line) ? "tab" : "space"
        end
        fields=split_fields(line,delimiter)
        if firstrow
            width=length(fields)
            require(max(time_column,value_column) <= width,"Selected column exceeds file width")
            hasheader=header == "yes" || (header == "auto" && any(s->tryparse(Float64,s)===nothing,fields))
            names=hasheader ? fields : ["Column $i" for i in 1:width]
            firstrow=false
            hasheader && continue
        end
        require(length(fields)==width,"Line $lineno: inconsistent column count")
        v=tryparse(Float64,fields[value_column])
        tv=time_column==0 ? length(x)*Float64(dt) : tryparse(Float64,fields[time_column])
        require(v!==nothing && isfinite(v),"Line $lineno: value is not finite numeric data")
        require(tv!==nothing && isfinite(tv),"Line $lineno: time is not finite numeric data")
        time_column!=0 && (tv*=time_scale)
        require(isfinite(tv) && (isempty(t) || tv>t[end]),"Line $lineno: time must be finite and strictly increasing")
        push!(t,tv); push!(x,v)
    end
    require(length(x)>=8,"At least 8 numerical rows are required")
    (;t,x,headers=names,delimiter)
end
function inspect_sampling(t)
    n=length(t); require(n>=8,"At least 8 samples are required")
    require(all(isfinite,t),"Times must be finite")
    mind=Inf; maxd=-Inf
    @inbounds for i in 2:n
        d=t[i]-t[i-1]; require(isfinite(d) && d>0,"Time must strictly increase")
        mind=min(mind,d); maxd=max(maxd,d)
    end
    duration=t[end]-t[1]; dt=duration/(n-1); fs=1/dt
    require(finitepositive(fs),"Invalid sample frequency")
    (;n,start=t[1],stop=t[end],duration,dt,min_dt=mind,max_dt=maxd,fs,nyquist=fs/2,resolution=fs/n,jitter=(maxd-mind)/dt)
end
function signal_stats(x)
    require(!isempty(x),"Empty signal")
    mean=0.0; m2=0.0; sq=0.0; lo=Inf; hi=-Inf
    for (i,v) in enumerate(x)
        require(isfinite(v),"Nonfinite value at sample $i")
        d=v-mean; mean+=d/i; m2+=d*(v-mean); sq+=v*v; lo=min(lo,v); hi=max(hi,v)
    end
    require(all(isfinite,(mean,m2,sq)),"Signal exceeds floating-point statistics range; rescale input")
    (;min=lo,max=hi,mean,std=sqrt(max(0,m2/length(x))),rms=sqrt(sq/length(x)))
end
function preprocess(x; remove_mean=true, detrend=false)
    y=Float64.(x); n=length(y)
    if remove_mean || detrend
        m=signal_stats(y).mean; y .-= m
    end
    if detrend
        center=(n-1)/2; numerator=0.0; denominator=0.0
        @inbounds for i in eachindex(y)
            q=i-1-center; numerator+=q*y[i]; denominator+=q*q
        end
        slope=numerator/denominator
        @inbounds for i in eachindex(y); y[i]-=slope*(i-1-center); end
    end
    y
end
function resample_signal(t,x)
    require(length(t)==length(x),"Time/value lengths differ")
    info=inspect_sampling(t); out=similar(x,Float64); times=collect(range(t[1],t[end],length=length(t))); j=1
    @inbounds for i in eachindex(times)
        while j<length(t)-1 && t[j+1]<times[i]; j+=1; end
        out[i]=x[j]+(x[j+1]-x[j])*(times[i]-t[j])/(t[j+1]-t[j])
    end
    times,out
end
function filter_signal(x,fs; kind="none", low_cut=1.0, high_cut=20.0)
    require(finitepositive(fs),"Invalid sampling frequency")
    kind=="none" && return Float64.(x)
    require(kind in ("lowpass","highpass","bandpass"),"Invalid filter")
    require((kind=="highpass" || 0<high_cut<fs/2) && (kind=="lowpass" || 0<low_cut<fs/2) && (kind!="bandpass" || low_cut<high_cut),"Cutoffs must satisfy 0 < low < high < Nyquist")
    n=length(x); require(n>=8,"At least 8 samples required")
    half=min(128,(n-1)÷2)
    function lowkernel(cut)
        kernel=[(j==0 ? 2cut/fs : sin(TAU*cut/fs*j)/(pi*j))*(.54+.46*cos(pi*j/half)) for j in -half:half]
        kernel ./= sum(kernel); kernel
    end
    h=kind=="highpass" ? lowkernel(low_cut) : lowkernel(high_cut)
    if kind=="highpass"; h .*= -1; h[half+1]+=1; end
    kind=="bandpass" && (h .-= lowkernel(low_cut))
    y=zeros(n)
    Threads.@threads for i in 1:n
        s=0.0
        @inbounds for j in -half:half
            k=i+j; k<1 && (k=2-k); k>n && (k=2n-k)
            s+=h[j+half+1]*x[k]
        end
        y[i]=s
    end
    y
end
function integrate_signal(t,x,count)
    require(count in (0,1,2),"Integration count must be 0, 1 or 2")
    require(length(t)==length(x),"Time/value lengths differ")
    y=Float64.(x)
    for _ in 1:count
        total=0.0; previous=y[1]; y[1]=0.0
        @inbounds for i in 2:length(y)
            current=y[i]; total+=(current+previous)*.5*(t[i]-t[i-1]); y[i]=total; previous=current
        end
    end
    y
end
