function window_values(n,kind)
    require(kind in ("none","rectangular","hann","hamming","blackman"),"Invalid window")
    [kind=="hann" ? .5-.5*cos(TAU*i/n) : kind=="hamming" ? .54-.46*cos(TAU*i/n) : kind=="blackman" ? .42-.5*cos(TAU*i/n)+.08*cos(2TAU*i/n) : 1.0 for i in 0:n-1]
end
onesided(k,n) = (k==1 || (iseven(n) && k==n÷2+1)) ? 1.0 : 2.0
function amplitude_spectrum(x,fs; window="hann", padding="none")
    require(finitepositive(fs),"Invalid sample frequency")
    require(padding in ("none","next","double"),"Invalid padding")
    n=length(x); require(n>=8,"At least 8 samples required")
    m=padding=="none" ? n : padding=="next" ? nextpow(2,n) : Base.checked_mul(2,nextpow(2,n))
    w=window_values(n,window); data=zeros(m)
    @inbounds for i in 1:n; data[i]=x[i]*w[i]; end
    f=rfft(data); sw=sum(w); values=[onesided(k,m)*abs(f[k])/sw for k in eachindex(f)]
    require(all(isfinite,values),"FFT overflow; rescale input")
    (;frequency=collect(0:length(f)-1).*(fs/m),values,nfft=m,bin_width=fs/m,resolution=fs/n,coherent_gain=sw/n)
end
function segment_settings(n,segment_length,overlap)
    require(segment_length isa Integer && segment_length>=8,"Segment length must be an integer >=8")
    require(isfinite(overlap) && 0<=overlap<1,"Overlap must be in [0,1)")
    len=min(n,segment_length)
    # Match JavaScript Math.round for positive half-integers (Julia round is ties-to-even).
    hop=max(1,floor(Int,len*(1-overlap)+.5))
    (;length=len,step=hop,segments=(n-len)÷hop+1)
end
"Call sink(time, frequency, PSD) once per full frame; PSD buffer is reused, copy only if retaining it."
function stft_each(sink,t,x,fs; segment_length=512,overlap=.75,window="hann")
    require(length(t)==length(x) && length(x)>=8,"Invalid time/value lengths")
    require(finitepositive(fs),"Invalid sample frequency")
    s=segment_settings(length(x),segment_length,overlap); w=window_values(s.length,window)
    input=zeros(s.length); output=zeros(ComplexF64,s.length÷2+1); values=zeros(length(output))
    plan=plan_rfft(input;flags=FFTW.ESTIMATE)
    scale=sqrt(fs)*sqrt(sum(abs2,w)); frequency=collect(0:length(values)-1).*(fs/s.length)
    for j in 0:s.segments-1
        start=j*s.step+1
        @inbounds for i in 1:s.length; input[i]=x[start+i-1]*w[i]; end
        mul!(output,plan,input)
        @inbounds for k in eachindex(values); values[k]=onesided(k,s.length)*abs2(output[k]/scale); end
        require(all(isfinite,values),"PSD overflow; rescale input")
        sink((t[start]+t[start+s.length-1])/2,frequency,values)
    end
    (;s...,frequency,bin_width=fs/s.length)
end
function welch_psd(x,fs; segment_length=1024,overlap=.5,window="hann")
    s=segment_settings(length(x),segment_length,overlap); values=zeros(s.length÷2+1)
    stft_each(range(0,step=1/fs,length=length(x)),x,fs;segment_length,overlap,window) do _,_,p
        @inbounds for k in eachindex(values); values[k]+=p[k]/s.segments; end
    end
    (;frequency=collect(0:length(values)-1).*(fs/s.length),values,segment_length=s.length,segments=s.segments,step=s.step,bin_width=fs/s.length)
end
function find_peaks(amplitude,psd; threshold=.02,distance=0.0,count=8)
    require(isfinite(threshold) && 0<=threshold<=1 && isfinite(distance) && distance>=0 && count isa Integer && count>=1,"Invalid peak settings")
    v=amplitude.values; f=amplitude.frequency; maxv=maximum(@view v[2:end])
    candidates=[i for i in 2:length(v) if v[i]>=threshold*maxv && v[i]>v[i-1] && (i==length(v) || v[i]>=v[i+1])]
    sort!(candidates;by=i->(-v[i],i)); chosen=NamedTuple[]
    for i in candidates
        if all(p->abs(p.frequency-f[i])>=distance,chosen)
            z=f[i]/psd.bin_width; k=floor(Int,z)+1
            p=k>=length(psd.values) ? psd.values[end] : psd.values[k]+(psd.values[k+1]-psd.values[k])*(z-(k-1))
            push!(chosen,(frequency=f[i],period=1/f[i],amplitude=v[i],psd=p))
        end
        length(chosen)>=count && break
    end
    chosen
end
