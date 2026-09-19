using Test, TimeSeriesSpectrumLab, TOML
fixtures=TOML.parse(read(`node $(joinpath(@__DIR__,"browser-fixture.cjs"))`,String))
maxerr=0.0
@testset "Browser engine parity" begin
    for c in values(fixtures)
        t=Float64.(c["t"]); x=Float64.(c["x"])
        c["resample"] && ((t,x)=resample_signal(t,x))
        fs=inspect_sampling(t).fs
        y=preprocess(x;remove_mean=c["remove_mean"],detrend=c["detrend"])
        y=filter_signal(y,fs;kind=c["filter"],low_cut=c["low_cut"],high_cut=c["high_cut"])
        y=integrate_signal(t,y,c["integrate"])
        a=amplitude_spectrum(y,fs;window=c["window"],padding="double")
        p=welch_psd(y,fs;segment_length=33,overlap=.5,window=c["window"])
        stft=Float64[]; times=Float64[]
        stft_each(t,y,fs;segment_length=32,overlap=.5,window="hann") do time,_,v; append!(stft,v); push!(times,time); end
        for (actual,key) in [(y,"processed"),(a.values,"fft"),(a.frequency,"frequency"),(p.values,"psd"),(stft,"stft"),(times,"stft_times")]
            @test actual≈c[key] rtol=5e-10 atol=1e-12
            global maxerr=max(maxerr,maximum(abs,actual.-c[key]))
        end
        peaks=find_peaks(a,p;distance=1.)
        @test length(peaks)==length(c["peaks"])
        for (actual,reference) in zip(peaks,c["peaks"]); @test collect(values(actual))≈reference rtol=5e-10 atol=1e-12; end
        if !isempty(c["frs"])
            o=FRSOptions(minimum=1.,maximum=25.,count=5,damping=[.005,.05,.2],steps_per_period=100,residual_cycles=2.,unit="gal")
            r=response_spectrum(y,1/fs;options=o)
            for (row,ref) in zip(r.rows,c["frs"])
                @test [row.sa,row.psa,row.sd,row.sv,row.psv]≈ref rtol=5e-10 atol=1e-12
            end
        end
    end
end
println("Maximum absolute error (history/FFT/frequencies/PSD/STFT): ",maxerr)
