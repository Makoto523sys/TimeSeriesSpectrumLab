# Exercise arbitrary-length FFT, Welch above 16M sample-visits, and streamed STFT.
using TimeSeriesSpectrumLab, TOML, Test
output=isempty(ARGS) ? joinpath(@__DIR__,"..","results","spectra-validation.toml") : ARGS[1]
n=1_000_001; fs=1000.0; t=collect(0:n-1)./fs
# Exactly bin-centred sine even at the odd arbitrary FFT length.
x=[sin(2pi*3000i/n)+.3cos(2pi*12000i/n) for i in 0:n-1]
amplitude_spectrum(x[1:129],fs); welch_psd(x[1:129],fs)
fft=@timed amplitude_spectrum(x,fs)
@test fft.value.nfft==n
@test fft.value.values[3001]≈1.0 rtol=1e-10
@test fft.value.values[12001]≈.3 rtol=1e-10
welch=@timed welch_psd(x,fs;segment_length=1024,overlap=.95)
@test welch.value.segments*welch.value.segment_length>16_000_000
frames=Ref(0); energy=Ref(0.)
stft=@timed stft_each(t,x,fs;segment_length=1024,overlap=.95) do _,_,p
    frames[]+=1; energy[]+=sum(p)*(fs/1024)
end
@test frames[]==welch.value.segments
@test energy[]/frames[]≈sum(welch.value.values)*welch.value.bin_width rtol=1e-10
@test length(stft.value.frequency)*frames[]>1_500_000
report=Dict("julia"=>string(VERSION),"samples"=>n,"fft_seconds"=>fft.time,"welch_seconds"=>welch.time,"stft_seconds"=>stft.time,"fft_allocated_bytes"=>fft.bytes,"welch_allocated_bytes"=>welch.bytes,"stft_allocated_bytes"=>stft.bytes,"welch_sample_visits"=>welch.value.segments*welch.value.segment_length,"stft_cells"=>length(stft.value.frequency)*frames[],"stft_frames"=>frames[],"note"=>"In-memory input; STFT callback accumulates energy without retaining the matrix; excludes output IO; warm-up excluded for FFT/Welch, STFT callback compilation included")
mkpath(dirname(abspath(output))); open(output,"w") do io; TOML.print(io,report;sorted=true); end
TOML.print(stdout,report;sorted=true)
