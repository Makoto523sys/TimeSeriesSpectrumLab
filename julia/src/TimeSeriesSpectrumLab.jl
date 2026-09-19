module TimeSeriesSpectrumLab
using FFTW, LinearAlgebra, Printf, SHA, TOML
export read_signal, inspect_sampling, signal_stats, preprocess, resample_signal,
       filter_signal, integrate_signal, amplitude_spectrum, welch_psd, stft_each,
       find_peaks, FRSOptions, response_spectrum, estimate_frs, run_analysis
const TAU = 2pi
require(ok, message) = ok || throw(ArgumentError(message))
finitepositive(v) = v isa Real && isfinite(v) && v > 0
include("signal.jl")
include("spectra.jl")
include("frs.jl")
include("pipeline.jl")
include("report.jl")
end
