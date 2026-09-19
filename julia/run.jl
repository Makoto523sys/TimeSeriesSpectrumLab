#!/usr/bin/env julia
using TimeSeriesSpectrumLab
if length(ARGS)!=1 || ARGS[1] in ("-h","--help")
    println("Usage: julia --project=julia --threads=auto julia/run.jl path/to/config.toml")
    println("Paths inside TOML are relative to the configuration file. See julia/README.md.")
    exit(isempty(ARGS) || ARGS[1] in ("-h","--help") ? 0 : 2)
end
try
    run_analysis(ARGS[1])
catch e
    showerror(stderr,e); println(stderr)
    exit(e isa InterruptException ? 130 : 1)
end
