# Nonzero input exceeding BOTH browser record length and FRS workload caps.
# Warm-up is excluded; compare serial/threaded values and measure allocations/RSS.
using TimeSeriesSpectrumLab, TOML, Test
output=isempty(ARGS) ? joinpath(@__DIR__,"..","results","large-validation.toml") : ARGS[1]
n=1_000_001; dt=.001
x=[sin(2pi*3.0*i*dt)+.25sin(2pi*11.0*i*dt) for i in 0:n-1]
o=FRSOptions(minimum=1.,maximum=20.,count=64,damping=[.02,.05],steps_per_period=100,residual_cycles=5.)
response_spectrum(x[1:100],dt;options=FRSOptions(minimum=1.,maximum=5.,count=2),parallel=false)
response_spectrum(x[1:100],dt;options=FRSOptions(minimum=1.,maximum=5.,count=2),parallel=true)
est=estimate_frs(n,dt,o)
@test est.work>100_000_000
println("Samples: ",n,"; FRS integration steps: ",est.work,"; threads: ",Threads.nthreads())
GC.gc(); serial=@timed response_spectrum(x,dt;options=o,parallel=false)
GC.gc(); parallel=@timed response_spectrum(x,dt;options=o,parallel=true)
@test serial.value.rows==parallel.value.rows
@test all(r->isfinite(r.sa)&&r.sa>0,parallel.value.rows)
# An independent time-step refinement check on the same million-point record.
checkopts=FRSOptions(minimum=3.,maximum=11.,count=2,damping=[.05],steps_per_period=200,residual_cycles=5.)
refined=response_spectrum(x,dt;options=checkopts,parallel=true)
coarse=response_spectrum(x,dt;options=FRSOptions(minimum=3.,maximum=11.,count=2,damping=[.05],steps_per_period=100,residual_cycles=5.),parallel=true)
errors=[abs(a.sa/b.sa-1) for (a,b) in zip(coarse.rows,refined.rows)]
@test maximum(errors)<.02
report=Dict("julia"=>string(VERSION),"os"=>string(Sys.KERNEL),"cpu"=>Sys.CPU_NAME,"threads"=>Threads.nthreads(),"samples"=>n,"integration_steps"=>string(est.work),"serial_seconds"=>serial.time,"threaded_seconds"=>parallel.time,"serial_allocated_bytes"=>serial.bytes,"threaded_allocated_bytes"=>parallel.bytes,"process_peak_rss_bytes"=>Sys.maxrss(),"serial_threaded_identical"=>true,"max_refinement_relative_change"=>maximum(errors),"scope"=>"FRS computation only; input generated in memory; compilation warm-up excluded; RSS includes Julia runtime, compilation and both runs; not a performance guarantee")
mkpath(dirname(abspath(output))); open(output,"w") do io; TOML.print(io,report;sorted=true); end
TOML.print(stdout,report;sorted=true)
