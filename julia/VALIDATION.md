# Julia implementation validation

Executed 2026-09-19 with 64-bit Julia 1.10.10, FFTW.jl 1.10.0, Linux, 4 Julia threads. These are numerical/CLI tests, not seismic-design certification or a prediction of performance on other hardware.

## Correctness

`julia --project=julia --threads=4 julia/test/runtests.jl`: **121 / 121 assertions passed**.

| Area | Assertions |
|---|---:|
| CSV/quoted header/whitespace parser, generated time, invalid rows, sampling and >524,288 points | 11 |
| FFT odd/even/arbitrary lengths, windows, DC/Nyquist, PSD energy, STFT, peaks, padding | 65 |
| Detrend, resampling, all three FIR filters, first/second integration | 10 |
| FRS analytic resonance/step/pulse response, units, deterministic threads, checkpoint reuse/partial completion/input mismatch/corruption | 21 |
| Complete pipeline, CSV and binary STFT, HTML/SVG, resume, integration and config errors | 14 |

`julia --project=julia --threads=4 julia/test/browser_parity.jl`: **132 / 132 assertions passed**.

This invokes the current repository's actual `index.html` numerical engine with Node and compares the Julia results for six cases: odd/even FFT lengths, four windows, double padding, all three FIR filters, detrending, one/two integrations, irregular sampling/resampling, Welch, STFT, peaks and FRS with 0.5%, 5%, 20% damping. Comparison tolerances are `rtol=5e-10`, `atol=1e-12`. Maximum absolute difference across history/FFT/PSD/STFT arrays and their tested axes was **5.035e-15**. FRS comparisons passed the same tolerance. This verifies numerical consistency, not equivalence of web UI behavior.

The initial parity tests detected a Julia closure-binding mistake in the band-pass kernel. It was fixed before the final passing tests. A standalone FIR pass-/stop-band test also covers this behavior.

## Actual workload above 100 million steps

Reproduce:

```sh
julia --project=julia --threads=4 julia/scripts/benchmark_frs.jl
```

- Nonzero synthetic acceleration: 1,000,001 samples at dt=0.001 s, with 3 Hz and 11 Hz components.
- 64 oscillator frequencies from 1 to 20 Hz, 2% and 5% damping, 100 minimum steps/period and 5 residual periods.
- **158,215,814 Newmark integration steps actually executed per run**. Both serial and threaded calculations completed.
- Serial: **3.196 s**; four threads: **1.044 s**.
- Serial/threaded spectra were exactly equal.
- Allocations measured inside the warmed FRS calls: approximately **8.10 MB** each, primarily the input acceleration buffer. Memory does not grow with integration-step count.
- Process peak RSS approximately **428.7 MB** includes the Julia runtime, compilation, original input and both calculations. It is not the FRS working-set allocation alone.
- An additional 100→200 steps/period comparison at 3 and 11 Hz on the same million-sample input changed Sa by at most **1.834e-5 relative**.

Raw measurement: [validation/frs-large.toml](validation/frs-large.toml). Compilation warm-up is excluded from timed FRS calls; CSV loading, checkpoint writes and report generation are not included. No acceleration histories were zeroed or skipped to shorten the benchmark.

## Large FFT/Welch/STFT

Reproduce:

```sh
julia --project=julia --threads=4 julia/scripts/benchmark_spectra.jl
```

| Calculation | Workload | Measured time |
|---|---|---:|
| Arbitrary-length real FFT | 1,000,001 samples | 0.321 s |
| Welch PSD | 20,058,112 segment-sample visits | 0.073 s |
| Streamed STFT callback | 19,588 frames / 10,048,644 PSD cells | 0.232 s |

FFT bin-centred amplitudes and Welch/STFT energy consistency passed. The STFT callback accumulated energy without retaining the full matrix. These timings exclude output IO; STFT callback compilation is included. Production CLI streams full STFT to binary or CSV. Binary file size, layout and values were checked in the end-to-end tests.

Raw measurement: [validation/spectra-large.toml](validation/spectra-large.toml).

## CLI example and output checks

`julia/run.jl julia/examples/full-analysis.toml` completed the 4,000-point example, including 600 FRS oscillator jobs / 8,615,031 steps. All 24 files listed by its completed run manifest existed; all 14 exported SVG files parsed as valid XML. The report uses bounded display envelopes and a bounded STFT preview; CSV/binary numerical outputs are complete.

Browser rendering and native Windows execution were not performed. This environment's cloud browser blocks local `file://` navigation. Generated HTML was checked through content assertions and SVG syntax; numerical and IO checks ran in Julia. Early development runs used `--compiled-modules=no` to avoid transient concurrent precompile locks; the final 121-assertion suite passed with normal compiled modules enabled.
