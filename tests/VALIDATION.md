# Numerical validation

Executed 2026-09-19T02:14:03.816Z using v24.19.0, linux/x64.

Source tested: `index.html#numerical-engine`.

22/22 checks passed. Tests use independent direct DFT references only for small vectors; production uses FFT.

| Check | Result | Detail |
|---|---|---|
| Independent DFT comparison: lengths 1–35, 63, 64, 65, 97 | PASS | maximum absolute complex-component error 6.040e-13 |
| CSV, tab, whitespace, quoted headers, no header, selected columns | PASS |  |
| Sampling metadata and population statistics | PASS |  |
| Known three-tone amplitude and frequency: all four windows | PASS |  |
| DC and even Nyquist are not doubled; odd final bin is doubled | PASS |  |
| Mean removal and least-squares linear detrend | PASS |  |
| Periodic window coherent gains and correction | PASS |  |
| Zero padding changes grid, preserves normalization and native resolution | PASS |  |
| PSD Parseval identity for odd/even lengths and every window | PASS |  |
| Welch segment averaging, overlap, tonal power, and PSD units normalization | PASS |  |
| Welch arbitrary signal equals independently averaged window-weighted power | PASS |  |
| Finite-data and settings failures are explicit | PASS |  |
| Large representable PSD avoids raw FFT-square overflow; invalid fs is rejected | PASS |  |
| Local peaks, period conversion, DC exclusion and separation | PASS |  |
| Nonuniform sampling blocks FFT; linear resampling preserves linear signal | PASS |  |
| STFT chirp ridge tracks instantaneous frequency and segment centers | PASS | maximum chirp ridge error 1.9531 Hz (bin width 3.90625 Hz) |
| Centered FIR low/high/band filters reject stop bands with interior phase retained | PASS |  |
| Trapezoidal integration: constant acceleration, zero initial conditions | PASS |  |
| All shipped example files parse and analyze as intended | PASS |  |
| Performance 10,000 samples: FFT + Welch + STFT | PASS | 37.4 ms; 18 STFT frames |
| Performance 100,000 samples: FFT + Welch + STFT | PASS | 363.8 ms; 194 STFT frames |
| Performance 300,000 samples: FFT + Welch + STFT | PASS | 1329.0 ms; 584 STFT frames |

## Performance

Wall-clock measurements are from this development container, not a browser guarantee. Each includes full-record arbitrary-length FFT, overlapping Welch PSD, peak detection, and STFT; timings exclude file parsing, Web Worker transfer, rendering, and user interaction.

| Samples | Elapsed (ms) | STFT frames |
|---:|---:|---:|
| 10000 | 37.4 | 18 |
| 100000 | 363.8 | 194 |
| 300000 | 1329.0 | 584 |

## Reproduce

Run `node tests/numerical.test.cjs` from a checkout. Node is used only for developer verification; end users open `index.html` directly. Regenerate sample files using `node tests/generate-examples.cjs`.

Tests prefer the shipped `index.html` numerical-engine script and otherwise use `dev/engine.js`. GUI/browser verification is documented separately, if performed. Population standard deviation is tested; segment detrending is not assumed by Welch tests. PSD integrals use discrete bin sums multiplied by bin width (no half-weight endpoints).
