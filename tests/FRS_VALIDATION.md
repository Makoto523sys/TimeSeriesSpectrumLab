# FRS numerical validation

10/10 checks passed. Tested the engine embedded in index.html using v24.19.0.

| Check | Result | Detail |
|---|---|---|
| Zero acceleration gives zero response for every measure and damping | PASS |  |
| Undamped resonant sine agrees with closed-form transient maxima | PASS | Sa relative error 8.226e-6 |
| Damped step agrees with analytic overshoot including near-critical damping | PASS |  |
| Independent RK4 reference for broadband piecewise-linear input and free tail | PASS | maximum relative error 4.299e-4 |
| Linearity, sign invariance and all supported input unit conversions | PASS |  |
| Absolute and pseudo acceleration differ at finite damping | PASS |  |
| Free vibration captures delayed peak; nonzero final sample resets equilibrium | PASS |  |
| Internal time-step refinement converges on the same sampled input | PASS | 20-step error 4.053e-1; 100-step error 9.751e-3 |
| FRS uses preprocessed acceleration without FFT window; integrated input rejected | PASS |  |
| Invalid ranges, damping, units, nonfinite data and work limits fail explicitly | PASS |  |

Reproduce: `node tests/frs.test.cjs`. References are analytic transient solutions and an independent fine-step RK4 integrator, not another call to Newmark. Numerical verification is not validation of input records, seismic design compliance, or full browser compatibility.
