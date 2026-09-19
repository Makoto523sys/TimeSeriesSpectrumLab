# DOM event integration validation

Executed 2026-09-19T02:13:17.192Z, Node v24.19.0. 16/16 checks passed against the scripts embedded in final `index.html`.

This is **not real-browser GUI verification**. It executes the actual application and numerical scripts in a VM using minimal DOM/canvas mocks. Worker messages execute in real Node worker threads using the actual Blob worker source. Browser local-file permissions, pixel rendering, browser download behavior, accessibility, and cross-browser compatibility remain unverified. No browser or alternative graphical surface was used.

| Check | Result |
|---|---|
| Load built-in sample through parse Worker and compute through analysis Worker | PASS |
| PSD toggle, FFT/PSD CSV payloads, peak and time-history exports | PASS |
| Selected value column changes detected frequency and blocks stale exports | PASS |
| Parameter change marks result dirty and export handler refuses stale data | PASS |
| Nonuniform input warns and blocks analysis until resampling enabled | PASS |
| Tab/no-header input, non-first time column, millisecond conversion | PASS |
| Malformed and missing selected numeric fields produce explicit errors | PASS |
| Invalid numeric settings and frequency range are reported | PASS |
| Cancel terminates real worker; late completion cannot restore result | PASS |
| Integration metadata separates processed and original physical units on dual axes | PASS |
| Chart event handlers zoom, pan, reset and populate hover text | PASS |
| FRS requires explicit acceleration units and rejects time integration | PASS |
| FRS Worker produces multiple damping curves, SI output and complete CSV metadata | PASS |
| FRS response measure, period ordering and log axis switch without recalculation | PASS |
| FRS setting changes block stale CSV and PNG, disabled FRS ignores invalid fields | PASS |
| Cancelling FRS terminates Worker and prevents stale result restoration | PASS |


Reproduce with `node tests/ui.test.cjs`. Canvas operations are mocked; chart tests verify event-driven state and metadata, not visible appearance.
