# Benchmarks

Time and peak memory of the solvers, to compare versions of the package on the same machine.
Not part of the test suite.

```
julia benchmark/run.jl 0.5.4              # a registered release
julia benchmark/run.jl /path/to/checkout   # a local checkout (records commit and dirty flag)
julia benchmark/compare.jl benchmark/results/<old>.csv benchmark/results/<new>.csv
julia benchmark/compare.jl benchmark/results/<a>.csv <b>.csv <c>.csv ...   # overview, >2 files
```

Options of `run.jl`: `--quick` (smallest size only), `--cases=isgl_cr,inverse` (subset).

- **Cases** (`cases.jl`): ISGL with real ranges, complex scaling and complex ranges; GEM3B1D with real
  and complex ranges; `inverse_solve` on ISGL matrices. Each at three sizes (two for `inverse`).
  Only the public API is used. Releases >= 0.3.1 are supported: the keywords of the old API
  (< 0.4.0) are translated in `cases.jl`, and a case is skipped for releases older than its
  `since` (complex ranges ISGL 0.5.1, GEM3B1D 0.5.2, `inverse_solve` 0.5.4). Releases < 0.5.4
  have no `*_matrices`, so `t_build` is missing. To change a case, add a new name instead,
  otherwise old results are no longer comparable.
- **Measurement** (`worker.jl`): every point runs in a fresh julia process, after a warm-up at a tiny
  size (compilation). `t_solve` is a full `*_solve` (energies only), `t_build` a `*_matrices` call
  (matrix build only). `dpeak_*` is the peak resident memory (VmHWM, reset before each call) above
  the resident memory before the call, in bytes. `E1..E3` are the lowest energies (or strengths for
  `inverse`), for a consistency check between versions.
- **Results** (`results/`): `<label>.csv` (one row per point), `<label>.toml` (package version and
  commit, julia version, CPU, julia and BLAS threads, load average), `<label>_Manifest.toml`,
  `<label>.log`. The label is `<version or commit>_<date>`. No host name or local paths are recorded.
- Timings are only comparable on the same machine with the same thread settings and an otherwise
  idle machine (check `loadavg_start` and the log). A full run takes about 35 min on a
  2x Xeon E5-2630 v3 (16 cores) and needs up to ~12 GB (v0.5.4, `isgl_cr` at N = 6912).
- Known differences: CSM energies of v0.3.1 are wrong (sign error, fixed by v0.3.9).
