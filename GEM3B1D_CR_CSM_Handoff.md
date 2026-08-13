# GEM3B1D CR/CSM Handoff Notes

Date: 2026-08-07

## Goal Completed
Implemented and stabilized complex-ranged support in GEM3B1D, including mixed CR+CSM mode, and validated the reduced 3B->2B spectator regression behavior.

## Main Issues Found and Resolved

1. CR bra-side consistency in matrix elements:
- Cause: bra-side range/norm values were not consistently conjugated in CR mode.
- Fix: use conjugated bra-side values in GEM3B1D matrix-element assembly.

2. Struct mismatch and field usage in fill path:
- Cause: temporary argument handling still used tuple-like assumptions after struct changes.
- Fix: use `TempArgs1D` field access consistently.

3. CR+CSM `BoundsError` in matrix-fill preparation:
- Cause: temporary argument array was preallocated for lower-triangular fill, but CR+CSM path uses full fill.
- Fix: allocate full `nbasis_total^2` capacity when full fill is required.

4. Baseline no-CR regression after introducing coordinate-level CR controls:
- Cause: coordinate-level CR flags could leak into non-CR runs.
- Fix: gate per-coordinate CR flags behind `complex_ranged=true` so default/no-CR behavior is unchanged.

5. Reduced-limit CR spectrum duplication concern:
- Observation: duplication was associated with applying CR in both Jacobi coordinates for the spectator setup.
- Resolution approach: allow coordinate-selective CR and use CR on `r` only in the reduced spectator regression.

## Functional Change Introduced

Added coordinate-selective complex-range control in GEM3B1D:
- `complex_ranged_r`
- `complex_ranged_R`

Behavior:
- If `complex_ranged=false`: both coordinates remain real-ranged.
- If `complex_ranged=true` and no coordinate flags provided: both coordinates use CR (backward-compatible).
- If `complex_ranged=true` and flags are provided: CR is applied per coordinate according to the flags.

## Files Updated

1. `src/GEM-3body-1D/GEM3B1D.jl`
- Threaded coordinate-level CR booleans through solver pipeline.
- Added safe gating so coordinate flags do not affect no-CR runs.

2. `src/GEM-3body-1D/preallocate.jl`
- Type and effective basis sizing now depend on `complex_ranged_r` and `complex_ranged_R`.
- Full temporary-argument capacity for CR+CSM full-fill mode.

3. `src/GEM-3body-1D/precomputation.jl`
- Range construction split by coordinate (`nu` vs `NU`) for selective CR.

4. `src/GEM-3body-1D/fillTVS.jl`
- Fill path now receives coordinate-level CR flags.
- `fill_full` logic aligned with mixed complex mode behavior.

5. `test/test3B1D.jl`
- Spectator-limit regression remains strict pairwise comparison:
  - `e2` vs `e3`
  - `e2cr` vs `e3cr`
  - `e2csmcr` vs `e3csmcr`
- Reduced-limit numerical setup uses `complex_ranged_r=true`, `complex_ranged_R=false`.

## Current Validation Snapshot (from user runs)

Representative outputs now show close agreement:
- no-CR: `e3` and `e2` match within tolerance
- CR: `e3cr` and `e2cr` match within tolerance
- CR+CSM: `e3csmcr` and `e2csmcr` match within tolerance, with small expected imaginary parts

## Suggested Commit Scope

Include:
- `src/GEM-3body-1D/GEM3B1D.jl`
- `src/GEM-3body-1D/preallocate.jl`
- `src/GEM-3body-1D/precomputation.jl`
- `src/GEM-3body-1D/fillTVS.jl`
- `test/test3B1D.jl`
- `GEM3B1D_CR_CSM_Handoff.md`

## Next Session: ISGL Port Plan

Because ISGL mirrors GEM3B1D structure, apply the same pattern:

1. API wiring:
- Add coordinate-level CR controls and thread through `size_estimate`, preallocation, precompute, and fill.

2. Type/effective-size consistency:
- Ensure basis lengths and matrix types follow coordinate-level CR state.

3. Full-fill logic for mixed complex mode:
- Recheck temporary-loop cardinality and matrix-fill assumptions under CR+CSM.

4. Validation:
- Add reduced-limit 3D-vs-2B checks in `test/testISGL.jl` for:
  - baseline
  - CR only
  - CSM only
  - CR+CSM

5. Guardrails:
- Keep no-CR default path behavior unchanged when CR options are not enabled.
