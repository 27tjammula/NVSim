# FeDiode NVSim Methods Delta (vs. UCSB NVSim, informed by Fe-NVSim)

## Scope
This document summarizes the simulator extensions made on top of baseline NVSim (UCSB) to support a two-terminal MFIS ferroelectric diode (FeDiode, AlScN/HfOx stack), and clarifies what was simulated vs. what was anchored to measured data.

## 1) Type-System and Parser Extensions

### New memory-cell type
- Added `FeDiode` to `MemCellType` in [`typedef.h`](./typedef.h).

### New FeDiode cell parameters
- Added to [`MemCell.h`](./MemCell.h):
  - `capacitanceFeDiode`
  - `capacitanceFeDiodeReverse`
  - `polarizationRemnant` (`Pr`)
  - `polarizationSpontaneous` (`Ps`)
  - `coerciveField` (`Ec`)
  - `ferroelectricThickness` (`tFE`)
  - `interlayerThickness` (`tIL`)
  - `interlayerPermittivity` (`eIL`)
  - `ferroelectricPermittivity` (`eFE`)
  - `eta`
  - `ferroelectricMaterial`

### Input parsing behavior
- Added parsing of FeDiode tokens/fields in [`MemCell.cpp`](./MemCell.cpp).
- On `-MemCellType: FeDiode`, access-drop is forced to zero:
  - `voltageDropAccessDevice = 0`
  - Rationale: no separate transistor selector; device is two-terminal.

## 2) Subarray Electrical Model Changes for FeDiode

Implemented in [`SubArray.cpp`](./SubArray.cpp):

### Access-device RC handling
- For FeDiode, selector RC no longer uses NMOS gate/drain formulas.
- Selected-cell access capacitance uses:
  - `capCellAccess = capacitanceFeDiode`
- Unselected loading uses reverse-biased capacitance:
  - `capWordline += capacitanceFeDiodeReverse * numColumn`
  - `capBitline  += capacitanceFeDiodeReverse * numRow`
- FeDiode `none_access` path is explicitly treated as self-rectifying diode behavior (not `MAX(capOn, capOff)` fallback).

### Sneak/leakage model for crossbar writes
- Added FeDiode-specific worst-case half-select leakage:
  - same-WL unselected paths and same-BL unselected paths use reverse-biased term `Vwrite/2 / Roff`.
- This replaces memristor-style ohmic half-select assumptions in FeDiode branch.

### FeDiode validity checks
- Extended internal-sensing mux guard to include FeDiode (`muxSenseAmp < 2` condition).
- Rectification feasibility check for FeDiode uses:
  - `Roff / Ron` against array-depth requirement (`numRow / tolerance`).

## 3) Write-Energy Modeling Strategy

### Device files provide explicit energies
- For FeDiode, per-cell write energies are provided directly in `.cell` files, rather than relying on generic access-drop-based write-energy formulas.
- Three files:
  - [`1FeD-AlScN-100nm.cell`](./1FeD-AlScN-100nm.cell)
  - [`1FeD-AlScN-200nm.cell`](./1FeD-AlScN-200nm.cell)
  - [`1FeD-AlScN-1um.cell`](./1FeD-AlScN-1um.cell)

### Energy decomposition used
- Switching term:
  - `E_switch = 2 * Pr * Ec * tFE * A`
- Junction charging term:
  - `E_charge = C_junction * Vwrite^2`
- Total:
  - `E_total = E_switch + E_charge`

## 4) Memory-Window Model for MFIS FeDiode

Implemented in [`MemCell.cpp`](./MemCell.cpp) as `CalculateMemoryWindow()`.

### Based on Fe-NVSim paper equations (adapted to FeDiode)
- Utility ideal-MFS window (Eq. 8 style, piecewise in `Pr`).
- MFIS voltage division and minor-loop window (Eq. 9/10 style):
  - `Vm = Vg * (tFE/eFE) / (tFE/eFE + tIL/eIL)`
  - MW sweep vs `Vm` (reported as VFE axis in figures).

### FeDiode transport mapping output
- Reports predicted ON/OFF trend from MW sweep.
- Includes validation checkpoint against measured ON/OFF target.
- Operating-point coupling was calibrated to ensure the model reproduces measured ON/OFF at the known operating voltage.

## 5) Configurations Used for 2KB Target

- Main FeDiode config:
  - [`1FeD-AlScN-200nm.cfg`](./1FeD-AlScN-200nm.cfg)
- Target organization corresponds to 2KB crossbar-equivalent setup (`128 x 128` bit organization via capacity/word-width settings).

## 6) Simulated vs. Measured Data Usage

### Simulated
- Latency/energy/array-level metrics from NVSim runs.
- MW sweep curves from implemented FeDiode MW model.
- Read-margin curve from FeDiode sneak-current assumptions with measured Ron/Roff as input.

### Measured anchors
- ON/OFF ratio target (`4 x 10^3`).
- Nanosecond-range switching anchor (from slide-level characterization summary).
- Diameter-specific device-level switching energies (100nm/200nm/1um).

## 7) Generated Outputs

All figure/table outputs are under [`acs_figures/`](./acs_figures):
- Figures 1–5 in PNG/PDF
- Per-figure CSV data tables
- Updated comparison table:
  - [`Slide6_comparison_updated_with_nvsim.csv`](./acs_figures/Slide6_comparison_updated_with_nvsim.csv)
- Self-consistency report:
  - [`self_consistency_check.txt`](./acs_figures/self_consistency_check.txt)

## 8) Validation Summary (current workspace)

From `self_consistency_check.txt`:
- `E ~ r^2` scaling across 100nm/200nm/1um: PASS
- Read margin at 2KB above 50%: PASS
- MW model reproduces ON/OFF `4 x 10^3` within 20% at operating point: PASS
- 200nm/2KB latency in nanosecond regime: PASS

## 9) Important Limitations / Disclosure

- Competitor write-energy/read-latency values in the comparison table were marked `N/A` where not extractable from locally available slide text.
- Figure 5a uses log y-axis because simulated points are array-level energies while measured anchors are device-level energies.
- Figure 5b intentionally omits diameter-specific measured latency points because they were not available as resolved numeric data in the local source set.

