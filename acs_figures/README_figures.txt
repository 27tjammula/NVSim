ACS Nano figure package generated from extended NVSim runs.

Files:
- Figure1_latency_vs_subarray.(png|pdf)
- Figure2_write_energy_vs_voltage.(png|pdf)
- Figure3_onoff_vs_vfe.(png|pdf)
- Figure4_read_margin_vs_size.(png|pdf)
- Figure5_energy_latency_vs_diameter.(png|pdf)
- Figure6_array_to_device_energy_ratio.(png|pdf)
- Figure1_latency_vs_subarray.csv
- Figure2_write_energy_vs_voltage.csv
- Figure3_onoff_vs_vfe.csv
- Figure4_read_margin_vs_size.csv
- Figure5_energy_latency_vs_diameter.csv
- Figure6_array_to_device_energy_ratio.csv
- Slide6_comparison_updated_with_nvsim.csv
- self_consistency_check.txt

Labeling convention:
- SIMULATED = NVSim output (or direct equations from implemented FeDiode model in SubArray/MemCell).
- MEASURED = anchors from device characterization/PPT summary (ON/OFF=4e3, nanosecond switching range, diameter-specific characterized devices).

Suggested interpretation:
- Figure 1 is best treated as a qualitative NVSim latency trend versus subarray size, not as a strong measurement-backed validation figure, because the current source set does not provide matched experimental read and write speed data.
- Figure 2 is best treated as a simulation/theory trend figure with per-cell analytical anchors at 3 V. The simulated array-access energy and the characterized single-device energy are intentionally different quantities, so the plot should not be over-claimed as a direct simulation-to-measurement match.
- Figure 3 is the strongest theory-to-measurement bridge in the current set: the MW model and WKB transport trend both connect cleanly to the measured ON/OFF anchor.
- Figure 4 is an array-feasibility/supporting figure showing how read margin evolves from the single-device limit through 8×8 and up to the 2 KB target array.
- Figure 5 is useful as a scaling/supporting figure, but its energy panel mixes array-level simulated energy with device-level measured anchors and should be described accordingly.
- Figure 6 is a useful supporting figure when discussing Figure 2, because it makes the array/device energy gap explicit instead of leaving it implicit.

Draft caption language:
- Figure 1. NVSim-projected read and write latency as a function of subarray size for the FeDiode crossbar. The write latency remains in the nanosecond regime across the evaluated organizations; this panel is intended to illustrate array-scaling behavior rather than provide a direct experiment-to-simulation comparison of both read and write speed.
- Figure 2. NVSim-projected write energy per access versus write voltage for 100 nm, 200 nm, and 1 um FeDiode devices. Cross markers at 3 V indicate per-cell analytical energy anchors derived from the extracted device parameters. Because the simulated curves are array-level access energies whereas the anchors are device-level quantities, the panel is intended to emphasize trend consistency across voltage and device size.
- Figure 3. ON/OFF ratio as a function of ferroelectric-layer voltage from the memory-window model at multiple coercive fields, together with a WKB transport estimate and the measured ON/OFF anchor. The close agreement at the operating point provides the most direct connection between the FeDiode model stack and the experimentally observed rectification level.
- Figure 4. NVSim-based read-margin projection for the FeDiode crossbar under a worst-case sneak-path condition. The single-device limit, an 8×8 array, and the 128×128 (2 KB) target organization are highlighted, showing that the projected read margin remains above the 50% criterion throughout the evaluated range.
- Figure 5. Projected write-energy and write-latency scaling with device diameter. The energy panel compares NVSim array-level write energy with device-level characterized energy anchors, while the latency panel shows the comparatively weak diameter dependence of the projected write latency in the present model.
- Figure 6. Ratio of NVSim array-level write energy to per-cell analytical write energy as a function of write voltage. This panel makes explicit the overhead between device-level and array-level energy and helps explain why the two quantities should not be interpreted as a direct one-to-one match in Figure 2.
