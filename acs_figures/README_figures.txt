ACS Nano figure package generated from extended NVSim runs.

Files:
- Figure1_latency_vs_subarray.(png|pdf)
- Figure2_write_energy_vs_voltage.(png|pdf)
- Figure3_onoff_vs_vfe.(png|pdf)
- Figure4_read_margin_vs_size.(png|pdf)
- Figure5_energy_latency_vs_diameter.(png|pdf)
- Figure1_latency_vs_subarray.csv
- Figure2_write_energy_vs_voltage.csv
- Figure3_onoff_vs_vfe.csv
- Figure4_read_margin_vs_size.csv
- Figure5_energy_latency_vs_diameter.csv
- Slide6_comparison_updated_with_nvsim.csv
- self_consistency_check.txt

Labeling convention:
- SIMULATED = NVSim output (or direct equations from implemented FeDiode model in SubArray/MemCell).
- MEASURED = anchors from device characterization/PPT summary (ON/OFF=4e3, nanosecond switching range, diameter-specific characterized devices).

Caption note for Figure 5:
- Figure 5a uses a logarithmic y-axis because SIMULATED values are array-level write energies while MEASURED anchors are device-level switching energies; this avoids visually collapsing the measured points near zero and preserves trend comparability.
