# ModAEM User's Manual

**Version 1.9.0**

*Vic Kelson, Ph.D. PE — WHPA Inc., Bloomington, IN*

ModAEM is a free, open-source analytic element groundwater flow model written in Fortran 95.
It uses the superposition of analytic functions — wells, line sinks, dipoles, and area sinks —
to compute heads and discharges in 2-D aquifer domains without the need for a numerical grid.

## Contents

- [Introduction](introduction.md) — History and philosophy of ModAEM
- [Groundwater Modeling Concepts](groundwater_modeling.md) — Analytic element background
- [Script Files](script_files.md) — Input file format and directives
- [Aquifer Specification](aquifer_specification.md) — The `AQU` module
- [Element Modules](element_modules.md) — Wells, line sinks, barriers, area sinks
- [Processing Directives](processing_directives.md) — Solution control directives
- [Analysis Modules](analysis_modules.md) — Grids, traces, observations
- [Validation](validation.md) — Test cases and benchmark results
- [Tools](tools.md) — Utility programs
- [Mathematical Formulation](mathematical_formulation.md) — Appendix: theory
