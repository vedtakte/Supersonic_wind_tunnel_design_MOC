# Supersonic Wind Tunnel Nozzle Design (MATLAB)

This repository contains a **MATLAB-based design and analysis of a supersonic wind tunnel nozzle**
using **compressible flow theory**. The project focuses on generating the nozzle contour required
to achieve a desired supersonic Mach number distribution in the test section.

---

## Objectives
- Design a converging–diverging nozzle for supersonic flow
- Compute nozzle area variation using isentropic relations
- Generate nozzle contour using analytical / numerical methods
- Visualize Mach number, pressure, and geometry

---

## Theory Background
The design is based on:
- Isentropic compressible flow relations
- Area–Mach number relationship
- Supersonic nozzle flow physics
- Method of Characteristics for contour generation

---

## Code Structure and Entry Point

The primary nozzle design workflow is implemented in:

- **`nozzle_design.m`**  
  Main script for designing the supersonic wind tunnel nozzle.
  This file computes the nozzle contour based on compressible flow
  relations and serves as the **entry point** for the project.

Supporting scripts include:

- `calculation_coordinates.m` – Computes characteristic line coordinates
- `moc_mesh_plot.m` – Plots Method of Characteristics mesh
- `new_kernel_mesh_plot.m` – Kernel-based mesh visualization
- `optimized_design.m` – Optimized nozzle geometry
- `new_diff_plate_size.m` – Diffuser / plate sizing calculations
- `hsa_table_data.m`, `table.m` – Tabulated flow and helper data
- `allinone.m` – Combined execution script (optional / legacy)

Users are recommended to start by running **`nozzle_design.m`**.

## How to Run

1. Open MATLAB
2. Navigate to the repository folder
3. Run the main nozzle design script:

```matlab
nozzle_design
