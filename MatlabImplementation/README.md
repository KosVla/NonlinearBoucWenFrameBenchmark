# Frame with Bouc-Wen hysteretic links as a multi-degree of freedom nonlinear response simulator - MATLAB Implementation

A MATLAB implementation for simulating the dynamic response of a shear frame with nonlinear hysteretic connections under ground motion. The framework uses Newmark implicit integration for structural dynamics coupled with explicit Runge-Kutta schemes for the Bouc-Wen hysteretic links.

---

## Table of Contents

- [Overview](#overview)
- [Usage](#usage)
  - [Input File & MODEL Struct Reference](#input-file-structure)
  - [Geometry and Connectivity](#geometry-and-connectivity)
  - [Material Properties](#material-properties)
  - [Hysteretic Link Configuration](#hysteretic-link-configuration)
  - [Boundary Conditions and Loads](#boundary-conditions-and-loads)
  - [Dynamic Analysis Parameters](#dynamic-analysis-parameters)
  - [System Matrices](#system-matrices)
  - [Results and History Variables](#results-and-history-variables)
- [Ground Motion Excitation](#excitation-and-boundary-conditions)
- [Citation](#citation)

---

## Overview

This benchmark multi degree of freedom simulator evaluates the numerical response of a **3D frame structure with hysteretic links** subjected to earthquake ground motion. The frame features:

- **Elastic beam elements** for structural members
- **Nonlinear hysteretic links** modeled with the Bouc-Wen law including stiffness degradation and strength deterioration effects
- **Plate elements**
- **Implicit Newmark integration** for structural dynamics
- **Explicit integration** (Euler/RK2/RK4) for hysteretic evolution
- **Flexible input specification** for geometry, materials, and structural elements

---

## Usage

The main workflow consists of:

1. **Define the structural model** via a dedicated input file
2. **Set the Bouc-Wen parameters for the hysteretic links**
3. **Specify excitation** 
4. **Run time integration**
5. **Extract and visualize results**

### Input File & MODEL Struct Reference

Input files are located in the `InputFiles/` directory. Each file defines a `MODEL` struct that contains all fields and variables necessary for assembling the benchmark frame. The input file offers the possibility to assemble nonlinear hysteretic links on the frame and includes advanced possibilities for spring elements, plate elements, and lumped masses. The following template input files are provided, all referring to a two-story frame:
| File | Location of Hysteretic Links |
|------|------------------------------|
| `InputFileLinksBeams.m` | Links in beams, all follow the same BW formulation |
| `InputFileLinksFirst.m` | Links in beams + basement columns, all follow the same BW formulation |
| `InputFileLinksAll.m`   | Links in beams + basement & first-story columns, all follow the same BW formulation |
| `InputFileLinksFirst_DifferentBW.m` | Links in beams + basement columns, chosen links follow a different BW formulation |

The `MODEL` struct is the central data structure containing all information about the structural system, from geometry and material properties to dynamic analysis parameters and simulation results. This section provides comprehensive documentation of all fields.

### Geometry and Connectivity

#### `MODEL.nodes`
**Type:** Matrix ∈ ℝ^(n_nodes × 3)

**Description:** This matrix contains the spatial coordinates of all nodes in the structural model, where n_nodes is the total number of nodes including both real structural nodes and virtual nodes created for hysteretic link modeling. Each row represents one node, and the three columns contain the x, y, and z coordinates in the global coordinate system.

**Details:**
- When hysteretic links are present, nodes are duplicated to create virtual links that serve as connection points for the Bouc-Wen elements.
- For example, if a hysteretic link connects structural nodes 1 and 2, you must create virtual nodes 19 and 20 at the same physical locations as nodes 1 and 2. The connectivity becomes: virtual_node_19 → BW link → structural_node_1 → [beam] → structural_node_2 → BW link → virtual_node_20.
- The first entries are the structural nodes, and the remaining nodes are virtual copies for link modeling.

**Example:**
```matlab
MODEL.nodes = [
    0.0   0.0   0.0    % Node 1: Origin
    7.5   0.0   0.0    % Node 2: X-direction
    0.0   0.0   3.2    % Node 3: Z-direction (1st floor)
    0.0   0.0   0.0    % Node 19: Virtual copy of Node 1 (for link)
];
```

#### `MODEL.beam_elements`
**Type:** Matrix ∈ ℝ^(n_beams × 3)

**Description:** This matrix defines the connectivity of all beam elements in the structure, where n_beams is the total number of beam elements. Each row represents one beam element with three entries: the element type, the starting node number, and the ending node number.

**Details:**
- Column 1: Element type (One element type is supported, for additional ones the user has to specify the beam assembly functions)
- Column 2: Starting node index (referring to row number in `MODEL.nodes`)
- Column 3: Ending node index (referring to row number in `MODEL.nodes`)
- Beam elements connect the real structural nodes and do not include the virtual nodes used for hysteretic links.

**Example:**
```matlab
MODEL.beam_elements = [
    0   1   7     % Beam from node 1 to node 7 (basement column)
    0   2   8     % Beam from node 2 to node 8
    0   7   13    % Beam from node 7 to node 13 (1st floor column)
];
```

#### `MODEL.nl_link_elements`
**Type:** Matrix ∈ ℝ^(n_links × 2)

**Description:** This matrix defines the connectivity of nonlinear hysteretic link elements, where n_links is the total number of hysteretic connections. Each row represents one link with two entries: the starting node and the ending node. These nodes must be the virtual duplicate nodes created specifically for hysteretic link modeling.

**Details:**
- The hysteretic links connect virtual nodes to real structural nodes, creating a series connection: virtual_node → [hysteretic_link] → structural_node → [beam] → structural_node → [hysteretic_link] → virtual_node.
- Each link element has 6 degrees of freedom at each end (12 total), but only the active DOFs specified in `MODEL.nl_link_flags` will exhibit hysteretic behavior.

**Example:**
```matlab
MODEL.nl_link_elements = [
    19   1    % Link from virtual node 19 to structural node 1
    20   2    % Link from virtual node 20 to structural node 2
    7    26   % Link from structural node 7 to virtual node 26
];
```

#### `MODEL.ndim`
**Type:** Scalar ∈ ℝ

**Description:** The total number of degrees of freedom in the structural system, computed as 6 × n_nodes, where each node has 6 DOFs (3 translations + 3 rotations).

**Details:**
- This value includes DOFs for both real structural nodes and virtual nodes.
- Used to initialize displacement, velocity, and acceleration vectors.
- **Critical:** This field must be defined in input files to avoid runtime errors.

---

### Material Properties

#### `MODEL.material_properties`
**Type:** Matrix ∈ ℝ^(n_materials × 3)

**Description:** This matrix contains the material properties for all distinct materials used in the structure, where n_materials is the number of different materials. Each row represents one material with three properties: elastic modulus (E), Poisson's ratio (ν), and mass density (ρ).

**Column Definitions:**
- Column 1: Elastic modulus E [Pa]
- Column 2: Poisson's ratio ν [dimensionless]
- Column 3: Mass density ρ [kg/m³]

**Example:**
```matlab
MODEL.material_properties = [
    210e9   0.30   8000    % Material 1: Steel
    30e9    0.15   2400    % Material 2: Reinforced concrete
];
```

#### `MODEL.beam_material_properties`
**Type:** Vector ∈ ℝ^(n_beams × 1)

**Description:** This vector assigns a material to each beam element by indexing into the `MODEL.material_properties` matrix. The i-th entry contains the material index (row number in `MODEL.material_properties`) to be used for the i-th beam element.

**Example:**
```matlab
MODEL.beam_material_properties = [
    1    % Beam 1 uses material 1 (steel)
    1    % Beam 2 uses material 1 (steel)
    2    % Beam 3 uses material 2 (concrete)
];
```

#### `MODEL.cross_sections`
**Type:** Matrix ∈ ℝ^(n_sections × 4)

**Description:** This matrix defines the geometric properties of all distinct cross-sectional shapes used in the structure, where n_sections is the number of different cross-section types. Each row represents one cross-section with four properties.

**Column Definitions:**
- Column 1: Cross-sectional area A [m²]
- Column 2: Torsional constant I₁ [m⁴]
- Column 3: Moment of inertia I₂ [m⁴]
- Column 4: Moment of inertia I₃ [m⁴]

**Example:**
```matlab
MODEL.cross_sections = [
    5.383e-3   204.3e-9   36.92e-6   13.36e-6    % Section 1: HEA 200 steel
    0.18       0.00381    0.00162    0.00486     % Section 2: 0.3m × 0.6m concrete
];
```

#### `MODEL.beam_cross_sections`
**Type:** Matrix ∈ ℝ^(n_beams × 2)

**Description:** This matrix assigns cross-sections to each beam element by indexing into the `MODEL.cross_sections` matrix. Each row corresponds to one beam element, with two entries allowing for tapered beams where the cross-section changes from one end to the other.

**Column Definitions:**
- Column 1: Cross-section index at the start node (left end)
- Column 2: Cross-section index at the end node (right end)

**Details:**
- For uniform beams, both columns have the same value.
- For tapered beams, different indices allow linear interpolation of properties.

**Example:**
```matlab
MODEL.beam_cross_sections = [
    1   1    % Beam 1: Uniform with section 1
    1   1    % Beam 2: Uniform with section 1
    2   2    % Beam 3: Uniform with section 2
];
```

---

### Hysteretic Link Configuration

#### Bouc-Wen model

The Bouc-Wen model captures rate-independent hysteretic behavior. A detailed description is found in the attached pdf. Below are the parameters used for its formulation:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `bw_a` | α | Elastic ratio (post-yield/initial stiffness)|
| `bw_k` | k | Link stiffness |
| `Alpha` | A | Hysteretic shape parameter |
| `Beta` | β | Hysteretic shape parameter |
| `Gamma` | γ | Hysteretic shape parameter |
| `N` | n | Smoothness exponent | 
| `deltav` | δ_ν | Strength deterioration rate | ≥ 0 |
| `deltan` | δ_η | Stiffness degradation rate | ≥ 0 |

#### `MODEL.nl_link_flags`
**Type:** Matrix ∈ {0,1}^(n_links × 6)

**Description:** This binary matrix controls which degrees of freedom of each hysteretic link are activated (bw_a < 1.0). Each row corresponds to one link element (from `MODEL.nl_link_elements`), and each column corresponds to one of the 6 local DOFs.

**Column Definitions:**
- Column 1: Axial DOF (translation along link axis)
- Column 2: Shear DOF 1 (translation perpendicular to axis)
- Column 3: Shear DOF 2 (translation perpendicular to axis)
- Column 4: Torsional DOF (rotation about link axis)
- Column 5: Bending DOF 1 (rotation about perpendicular axis)
- Column 6: Bending DOF 2 (rotation about perpendicular axis)

**Values:**
- 1 = Hysteretic behavior activated
- 0 = Effectively sets bw_a = 1.0

**Example:**
```matlab
MODEL.nl_link_flags = [
    1   1   1   0   0   0    % Link 1: Hysteretic in translations only
    1   1   1   1   1   1    % Link 2: Fully hysteretic (all 6 DOFs)
    0   0   0   0   0   0    % Link 3: No hysteresis
];
```

#### `MODEL.nl_link_bw_properties`
**Type:** Vector ∈ ℝ^8

**Description:** This vector contains the default Bouc-Wen parameters that will be applied to all hysteretic links unless overridden by `MODEL.nl_links_alternate`. The parameters are set in the main script via `MODEL.BW` and then extracted into this vector during model setup. The notation is: [bw_a, bw_k, Alpha, Beta, Gamma, N, deltav, deltan]

**Example:**
```matlab
MODEL.nl_link_bw_properties = [0.25  2e8  1.0  3.0  2.0  1.0  0.0  0.0];
```

#### `MODEL.nl_links_alternate`
**Type:** Matrix ∈ ℝ^(n_alternate × 9)

**Description:** This matrix allows specification of custom Bouc-Wen parameters for individual links that differ from the default parameters in `MODEL.nl_link_bw_properties`. Each row defines alternate parameters for one specific link element, where n_alternate is the number of links with custom properties.

**Column Definitions:**
- Column 1: Link element index (row number in `MODEL.nl_link_elements`)
- Columns 2-9: Bouc-Wen parameters [bw_a, bw_k, Alpha, Beta, Gamma, N, deltav, deltan]

**Details:**
- Only links that require non-default parameters need to be listed.
- The parameters in columns 2-9 completely override the default parameters for the specified link.
- All 6 active DOFs of the link will use these parameters.

**Example:**
```matlab
MODEL.nl_links_alternate = [
    5   0.4   2e8   1.0   3.0   2.0   1.0   0.0   0.0    % Link 5 has α=0.4
    7   0.3   1.5e8  1.2   2.5   2.5   1.0   0.02  0.01   % Link 7 with degradation
];
```

#### `MODEL.BW`
**Type:** Struct

**Description:** This structure contains all Bouc-Wen model parameters and configuration options. It is defined in the main script and passed to the model setup functions. The parameters have been already describe. Only the additional fields are described here. 

**Fields:**

- **`MODEL.BW.bw_a`** 
- **`MODEL.BW.bw_k`** 
- **`MODEL.BW.Alpha`**
- **`MODEL.BW.Beta`**
- **`MODEL.BW.Gamma`**
- **`MODEL.BW.N`** 
- **`MODEL.BW.deltav`**
- **`MODEL.BW.deltan`**
- **`MODEL.BW.integration_method`** (String): Integration scheme for the Bouc-Wen evolution equation. Options are `'Euler'`, `'RK2'` (2nd order midpoint), or `'RK4'` (4th order).

#### `MODEL.BW.HystLinks`
**Type:** Matrix ∈ ℝ^(n_dofs_active × 10)

**Description:** This internal matrix is assembled automatically by `AssembleHystereticLinks()` and contains the DOF-level properties for all active hysteretic degrees of freedom, where n_dofs_active is the total number of activated hysteretic DOFs across all links.

**Column Definitions:**
- Columns 1-2: [start_DOF, end_DOF] - Global DOF indices
- Columns 3-10: Bouc-Wen parameters for this DOF

**Details:**
- Each row represents one active hysteretic DOF (not one link element).
- A link with 3 active DOFs contributes 3 rows to this matrix.
- This vectorized format enables efficient parallel evaluation of all hysteretic links.

#### `MODEL.BW.HistBW`
**Type:** Struct

**Description:** This structure contains the current state variables for all active hysteretic DOFs. It is updated at each time step during the simulation and represents the "memory" of the hysteretic system.

**Fields:**

- **`MODEL.BW.HistBW.R`** (Vector ∈ ℝ^(n_dofs_active × 1)): Current restoring forces at each active hysteretic DOF in [N] or [N·m].

- **`MODEL.BW.HistBW.Um`** (Vector ∈ ℝ^(n_dofs_active × 1)): Current relative displacements at each active hysteretic DOF in [m] or [rad].

- **`MODEL.BW.HistBW.Zeta`** (Vector ∈ ℝ^(n_dofs_active × 1)): Current values of the hysteretic variable z (dimensionless), representing the internal state of the Bouc-Wen model.

- **`MODEL.BW.HistBW.E`** (Vector ∈ ℝ^(n_dofs_active × 1)): Cumulative dissipated energy at each active hysteretic DOF in [J], used for degradation calculations.

---

### Boundary Conditions and Loads

#### `MODEL.nodal_displacements`
**Type:** Matrix ∈ ℝ^(n_bc × 13)

**Description:** This matrix specifies displacement boundary conditions (supports and prescribed displacements) for nodes in the structure, where n_bc is the number of nodes with at least one constrained DOF. Each row defines the boundary conditions for one node.

**Column Definitions:**
- Column 1: Node index
- Columns 2-7: Constraint flags [fx, fy, fz, frx, fry, frz] where 1=constrained, 0=free
- Columns 8-13: Prescribed displacement values [ux, uy, uz, θx, θy, θz]

**Details:**
- For fixed supports, set flags to 1 and values to 0.
- For prescribed non-zero displacements (e.g., support settlements), set flags to 1 and specify values.
- Boundary conditions are enforced using direct elimination (zeroing matrix rows/columns, setting diagonal to 1).
- For hysteretic base isolation, apply BCs to virtual ground nodes, not structural base nodes.

**Example:**
```matlab
MODEL.nodal_displacements = [
    19  1  1  1  1  1  1  0  0  0  0  0  0    % Node 19: Fully fixed
    20  1  0  1  0  0  0  0  0  0  0  0  0    % Node 20: Fixed in x,z; free in y
    21  1  1  0  0  0  0  0  0  0.01  0  0  0 % Node 21: Prescribed uz=0.01m
];
```

#### `MODEL.nodal_loads`
**Type:** Matrix ∈ ℝ^(n_loads × 7)

**Description:** This matrix specifies concentrated forces and moments applied directly at nodes, where n_loads is the number of loaded nodes. 

**Column Definitions:**
- Column 1: Node index
- Columns 2-7: Load components [Fx, Fy, Fz, Mx, My, Mz] in [N] and [N·m]

**Example:**
```matlab
MODEL.nodal_loads = [
    12   30e3   0   0   0   0   0     % Node 12: 30 kN in x-direction
    8    0   -50e3   0   0   0   0    % Node 8: 50 kN downward (y)
];
```

#### `MODEL.beam_loads`
**Type:** Matrix ∈ ℝ^(n_beams × 6)

**Description:** This matrix specifies distributed loads along beam elements, where each row corresponds to one beam element. The loads can vary linearly from the start to the end of the beam.

**Column Definitions:**
- Columns 1-3: Load at start node [px1, py1, pz1] in [N/m]
- Columns 4-6: Load at end node [px2, py2, pz2] in [N/m]

**Details:**
- Positive values follow the global coordinate directions.
- For uniform loads, set columns 1-3 equal to columns 4-6.
- Self-weight is included via mass density in material properties.

**Example:**
```matlab
MODEL.beam_loads = [
    0   0   -1e3   0   0   -1e3    % Beam 1: Uniform -1 kN/m in z
    0   0   -2e3   0   0   -1e3    % Beam 2: Varying load (2→1 kN/m)
];
```

#### `MODEL.springs`
**Type:** Matrix ∈ ℝ^(n_springs × 7)

**Description:** This matrix specifies additional elastic spring elements attached to nodes, where n_springs is the number of spring-supported nodes. Each row adds springs to the 6 DOFs of one node.

**Column Definitions:**
- Column 1: Node index
- Columns 2-7: Spring stiffnesses [kx, ky, kz, krx, kry, krz] in [N/m] or [N·m/rad]

**Details:**
- Spring stiffnesses are added directly to the global stiffness matrix.
- Set to zero for DOFs without springs.

**Example:**
```matlab
MODEL.springs = [
    1   1e6   1e6   1e7   0   0   0    % Node 1: Elastic foundation
];
```

#### `MODEL.masses`
**Type:** Matrix ∈ ℝ^(n_masses × 7)

**Description:** This matrix specifies additional lumped masses attached to nodes, where n_masses is the number of nodes with added masses. Each row adds masses to the 6 DOFs of one node.

**Column Definitions:**
- Column 1: Node index
- Columns 2-7: Added masses [mx, my, mz, Jx, Jy, Jz] in [kg] or [kg·m²]

**Details:**
- Added masses are incorporated into the diagonal of the global mass matrix.
- Useful for modeling equipment, cladding, or other non-structural mass.

**Example:**
```matlab
MODEL.masses = [
    13   5000   5000   5000   0   0   0    % Node 13: 5000 kg added mass
];
```

---

### Dynamic Analysis Parameters

#### `MODEL.dyn`
**Type:** Struct

**Description:** This structure contains all parameters controlling the dynamic time-history analysis.

**Fields:**

- **`MODEL.dyn.dt`** (Scalar ∈ ℝ⁺): Time step size in [s]. This is the temporal discretization used for both the Newmark integration of structural dynamics and the Runge-Kutta integration of the Bouc-Wen model.

- **`MODEL.dyn.nt`** (Scalar ∈ ℕ): Total number of time steps in the analysis. The simulation duration is `T_total = dt × (nt - 1)`. This value must match the length of the excitation signal.

- **`MODEL.dyn.a`** (Scalar ∈ ℝ⁺): Rayleigh damping coefficient α for mass-proportional damping in [1/s]. The damping matrix is computed as `C = α·M + β·K`. Can be specified directly or computed from modal damping ratios.

- **`MODEL.dyn.b`** (Scalar ∈ ℝ⁺): Rayleigh damping coefficient β for stiffness-proportional damping in [s]. Can be specified directly or computed from modal damping ratios.

**Damping Computation:** If `MODEL.zeta` and `MODEL.OmegaIndexes` are provided in the main script, the Rayleigh coefficients are computed automatically to achieve the specified modal damping ratios at the given mode indices.

---

### System Matrices

These matrices are generated internally during the simulation and represent the assembled finite element system.

#### `MODEL.K`
**Type:** Sparse Matrix ∈ ℝ^(ndim × ndim)

**Description:** Global stiffness matrix of the assembled structure. 

#### `MODEL.M`
**Type:** Sparse Matrix ∈ ℝ^(ndim × ndim)

**Description:** Global mass matrix of the assembled structure.

#### `MODEL.Mall`
**Type:** Sparse Matrix ∈ ℝ^(ndim × ndim)

**Description:** Unmodified global mass matrix before boundary condition application. 

#### `MODEL.C`
**Type:** Sparse Matrix ∈ ℝ^(ndim × ndim)

**Description:** Global damping matrix computed as a linear combination of mass and stiffness matrices following Rayleigh damping: `C = α·M + β·K`, where α and β are the damping coefficients.

#### `MODEL.fint`
**Type:** Vector ∈ ℝ^(ndim × 1)

**Description:** Internal force vector representing elastic restoring forces from beam elements and hysteretic restoring forces from active Bouc-Wen links at the current configuration.

#### `MODEL.u`
**Type:** Vector ∈ ℝ^(ndim × 1)

**Description:** Current displacement vector containing the displacements and rotations of all DOFs at the current time step. This represents the instantaneous configuration of the structure.

#### `MODEL.Rmatrix`
**Type:** Matrix ∈ ℝ^(ndim × nt)

**Description:** Earthquake excitation force time history for all DOFs. This matrix contains the inertia forces due to ground motion computed as `F_eq = -S·a_ground(t)`, where S is the lumped mass distribution matrix.

**Details:**
- Column i contains the force vector at time step i.
- Only non-zero at DOFs corresponding to real structural nodes (not virtual nodes).

#### `MODEL.bc_dofs`
**Type:** Vector ∈ ℕ^(n_constrained × 1)

**Description:** Indices of all constrained degrees of freedom based on `MODEL.nodal_displacements`. This vector lists the global DOF numbers that have prescribed displacements (boundary conditions).

#### `MODEL.freedofs`
**Type:** Vector ∈ ℕ^(n_free × 1)

**Description:** Indices of all free (unconstrained) degrees of freedom. 

### Results and History Variables

These fields store the simulation results and time histories of all system quantities.

#### `MODEL.U`
**Type:** Matrix ∈ ℝ^(ndim × nt)

**Description:** Displacement time history for all degrees of freedom in the model. Column i contains the displacement vector at time step i, including translations [m] and rotations [rad] for all nodes (real and virtual).

**Indexing:** To extract displacement of node j, DOF k at time step i:
```matlab
dof_index = 6*(j-1) + k;
displacement = MODEL.U(dof_index, i);
```

#### `MODEL.V`
**Type:** Matrix ∈ ℝ^(ndim × nt)

**Description:** Velocity time history for all degrees of freedom. Column i contains the velocity vector at time step i, in [m/s] for translations and [rad/s] for rotations.

#### `MODEL.A`
**Type:** Matrix ∈ ℝ^(ndim × nt)

**Description:** Acceleration time history for all degrees of freedom. Column i contains the acceleration vector at time step i, in [m/s²] for translations and [rad/s²] for rotations.

#### `MODEL.Uorig`
**Type:** Matrix ∈ ℝ^(ndim_real × nt)

**Description:** Displacement time history for real structural nodes only, excluding virtual nodes used for hysteretic link modeling, where ndim_real = 6 × n_real_nodes. 

**Computation:**
```matlab
n_real = size(MODEL.nodes, 1) - size(MODEL.nl_link_elements, 1);
MODEL.Uorig = MODEL.U(1:(6*n_real), :);
```

#### `MODEL.BW.HistR`
**Type:** Matrix ∈ ℝ^(n_dofs_active × nt)

**Description:** Time history of hysteretic restoring forces for all active hysteretic DOFs. Row i contains the force time history for the i-th active DOF in [N] or [N·m].

**Details:**
- Each active DOF (flag=1 in `MODEL.nl_link_flags`) contributes one row.
- Forces are in the local coordinate system of each link.

#### `MODEL.BW.HistU`
**Type:** Matrix ∈ ℝ^(n_dofs_active × nt)

**Description:** Time history of relative displacements for all active hysteretic DOFs. Row i contains the displacement time history for the i-th active DOF in [m] or [rad].

#### `MODEL.time`
**Type:** Scalar ∈ ℝ⁺

**Description:** Wall-clock time in seconds required to complete the time integration simulation. This measures computational performance and is output by the `tic`/`toc` commands.

---

### DOF Indexing Convention

Each node in the model has **6 degrees of freedom** in the global coordinate system. The DOF numbering follows a sequential pattern:

**Global DOF Formula:**
```
DOF_index = 6 × (node_number - 1) + local_DOF
```

**Local DOF Mapping:**

| Local DOF | Global Index | Description | Units |
|-----------|--------------|-------------|-------|
| 1 | `6(j-1)+1` | Translation in x-direction | [m] |
| 2 | `6(j-1)+2` | Translation in y-direction | [m] |
| 3 | `6(j-1)+3` | Translation in z-direction | [m] |
| 4 | `6(j-1)+4` | Rotation about x-axis (torsion) | [rad] |
| 5 | `6(j-1)+5` | Rotation about y-axis | [rad] |
| 6 | `6(j-1)+6` | Rotation about z-axis | [rad] |

**Example:** For node 5, the y-direction displacement DOF is: `6×(5-1)+2 = 26`

## Ground Motion Excitation

Ground motion is specified as an **acceleration time history** in the global coordinate system.

**Input Fields:**
```matlab
Input.SynthesizedAccelerogram  % Vector ∈ ℝ^(nt × 1) in [m/s²]
Input.angle                    % Scalar in [rad], angle with x-axis
```

**Important Notes:**
- Earthquake forces are applied only to real structural nodes, not virtual nodes.
- Lumped mass is computed from the unmodified consistent mass matrix (`MODEL.Mall`).

## Citation

If you use this code in your research, please cite:

```bibtex
@article{Vlachas2021,
  title={A local basis approximation approach for nonlinear parametric model order reduction},
  author={Vlachas, Konstantinos and Tatsis, Konstantinos and Agathos, Konstantinos and Brink, Adam R and Chatzi, Eleni},
  journal={Journal of Sound and Vibration},
  volume={502},
  pages={116055},
  year={2021},
  publisher={Elsevier},
  doi={10.1016/j.jsv.2021.116055}
}
```

---

## License

[Specify your license here - e.g., MIT, GPL, Apache 2.0, etc.]

---

## Support

For questions, bug reports, or feature requests:

- **Issues:** Open an issue on GitHub with detailed description
- **Discussions:** Use GitHub Discussions for questions
- **Email:** vlachas@ibk.baug.ethz.ch
---

**Repository:** 

**Last Updated:** December 2025

**Version:** 1.0
