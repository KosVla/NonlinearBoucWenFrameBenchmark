# Frame with Bouc-Wen Hysteretic Links - Python Implementation

A Python implementation for simulating the dynamic response of a shear frame with nonlinear hysteretic connections under ground motion excitation. The framework uses Newmark implicit integration for structural dynamics coupled with explicit Runge-Kutta schemes for the Bouc-Wen hysteretic links.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Usage](#usage)
- [Code Architecture](#code-architecture)
  - [Core Modules](#core-modules)
  - [Input Files](#input-files)
- [Model Reference](#model-reference)
  - [Geometry and Connectivity](#geometry-and-connectivity)
  - [Material Class](#material-class)
  - [Excitation Class](#excitation-class)
  - [BoucWenModel Class](#boucwenmodel-class)
  - [ModelAssembly Class](#modelassembly-class)
  - [Newmark Class](#newmark-class)
- [Results Dictionary](#results-dictionary)
- [Plotting Utilities](#plotting-utilities)
- [Citation](#citation)
- [License](#license)

---

## Overview

This benchmark multi-degree-of-freedom simulator evaluates the numerical response of a **3D frame structure with hysteretic links** subjected to earthquake ground motion. The frame features:

- **Elastic beam elements** for structural members (Euler-Bernoulli formulation)
- **Nonlinear hysteretic links** modeled with the Bouc-Wen law including stiffness degradation and strength deterioration effects
- **Implicit Newmark integration** for structural dynamics (average acceleration method)
- **Explicit integration** (Euler/RK2/RK4) for hysteretic evolution
- **Modular object-oriented design** for easy extension and modification

---

## Requirements

```python
numpy >= 1.20.0
scipy >= 1.7.0
matplotlib >= 3.3.0
```

---

## Usage

### Quick Start

The main workflow is demonstrated in `Main_ExampleEarthquakes.py` and consists of:

1. **Defining the system & excitation parameters** via dedicated input file and classes
2. **Assembling all component to a dedicated object** for handling the system assembly
3. **Passing the Model to an integrator** 
4. **Running time integration**
5. **Extracting and visualizing results**

---

## Code Architecture

### Core Modules

Located in `core/` directory:

| Module | Description |
|--------|-------------|
| `Assembly.py` | `ModelAssembly` class - system matrix assembly, BC application, FE formulation |
| `Newmark.py` | `Newmark` class - implicit time integration with Newton-Raphson |
| `BoucWenModel.py` | `BoucWenModel` class - hysteretic model with explicit integration |
| `Excitation.py` | `Excitation` class - ground motion excitation |
| `Utils.py` | Plotting utilities for time histories, hysteresis, and phase space |

### Input Files

Located in `InputFiles/` directory:

| File | Content |
|------|---------|
| `InputFile.py` | Geometry definitions (`AssembleFrame_*` functions) and `Material` class |
| `*.txt` | Ground motion accelerograms (El Centro, Kobe, Morgan Hill) |

Input setups are located in `InputFiles/InputFile.py`. Two example configurations are provided:

| Function | Description | Hysteretic Link Locations |
|----------|-------------|---------------------------|
| `AssembleFrame_BeamsExample` | Simplified frame | Links in beams only |
| `AssembleFrame_FloorsExample` | Full two-story frame | Links in basement columns and beams |

The input populates the `MODEL` object with:
- `MODEL.Nodes`: Nodal coordinates (physical + virtual nodes)
- `MODEL.BeamElements`: Beam connectivity
- `MODEL.nl_link_elements`: Hysteretic link connectivity
- `MODEL.nl_link_flags`: Active DOFs for each link
- `MODEL.NodalDisplacements`: Boundary conditions
- `MODEL.ICs`: Initial conditions and damping
- `MODEL.n_dofs`: Total degrees of freedom

A more detailed presentation of these fields are given later. 
---

## Model Reference

### Geometry and Connectivity

#### Node Numbering Convention

Nodes are numbered using **1-based indexing** in the input files (for consistency with MATLAB), but are converted to **0-based indexing** internally for NumPy array operations.

**Virtual Node Strategy:**
- Physical structural nodes come first (e.g., nodes 1-18)
- Virtual nodes follow (e.g., nodes 19-52 for `BeamsExample`)
- Virtual nodes duplicate physical node coordinates to create connection points for hysteretic links
- Series connection: `virtual_node → [BW_link] → physical_node → [beam] → physical_node → [BW_link] → virtual_node`

#### `MODEL.Nodes`
**Type:** `numpy.ndarray`, shape `(n_nodes, 3)`

**Description:** Spatial coordinates `[x, y, z]` of all nodes.

**Example:**
```python
MODEL.Nodes = np.array([
    [0.0,   0.0,   0.0],   # Node 0 (1 in input): Origin
    [7.5,   0.0,   0.0],   # Node 1 (2 in input)
    [0.0,   0.0,   3.2],   # Node 2 (3 in input): First floor
    [0.0,   0.0,   0.0],   # Node 18 (19 in input): Virtual copy of node 0
])
```

#### `MODEL.BeamElements`
**Type:** `numpy.ndarray`, shape `(n_beams, 3)`

**Description:** Beam element connectivity.

**Columns:**
- Column 0: Element type index (references material properties)
- Column 1: Starting node (1-indexed in input)
- Column 2: Ending node (1-indexed in input)

**Example:**
```python
MODEL.BeamElements = np.array([
    [0, 1, 7],   # Type 0 beam from node 1 to node 7
    [0, 2, 8],   # Type 0 beam from node 2 to node 8
    [1, 7, 13],  # Type 1 beam from node 7 to node 13
])
```

#### `MODEL.nl_link_elements`
**Type:** `numpy.ndarray`, shape `(n_links, 3)`

**Description:** Hysteretic link connectivity.

**Columns:**
- Column 0: Starting node (1-indexed in input)
- Column 1: Ending node (1-indexed in input)
- Column 2: Link type index (references Bouc-Wen parameter sets)

**Example:**
```python
MODEL.nl_link_elements = np.array([
    [19, 1, 0],  # Link from virtual node 19 to physical node 1, type 0
    [20, 2, 0],  # Link from virtual node 20 to physical node 2, type 0
    [7, 26, 2],  # Link from physical node 7 to virtual node 26, type 2
])
```

#### `MODEL.nl_link_flags`
**Type:** `numpy.ndarray`, shape `(n_links, 6)`

**Description:** Binary flags indicating which DOFs have hysteretic behavior (1=active, 0=elastic).

**DOF Order:** `[ux, uy, uz, θx, θy, θz]`

**Example:**
```python
MODEL.nl_link_flags = np.ones((n_links, 6))  # All DOFs hysteretic
```

#### `MODEL.NodalDisplacements`
**Type:** `numpy.ndarray`, shape `(n_bc_nodes, 13)`

**Description:** Boundary condition specification.

**Columns:**
- Column 0: Node index (1-indexed in input)
- Columns 1-6: Constraint flags `[fx, fy, fz, θx, θy, θz]` (1=fixed, 0=free)
- Columns 7-12: Prescribed values `[ux, uy, uz, θx, θy, θz]`

**Example:**
```python
MODEL.NodalDisplacements = np.array([
    [19, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],  # Node 19: fully fixed
    [20, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # Node 20: fixed in x,z only
])
```

#### `MODEL.ICs`
**Type:** `dict`

**Description:** Initial conditions and damping parameters.

**Keys:**
- `'Disps'`: Initial displacement vector, shape `(n_dofs, 1)`
- `'Velocities'`: Initial velocity vector, shape `(n_dofs, 1)`
- `'Accelerations'`: Initial acceleration vector, shape `(n_dofs, 1)`
- `'zeta'`: List of modal damping ratios `[ζ1, ζ2]`
- `'OmegaIndexes'`: List of mode indices for Rayleigh damping `[mode1, mode2]`

**Example:**
```python
MODEL.ICs = {
    'Disps': np.zeros((ndim, 1)),
    'Velocities': np.zeros((ndim, 1)),
    'Accelerations': np.zeros((ndim, 1)),
    'zeta': [0.02, 0.02],  # 2% damping in modes 0 and 1
    'OmegaIndexes': [0, 1]
}
```

---

### Material Class

**Location:** `InputFiles/InputFile.py`

**Constructor:**
```python
Material(Eyoung, nee=0.3, rho=8000)
```

**Parameters:**
- `Eyoung`: Elastic modulus [Pa]. Can be scalar or list for multiple materials
- `nee`: Poisson's ratio (default: 0.3)
- `rho`: Mass density [kg/m³] (default: 8000)

**Attributes:**
- `mat_properties`: Array of shape `(n_materials, 3)` containing `[E, ν, ρ]` for each material
- `cross_sections`: Array of shape `(1, 4)` containing `[A, I1, I2, I3]` for beam cross-sections
  - `A`: Cross-sectional area [m²]
  - `I1`: Torsional constant [m⁴]
  - `I2`, `I3`: Moments of inertia [m⁴]

**Example:**
```python
# Single material (steel)
Mat = Material(210e9, nee=0.3, rho=8000)

# Multiple materials
Mat = Material([210e9, 210e9, 30e9, 30e9], nee=0.3, rho=8000)
```

---

### Excitation Class

**Location:** `core/Excitation.py`

**Constructor:**
```python
Excitation(fintegration=1000, angle=np.pi/4, Amp=10, type="Kobe")
```

**Parameters:**
- `fintegration`: Integration frequency [Hz] (default: 1000)
- `angle`: Ground motion direction angle [rad] from x-axis (default: π/4)
- `Amp`: Amplitude scaling factor (default: 10)
- `type`: Earthquake record type - `"Kobe"`, `"ElCentro"`, or `"Morgan"` (default: "Kobe")

**Attributes:**
- `dt`: Time step [s] = `1/fintegration`
- `angle`: Direction angle [rad]
- `Amp`: Amplitude scaling
- `SynthesizedAccelerogram`: Interpolated ground acceleration array, shape `(nt, 1)` [m/s²]
- `time`: Time vector, shape `(nt,)` [s]

**Ground Motion Files:**
- `InputFiles/KobeAccelNoScaling.txt`
- `InputFiles/ElCentroAccelNoScaling.txt`
- `InputFiles/MorganAccelNoScaling.txt`

Format: Two rows, comma-separated values
- Row 0: Time points [s]
- Row 1: Acceleration values [m/s²]

**Example:**
```python
# Kobe earthquake at 45° with 10× amplification
Forcing = Excitation(fintegration=1000, Amp=10.0, angle=np.pi/4, type="Kobe")
```

---

### BoucWenModel Class

**Location:** `core/BoucWenModel.py`

**Constructor:**
```python
BoucWenModel(bw_k=[2.0e8,2.0e8], bwa=0.25, Beta=3.0, Gamma=2.0,
             Alpha=1.0, N=1.0, deltav=0.0, deltan=0.0,
             integration_method='RK4')
```

**Bouc-Wen Parameters:**

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| `bw_a` | α | Elastic ratio (post-yield/initial stiffness)|
| `bw_k` | k | Link stiffness |
| `Alpha` | A | Hysteretic shape parameter |
| `Beta` | β | Hysteretic shape parameter |
| `Gamma` | γ | Hysteretic shape parameter |
| `N` | n | Smoothness exponent | 
| `deltav` | δ_ν | Strength deterioration rate |
| `deltan` | δ_η | Stiffness degradation rate |
| `integration_method` | - | Integration scheme: 'Euler', 'RK2', or 'RK4' |

**Multiple Material Types:**
- If `bw_k` is a list of length `n`, creates `n` different Bouc-Wen material types
- Other parameters are replicated for all types unless also provided as lists
- Link type is specified in `MODEL.nl_link_elements[:, 2]`

**Attributes:**
- `Models`: Array of shape `(n_types, 8)` containing parameters for each material type
- `Properties`: Array of shape `(n_active_dofs, 10)` assembled by `AssembleProperties()`
  - Columns 0-1: `[start_DOF, end_DOF]` global DOF indices
  - Columns 2-9: Bouc-Wen parameters `[bwa, bw_k, Alpha, Beta, Gamma, N, deltav, deltan]`
- `Hist`: Dictionary containing history variables for all active hysteretic DOFs
  - `'HistR'`: Restoring forces, shape `(n_active_dofs,)`
  - `'HistU'`: Relative displacements, shape `(n_active_dofs,)`
  - `'HistZeta'`: Hysteretic variable ζ, shape `(n_active_dofs,)`
  - `'HistE'`: Dissipated energy, shape `(n_active_dofs,)`
  - `'HystereticCurvesR'`: Force time history, shape `(n_active_dofs, nt)`
  - `'HystereticCurvesU'`: Displacement time history, shape `(n_active_dofs, nt)`

**Methods:**
- `AssembleProperties(nl_links)`: Assembles DOF-level properties from link-level definitions
- `evaluateBoucWen(U, BWhistory)`: Evaluates hysteretic forces and tangent stiffnesses
  - Returns: `(Rs, Ks, updated_history)` as tuple

**Example:**
```python
# Single material type with RK4 integration
BW = BoucWenModel(bw_k=2e8, bwa=0.25, integration_method='RK4')

# Four different material types
bw_ks = [1e8, 2e8, 1.5e8, 2.5e8]
BW = BoucWenModel(bw_k=bw_ks, bwa=0.25, integration_method='RK4')
```

---

### ModelAssembly Class

**Location:** `core/Assembly.py`

**Constructor:**
```python
ModelAssembly(Material, Excitation, BoucWenModel, perturb=np.ones((1,1)),
              InputFile='FloorsExample')
```

**Parameters:**
- `Material`: Material object defining properties and cross-sections
- `Excitation`: Excitation object defining ground motion
- `BoucWenModel`: BoucWenModel object defining hysteretic behavior
- `perturb`: Perturbation factors for beam stiffnesses, shape `(n_beams, 1)` (default: all ones)
- `InputFile`: Name of input configuration - `'FloorsExample'` or `'BeamsExample'` (default: 'FloorsExample')

**Key Attributes:**
- `Nodes`, `BeamElements`, `nl_link_elements`: Geometry (populated by input function)
- `K`: Global stiffness matrix (after BC application), shape `(n_dofs, n_dofs)`
- `M`: Global mass matrix (after BC application), shape `(n_dofs, n_dofs)`
- `C`: Rayleigh damping matrix, shape `(n_dofs, n_dofs)`
- `fint`: Internal force vector, shape `(n_dofs, 1)`
- `RHS`: Ground motion force matrix, shape `(n_dofs, nt)`
- `Displacements`, `Velocities`, `Accelerations`: State vectors, shape `(n_dofs, nt)`
- `BoucWenModel`: Reference to BoucWenModel object (contains `Hist` dictionary)

**Key Methods:**
- `initialize()`: Assembles system matrices, applies BCs, computes damping, assembles excitation
- `getSystemMatrices(DispsU, HistDict, Assemble)`: Assembles K, M, fint for current state
- `assembleBeams(DispsU, Kstiff, Mmass, fint)`: Assembles beam element contributions
- `applyBCs(Kstiff, Mmass, fint)`: Enforces boundary conditions by penalty method
- `AssembleSandZ(Mass)`: Creates ground motion loading vector using lumped masses
- `GetLumpedMass(Mass, translational_dofs)`: Computes lumped mass using HRZ scheme
- `RayleighDamping(Kstiff, Mmass, zetas, OmegasS)`: Computes α, β damping coefficients

**Example:**
```python
MODEL = ModelAssembly(Mat, Forcing, BW, perturb=np.ones((26,1)),
                      InputFile='BeamsExample')
# MODEL is now ready for time integration
```

---

### Newmark Class

**Location:** `core/Newmark.py`

**Constructor:**
```python
Newmark(Assembly, Parameters=None)
```

**Parameters:**
- `Assembly`: ModelAssembly object containing the structural system
- `Parameters`: Optional custom Newmark parameters

**Methods:**

#### `simulation(MODEL=None, Parameters=None)`

Performs implicit Newmark time integration with Newton-Raphson iteration.

**Algorithm:**
- Average acceleration method: `γ = 1/2`, `β = 1/4`
- Predictor-corrector scheme
- Newton-Raphson iteration for nonlinear equilibrium

**Returns:** Dictionary with keys:
- `'Displacements'`: Shape `(n_dofs, nt)` [m or rad]
- `'Velocities'`: Shape `(n_dofs, nt)` [m/s or rad/s]
- `'Accelerations'`: Shape `(n_dofs, nt)` [m/s² or rad/s²]
- `'HystereticR'`: Hysteretic forces, shape `(n_active_dofs, nt)` [N or N·m]
- `'HystereticU'`: Hysteretic displacements, shape `(n_active_dofs, nt)` [m or rad]
- `'OriginalDisps'`: Displacements of physical nodes only, shape `(n_physical_dofs, nt)`

**Example:**
```python
Solver = Newmark(MODEL)
Results = Solver.simulation()

# Extract displacement of DOF 36 over time
disp_36 = Results['Displacements'][36, :]
```

---

## Results Dictionary

The `Newmark.simulation()` method returns a dictionary with the following keys:

| Key | Shape | Description | Units |
|-----|-------|-------------|-------|
| `'Displacements'` | `(n_dofs, nt)` | Displacement time history (all nodes) | [m], [rad] |
| `'Velocities'` | `(n_dofs, nt)` | Velocity time history | [m/s], [rad/s] |
| `'Accelerations'` | `(n_dofs, nt)` | Acceleration time history | [m/s²], [rad/s²] |
| `'HystereticR'` | `(n_active_dofs, nt)` | Hysteretic restoring forces | [N], [N·m] |
| `'HystereticU'` | `(n_active_dofs, nt)` | Hysteretic relative displacements | [m], [rad] |
| `'OriginalDisps'` | `(n_physical_dofs, nt)` | Displacements of physical nodes only | [m], [rad] |

**DOF Indexing:**
Each node has 6 DOFs: `[ux, uy, uz, θx, θy, θz]`

Global DOF index for node `j` (0-indexed), local DOF `k` (0-5):
```python
dof_index = 6*j + k
```

**Example:**
```python
# Extract y-displacement of node 5 (0-indexed)
node_idx = 5
dof_y = 1  # y-direction translation
global_dof = 6*node_idx + dof_y
uy_node5 = Results['Displacements'][global_dof, :]
```

---

## Plotting Utilities

**Location:** `core/Utils.py`

### `Plot_TH_Response(Field, t, ndof_plot=None, title="Time History Response")`
Plots time history of a field quantity.

**Parameters:**
- `Field`: Array of shape `(n_dofs, nt)`
- `t`: Time vector, shape `(nt,)`
- `ndof_plot`: DOF to plot (default: DOF with maximum response)
- `title`: Plot title

### `Plot_Hysteretic_Curve(HistU, HistR, ndof_plot=None)`
Plots hysteretic force-displacement curve.

**Parameters:**
- `HistU`: Displacement array, shape `(n_active_dofs, nt)`
- `HistR`: Force array, shape `(n_active_dofs, nt)`
- `ndof_plot`: Link DOF to plot (default: DOF with maximum force)

### `Plot_Phase_Space(U, V, ndof_plot=None)`
Plots phase portrait (displacement vs. velocity).

**Parameters:**
- `U`: Displacement array, shape `(n_dofs, nt)`
- `V`: Velocity array, shape `(n_dofs, nt)`
- `ndof_plot`: DOF to plot (default: DOF with maximum velocity)

---

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

Apache-2.0 license

---

## Support

For questions, bug reports, or feature requests:

- **Issues:** Open an issue on GitHub with detailed description
- **Discussions:** Use GitHub Discussions for questions
- **Email:** vlachask@ibk.baug.ethz.ch

---

**Repository:** https://github.com/KosVla/NonlinearBoucWenFrameBenchmark

**Last Updated:** December 2025

**Version:** 1.0
