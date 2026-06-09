# BESFEM: Battery Electrode Simulation using MFEM (Multi-Particle)

BESFEM is a high-performance, MPI-enabled electrochemical simulation framework for lithium-ion battery electrodes.  
It supports **half-cell** and **full-cell** simulations, **Cahn–Hilliard** and **diffusion-based** transport models, and the **Smoothed Boundary Method (SBM)** for diffuse interfaces.  
The code is written in C++ and built on top of the **MFEM** finite-element library.

This branch supports the development of utilizing BESFEM for multi-particle applications. This is useful to understand the impacts of an electrode with blended chemistry. For example: a silicon-carbon anode. 


---

## Project Structure

```
BESFEM/
│
├── include/             # Header files for physics modules
├── src/                 # Source files for all simulations
├── inputs/
│   ├── mesh/            # Mesh files
│   ├── distance/        # SBM distance fields
│   ├── materials/       # Material Property Data
│   ├── Constants.cpp    # Constants
│   └── run_config.txt   # Simulation config file
│
├── outputs/
│   └── Results/         # Auto-generated simulation outputs
│
├── tests/               # Unit tests
├── plotting/            # Plotting files
└── bin/                 # Compiled executables

```

---

## Building BESFEM

Ensure MFEM, HYPRE, and MPI (OpenMPI or MPICH) are installed and available.

```bash
# clone the repository
git clone https://github.com/HCY-Group/BESFEM.git

# enter into besfem folder
cd BESFEM 

# On the HPCC 
module load MFEM

# compile all of the code - you may need to update  the makefile MFEM/HYPRE include + library paths as needed
make 

# enter folder with executable file
cd bin
```

---

## Running BESFEM

### Configuration File

BESFEM simulations are configured through a user-editable text file:

```text
inputs/run_config.txt
```

All simulation settings, including geometry, materials, initial lithiation, boundary conditions, and runtime options, are specified in this file.


### Example Configuration

```ini
mode = half
electrode = cathode

mesh_file = ../inputs/colored_labels_labels.tif
mesh_type = v

cathode_distance = ../inputs/dummy.gf
anode_distance = ../inputs/dummy.gf

num_steps = 1000
combine_particles = false

cathode_materials = NMC,NMC,NMC
anode_materials = Graphite,Graphite,Graphite

init_cathode_particles = 0.30,0.30,0.30
init_anode_particles = 0.20,0.15,0.10

init_CnE = 0.001
init_BvA = -0.01
init_BvC = 3.96
init_BvE = -0.10
```

---

### Configuration Options

#### Simulation Mode

| Parameter   | Options                                          |
| ----------- | ------------------------------------------------ |
| `mode`      | `half`, `full`                                   |
| `electrode` | `anode`, `cathode` (used only in half-cell mode) |

Example:

```ini
mode = half
electrode = cathode
```

---

#### Mesh Settings

| Parameter   | Description                                 |
| ----------- | ------------------------------------------- |
| `mesh_file` | Mesh or TIFF file                           |
| `mesh_type` | `ml` (MATLAB mesh) or `v` (voxel/TIFF mesh) |

Example:

```ini
mesh_file = ../inputs/colored_labels_labels.tif
mesh_type = v
```

---

#### Distance Functions

| Parameter          | Description            |
| ------------------ | ---------------------- |
| `cathode_distance` | Cathode distance field |
| `anode_distance`   | Anode distance field   |

For half-cell cathode simulations:

```ini
cathode_distance = ../inputs/dummy.gf
anode_distancce = ../inputs/distance/dsF_3x90_r.txt
```

---

#### Time Stepping

```ini
num_steps = 1000
```

---

#### Particle Grouping

If all particles should behave as a single particle:

```ini
combine_particles = true
```

Otherwise:

```ini
combine_particles = false
```

---

#### Material Assignment

Materials are assigned per particle group.

Example:

```ini
cathode_materials = LFP,LFP,NMC
```

This corresponds to:

| Particle | Material |
| -------- | -------- |
| 0        | LFP      |
| 1        | LFP      |
| 2        | NMC      |

Supported cathode materials:

* `LFP`
* `NMC`

Supported anode materials:

* `Graphite`

Example:

```ini
anode_materials = Graphite,Graphite,Graphite
```

---

#### Initial Particle Lithiation

Initial lithiation is specified per particle.

Example:

```ini
init_cathode_particles = 0.15,0.20,0.10
```

corresponds to:

| Particle | Initial Lithiation |
| -------- | ------------------ |
| 0        | 0.15               |
| 1        | 0.20               |
| 2        | 0.10               |

---

#### Electrolyte Initial Condition

```ini
init_CnE = 0.001
```

---

#### Boundary Conditions

Initial potentials are specified using:

```ini
init_BvA = -0.01
init_BvC = 3.96
init_BvE = -0.10
```

where:

* `init_BvA` = anode potential
* `init_BvC` = cathode potential
* `init_BvE` = electrolyte potential

---

### Running a Simulation

#### Default Configuration

If using the default configuration file:

```bash
mpirun -np 8 ./battery_simulation
```

---

#### Specifying a Configuration File

```bash
mpirun -np 8 ./battery_simulation \
    -cfg ../inputs/run_config.txt
```

---

#### Specifying an Output File

To save the outputs of the simulation, redirect to a separate file: 

```bash
mpirun -np 8 ./battery_simulation \
    -cfg ../inputs/run_config.txt > output_file.txt
```

To filter to just the timestepping data, use `grep` :
```bash
cat output_file.txt | grep -e timestep > timestep.txt
```

---

#### Example: Mixed LFP/NMC Cathode

```ini
mode = half
electrode = cathode

mesh_file = ../inputs/colored_labels_labels.tif
mesh_type = v

cathode_distance = ../inputs/dummy.gf

num_steps = 3200
combine_particles = false

cathode_materials = LFP,LFP,NMC
anode_materials = Graphite,Graphite,Graphite

init_cathode_particles = 0.15,0.20,0.10
init_anode_particles = 0.02, 0.02, 0.02

init_CnE = 0.001
init_BvA = -0.1
init_BvC = 3.40
init_BvE = -0.10
```

Run:

```bash
mpirun -np 8 ./battery_simulation \
    -cfg ../inputs/run_config.txt
```






---

## Generating Doxygen Documentation

```bash
module load Doxygen
module load Graphviz
doxygen Doxyfile
cd html
```
Preview doxygen locally:
```bash
python3 -m http.server 8001 --bind 127.0.0.1
```

---

## Plotting Using PyGLVis

```bash
pip install glvis
```

To plot, please reference the `pyglivs.ipynb` file within the `plotting` folder. 
You will need to adjust the input files of `mesh` and the `GridFunction x` that you are plotting. 

---

## Core BESFEM Equations

**Cathode concentration**
```math
\frac{\partial C_c}{\partial t}
=
\frac{1}{\psi_c}
\nabla \cdot \left( \psi_c D_c \nabla C_c \right)
-
\frac{|\nabla \psi_c|}{\psi_c} r_c
```

**Cathode potential**
```math
\nabla \cdot \left( \psi_c \kappa_c \nabla \phi_c \right)
-
|\nabla \psi_c| z_- F r_c
=
0
```

**Cathode reaction rate**
```math
r_c
=
k_c^f C_+
\exp\left(
\frac{-\alpha z_+ F [\phi]_e^c}{RT}
\right)
-
k_c^b C_c
\exp\left(
\frac{(1-\alpha) z_+ F [\phi]_e^c}{RT}
\right)
```

**Electrolyte concentration**
```math
\frac{\partial C_e}{\partial t}
=
\frac{1}{\psi_e}
\nabla \cdot \left( \psi_e D_e \nabla C_e \right)
+
\frac{|\nabla \psi_c|}{\psi_e}\frac{r_c t_-}{\nu_+}
+
\frac{|\nabla \psi_a|}{\psi_e}\frac{r_a t_-}{\nu_+}
```

**Electrolyte potential**
```math
\nabla \cdot
\left[
\psi_e (z_+ m_+ - z_- m_-) F C_e \nabla \phi_e
\right]
+
|\nabla \psi_c| \frac{r_c}{\nu_+}
+
|\nabla \psi_a| \frac{r_a}{\nu_+}
=
\nabla \cdot
\left[
\psi_e (D_- - D_+) \nabla C_e
\right]
```

**Anode reaction rate**
```math
r_a
=
k_a^f C_+
\exp\left(
\frac{-\alpha z_+ F [\phi]_e^a}{RT}
\right)
-
k_a^b C_a
\exp\left(
\frac{(1-\alpha) z_+ F [\phi]_e^a}{RT}
\right)
```

**Anode concentration**
```math
\frac{\partial C_a}{\partial t}
=
\frac{1}{\psi_a}
\nabla \cdot
\left[
\psi_a M_a
\nabla
\left(
\frac{\partial f_G}{\partial C_a}
-
\varepsilon \nabla^2 C_a
\right)
\right]
-
\frac{|\nabla \psi_a|}{\psi_a} r_a
```

**Anode potential**
```math
\nabla \cdot \left( \psi_a \kappa_a \nabla \phi_a \right)
-
|\nabla \psi_a| z_- F r_a
=
0
```

**Mass Matrix:** used for any term with a time derivative.
```bash
AddDomainIntegrator(new mfem::MassIntegrator());
```

**Stiffness Matrix:** used for PDE terms involving ∇.
```bash
AddDomainIntegrator(new mfem::DiffusionIntegrator());
```

**Linear System Assembly:** Produces `A * X = B`
```bash
FormLinearSystem(ess_tdof_list, x, b, A, X, B);
```


<!-- 
# Running BESFEM

## Configuration File

BESFEM simulations are configured through a user-editable text file:

```text
inputs/run_config.txt
```

All simulation settings, including geometry, materials, initial lithiation, boundary conditions, and runtime options, are specified in this file. -->
