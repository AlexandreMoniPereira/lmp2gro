# Case 1 Tutorial: MOF-5 Conversion and Simulation

This tutorial provides a step-by-step guide for converting and simulating MOF-5. The necessary files are located in the `case1_files` folder:

- This README guide.
- LAMMPS data files for MOF-5 (neutral and charged variants).
- Necessary Hydrogen topology (.itp) files
- GROMACS parameter files (.mdp) for energy minimization and NPT ensembles.

---

### CIF2LAMMPS Pre-processing Notes
The LAMMPS data files for MOF-5 were generated using `cif2lammps`. To force the Fourier angle style within the topology generation, the `UFF4MOF_construction.py` script was modified.

In the `angle_parameters` function (line 465), the original conditional:

`if theta0_j not in (90.0, 120.0, 180.0):`

was replaced with:

`if True:`

This modification ensures a superior fit with the GROMOS angle potential, enhancing the numerical stability of the simulation.

---

## Step 1 - Utilizing lmp2gro

Generate the GROMACS topology and structural files from the uncharged data file (`data.MOF_5_opt_2x2x2`). From the `lmp2gro` root directory, execute:

```bash
python3 lmp2gro.py case1_files/data.MOF_5_opt_2x2x2 -r MOF --folder case1
```

This command populates the `case1` directory with the following files:
`atomtypes.itp`, `conf.gro`, `conf.itp`, `ffbonded.itp`, and `topol.top`.

**Manual Topology Adjustment:** Open `topol.top` and insert a semicolon (`;`) at the start of line 10 to comment out the `case2` defaults. Then, navigate into the case1 directory.


## Step 2 - Energy Minimization (EM)

A standard GROMACS execution requires three components:

1. **Topology** (`.top`/`.itp` files).
2. **Structure/Coordinates** (`.gro` file).
3. **Simulation Parameters** (`.mdp` file).

Copy `em.mdp` from `case1_files` into your current directory. Generate the binary run input (`.tpr`) using:

```bash
gmx grompp -f em.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 3
```
*Note: The `-maxwarn 3` flag is used to bypass harmless warnings.*

Run the minimization using `mdrun`:

```bash
gmx mdrun -v -s em.tpr -deffnm em_out
```
*`-v` enables verbose output; `-deffnm` sets the prefix for all output files.*

## Step 3 - NPT Simulation (Equilibration & Production)

Copy `npt.mdp` from `case1_files`. This setup follows the 20 ns protocol (5 ns equilibration + 15 ns production) described in the original paper.

Assemble the `.tpr` file:
```bash
gmx grompp -f npt.mdp -c em_out.gro -p topol.top -o npt.tpr -maxwarn 3
```

Execute the MD production run:
```bash
gmx mdrun -v -s npt.tpr -deffnm npt_out
```

## Step 4 - Simulation with Guest Molecules (H2)

### 4.1 Initializing the Working Directory
Return to the `lmp2gro` main folder and create a new instance:

```bash
python3 lmp2gro.py case1_files/data.MOF_5_opt_2x2x2 -r MOF --folder case1.1
```

Enter the new case1.1 folder.

Copy `HYD.itp` and `HYD.gro` (hydrogen parameters generated via `obgmx`) from the `case1_files` folder.

**Update `topol.top`:**

1. Comment out the `case2` defaults line.
2. Add `#include "HYD.itp"` alongside the other include directives.
3. Append `HYD 10` to the `[ molecules ]` section to account for 10 hydrogen molecules.

**Update `atomtypes.itp`:**
Append the non-bonded parameters for hydrogen at the end of the file:
```text
H_            1    1.0079  0.0000     A        0.2571              0.1842
```

### 4.2 Molecular Insertion
Insert the hydrogen molecules into the framework coordinates:

```bash
gmx insert-molecules -f conf.gro -ci HYD.gro -nmol 10 -o conf.gro
```

### 4.3 Execution
Repeat **Step 2** and **Step 3** using the updated `conf.gro` file to perform the guest-host simulation.
