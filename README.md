# Elastic Interaction of Pressurized Cavities in Hyperelastic Media

This repository contains the computational and semi-analytical codes associated with the paper:

> **Saeedi, A., & Kothari, M. (2025). Elastic Interaction of Pressurized Cavities in Hyperelastic Media: Attraction and Repulsion.** *Journal of Applied Mechanics, 92*(5), 051008. [DOI: 10.1115/1.4067855](https://doi.org/10.1115/1.4067855)

---

## Overview

This work investigates the elastic interaction of two pressurized cylindrical cavities embedded in an elastic medium. The key findings are:

- In **linear elastic** media, two pressurized cavities **always attract** each other, regardless of separation.
- In **nonlinear (hyperelastic)** media (neo-Hookean, Mooney–Rivlin, Arruda–Boyce), the interaction **transitions from attraction to repulsion** above a critical pressure-to-shear-modulus ratio (P/μ ≳ 1).
- A critical inter-cavity separation distance (η/R)_critical governs where the driving force changes sign.

The interaction is quantified through the **driving force**, defined as the negative derivative of total potential energy with respect to the center-to-center separation η:

$$\mathcal{F} = -\frac{d\Pi^{\eta}_{\text{eqm}}}{d\eta}$$

---

## Repository Structure

```
.
├── README.md
│
├── Semi-Analytical (MATLAB)
│   ├── Analytical_LinearElastic.m         # Main script: computes strain energy, stresses, and driving force for two pressurized holes in a linear elastic medium using bipolar coordinates
│   ├── stress_cc.m                        # Function: radial stress σ_χχ in bipolar coordinates
│   ├── stress_xx.m                        # Function: hoop stress σ_ξξ in bipolar coordinates
│   └── stress_cx.m                        # Function: shear stress σ_χξ in bipolar coordinates
│
├── FEM (Python / ABAQUS)
│   ├── Two_pressurized_holes_ABAQUS_script_2D.py   # Full pipeline for two-cavity simulations: model creation → job submission → post-processing
│   └── One_pressurized_hole_ABAQUS_script_2D.py    # Baseline single-cavity simulation for validation (Appendix B of paper)
│
└── Post-processing (MATLAB)
    └── Deformed_area_calculator_2D.m      # Computes deformed cavity area via Delaunay triangulation of extracted node coordinates
```

---

## Dependencies

### MATLAB scripts
- MATLAB R2020a or later (older versions likely work)
- No additional toolboxes required

### Python / ABAQUS scripts
- ABAQUS 2022 or compatible version (scripts use the ABAQUS Python scripting interface)
- NumPy (bundled with ABAQUS's Python environment)
- Scripts must be run from within the ABAQUS/CAE through 'File > Run Script'

---

## Usage

### Semi-Analytical Solution (Linear Elasticity)

1. Open MATLAB and navigate to the repository directory.
2. Ensure `stress_cc.m`, `stress_xx.m`, and `stress_cx.m` are in the same directory as `Analytical_LinearElastic.m`.
3. Set parameters at the top of `Analytical_LinearElastic.m`:
   - `R`: hole radius
   - `E`, `nu`: Young's modulus and Poisson's ratio
   - `P`: applied internal pressure
   - `case_num`: `1` for a single separation distance, `2` to sweep over a range of `c` values (and thus separation distances)
4. Run `Analytical_LinearElastic.m`. Output includes:
   - Stress plots (radial, hoop, shear) vs. angular position on a test circle
   - Non-dimensionalized strain energy and potential energy vs. η/R
   - Non-dimensionalized driving force vs. η/R

> **Note on notation:** The scripts use the bipolar coordinate notation from Ling (1948) and Davanas (1992), with `chi` (χ) as the radial-like coordinate and `xi` (ξ) as the angular-like coordinate. The parameter `c` is the focal distance of the bipolar system; `chi0` is the χ-value corresponding to the hole boundary. The edge-to-edge distance is `L = 2*(sqrt(c^2 + R^2) - R)`, and the center-to-center distance is `eta = L + 2R`.

### Finite Element Simulations (Two Cavities)

**Step 1: Generate input files and models**

Open `Two_pressurized_holes_ABAQUS_script_2D.py` and configure:
- `base_dir`: path where input/output files will be saved
- `P`: applied pressure
- `d_min`, `d_max`, `d_increment`: range of center-to-center distances to sweep
- `R`, `R_m`: hole radius and domain radius
- `material_model`: `'NH'` (neo-Hookean), `'LE'` (linear elastic), `'MR'` (Mooney–Rivlin), or `'AB'` (Arruda–Boyce)
- Material constants for the chosen model (see parameter block in script)

Section 1 of the script inside ABAQUS/CAE creates models and write `.inp` files.

**Step 2: Submit jobs**

Section 2 submits all `.inp` files and wait for completion. Jobs are named `P{pressure}_dR{eta_over_R}`.

**Step 3: Extract results**

Section 3 extracts from each `.odb` file:
- Total strain energy → written to `StrainEnergy.rpt`
- Deformed nodal coordinates of the left hole boundary → written to `deformed_coordinates-{index}.csv`

### Deformed Area Calculation

After running the FEM simulations, configure `Deformed_area_calculator_2D.m`:
- `base_dir`: path to the CSV files from Step 3
- `P`, `d_min`, `d_max`, `d_increment`, `R`: must match the FEM simulation parameters

Run the script to compute the deformed cavity area at each separation distance using Delaunay triangulation. The area values are stored in `DeformedAreas`. These are needed to compute the potential energy term −P·ΔA (Eq. 9 of the paper).

### Single Cavity Baseline

`One_pressurized_hole_ABAQUS_script_2D.py` follows the same three-step structure as the two-cavity script but simulates only a single pressurized hole. This is used to validate the computational framework against the analytical solution for a single cavity (Appendix B of the paper) and to obtain the reference potential energy PE₀.

---

## Coordinate System and Key Quantities

The semi-analytical solution uses **bipolar coordinates** (χ, ξ), where:
- Circles of constant χ correspond to the hole boundaries and concentric circles around them
- The hole boundary is located at χ = χ₀ = sinh⁻¹(c/R)
- The Jacobian of the transformation is J = c / (cosh(χ) − cos(ξ))

Key non-dimensionalizations used throughout:
| Quantity | Non-dimensionalization |
|---|---|
| Separation distance | η/R (center-to-center / hole radius) |
| Potential energy | PE / \|PE∞\| |
| Driving force | 𝓕 / μ |
| Pressure | P / μ |

where PE∞ = −πP²R²/μ is the potential energy of two non-interacting cavities.

---

## Citation

If you use these codes in your work, please cite:

```bibtex
@article{saeedi2025elastic,
  title   = {Elastic Interaction of Pressurized Cavities in Hyperelastic Media: Attraction and Repulsion},
  author  = {Saeedi, Ali and Kothari, Mrityunjay},
  journal = {Journal of Applied Mechanics},
  volume  = {92},
  number  = {5},
  pages   = {051008},
  year    = {2025},
  doi     = {10.1115/1.4067855}
}
```

---

## Acknowledgments

Support from NH BioMade under the National Science Foundation EPSCoR award #1757371 is acknowledged. Insightful discussions with Professor Haneesh Kesari (Brown University) are acknowledged.

---

## Contact

- Ali Saeedi — ali.saeedi@unh.edu  
- Mrityunjay Kothari — mrityunjay.kothari@unh.edu  
  Department of Mechanical Engineering, University of New Hampshire
