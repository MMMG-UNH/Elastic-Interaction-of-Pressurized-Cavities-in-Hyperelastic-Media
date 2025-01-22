Multiphysics-Mechanics-of-Materials-Group-at-UNH

This repository contains the code and data accompanying the research paper:
"Elastic Interaction of Pressurized Cavities in Hyperelastic Media: Attraction and Repulsion"
Ali Saeedi, Mrityunjay Kothari
Published in Journal of Applied Mechanics, 2025

Abstract
[This study computationally investigates the elastic interaction of two pressurized cylindrical cavities in a 2D hyperelastic medium. 
Unlike linear elasticity, where interactions are exclusively attractive, nonlinear material models (neo-Hookean, Mooney-Rivlin, Arruda-Boyce) exhibit both attraction and repulsion between the cavities. 
A critical pressure-shear modulus ratio governs the transition, offering a pathway to manipulate cavity configurations through material and loading parameters.
At low ratios, the interactions are always attractive; at higher ratios, both attractive and repulsive regimes exist depending on the separation between the cavities.
Effect of strain stiffening on these interactions are also analyzed. These insights bridge theoretical and applied mechanics, with implications for soft material design and subsurface engineering.]

Citation
If you find this code useful in your research or work, please consider citing our paper:
@article{Saeedi2025,
  title     = {Elastic Interaction of Pressurized Cavities in Hyperelastic Media: Attraction and Repulsion},
  author    = {Saeedi, A. and Cothari, M.},
  journal   = {Journal of Applied Mechanics},
  year      = {2025},
  volume    = {XX},
  number    = {X},
  pages     = {XX-XX},
  doi       = {10.xxxx/xxxxxx},
}

Usage
[You can run the Python scripts in Abaqus by navigating to File > Run Script. Before running the scripts, ensure you update the parameters at the beginning of the file to match your specific problem.
The MATLAB file, Deformed_area_calculator_2D.m, can only be used after completing the simulations associated with the Python scripts and exporting the deformed coordinates to Excel files.]

Acknowledgments
We appreciate your interest in this work and welcome any feedback or contributions.
