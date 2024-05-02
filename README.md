# umat_finite_viscoelasticity

This repository consists the work of my Mini Thesis titled _"Implementation of a Viscoelastic Model for Hydrogels in ABAQUS"_ for the M.Sc. Computer-Aided Mechanical Engineering course at RWTH Aachen University. The outcome of this work was a User Material Subroutine based on the Finite Viscoelasticity Theory with multiple relaxation mechanisms and 1-parameter Ogden Models for the strain energy function.

Contents:
- UMAT          : ABAQUS User Subroutine FORTRAN Files (UMAT) 
- optimization  : MATLAB Live Scripts and Functions for parameter optimization using Levenberg-Marquardt optimization
- report        : LaTeX source and figures used for generating reportAJC.pdf


Changelog:
- 09.02.24      : Update feedback from XF Cheng
  - 1. Fix incorrect coefficient in UMAT files (residual calculation, c.f. Eq. 2.154 & 2.155 -> 1/(3*4) instead of 1/6)
  - 2. Fix incorrect sign in Report Eq. 2.158 and correspondingly update UMAT files
- 02.05.24      : Update feedback from Shields602
  - 1. Delete unused computation of `BTOTV`, and `SPRIND` call in UMAT files. 
