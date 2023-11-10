# Decentralized Target Tracking

Source code, posters, and other material related to PhD thesis 

> Robin Forsling, *The Dark Side of Decentralized Target Tracking: Unknown Correlations and Communication Constraints*, Linköping University, Linköping, Sweden, Dec. 2023.

## Thesis Source Code

The thesis source code is written in Matlab and is contained in `src/`

Comments about the source code:

* Make sure to add `src/thesis_lib` to your Matlab path to have access to all required functionality

* Some of the simulations involve optimization problems solved using YALMIP and the Mosek solver with a valid license. These can be found here:
  
  * YALMIP: https://yalmip.github.io
  
  * Mosek: https://www.mosek.com

* All numerical evaluations and all relevant thesis examples are included

* It is possible to add new scenarios by modifying `src/thesis_lib/scenario/load_scen.m` 

## Thesis Figures

All figures included in the thesis are available as pdf files in `figure/` 

Note, some of the figures are under &copy; IEEE.

## Thesis Posters

Posters related to the thesis are provided in `poster/` 
