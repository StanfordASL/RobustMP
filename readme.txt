# RobustMP: Robust Online Motion Planning 

Code to accompany: [Robust Online Motion Planning via Contraction Theory and Convex Optimization](http://asl.stanford.edu/wp-content/papercite-data/pdf/Singh.Majumdar.Slotine.Pavone.ICRA17.pdf)

## Getting Started 

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. All code is written in MATLAB.

### Prerequisites

There are five packages required: 
* [spotless](https://github.com/spot-toolbox/spotless) : parsing the CCM optimization problems (sum-of-squares programming), 
* [Mosek](https://www.mosek.com/) : SDP solver for the CCM synthesis problems, 
* [TOMLAB](https://tomopt.com/tomlab/) : parsing the trajectory optimization and CCM controller problems (non-convex), 
* [SNOPT, NPSOL](https://ccom.ucsd.edu/~optimizers/) : solvers for the non-convex problems,
* [MPT3.0](https://www.mpt3.org/) : used for auxiliary plotting.

spotless and MPT are freely available. Mosek is freely available on an academic license. TOMLAB must be purchased (a 21-day trial license is available on their website). Evaluation copies of SNOPT and NPSOL are included with the trial TOMLAB download. In the near future, we will also release open-source alternatives for TOMLAB. 

### Installing

Having installed the prerequisites (and adding them to the MATLAB path), download the repo, navigate to the repo in MATLAB, and execute the following in the command window:

```
run startup_nmpc.m
```

This will add all necessary sub-directories of the repo to the MATLAB path.

### Summary

There are four main examples: PVTOL (planar quadrotor), FLR (flexible link robot), Quadrotor (3D), and a synthetic Tube MPC problem. To familiarize yourself with TOMLAB, we also provide code for solving a toy problem (The Brachistochrone Curve) within the folder: "Toy_Brach."

The PVTOL features global trajectory optimization with the following options for tracking: (i) track the global trajectory, or (i) update reference trajectory using local re-planning (i.e., the receding-horizon algorithm in the paper) and track this updating trajectory. 

The FLR features global trajectory optimization and tracking of this fixed trajectory only. 

The Quadrotor features global trajectory optimization using (i) waypoint generation via geometric FMT*, and (ii) polynomial spline smoothing. This trajectory is then held fixed during tracking. 

TubeMPC features the traditional NMPC algorithm (no global optimization, only receding-horizon optimization).

The different planning + tracking combinations aim to demonstrate the modularity of incorporating robustness using tubes. 

## Typical Workflow

The basic workflow is (i) Synthesize CCM using sum-of-squares optimization (offline), (ii) Choose an online planning scheme (like the examples above), (iii) Simulate with disturbances and track with the CCM controller. Below is a (very) quick guide on the PVTOL example. Additional details are coming as I write them. 

### CCM Computation

CCM synthesis methodology is pretty much the same for all systems: (i) define state-space constraints, (ii) do the global optimization to figure out optimal lambda and condition number, and (iii) compute and save this CCM using anonymous functions. The actual spotless code for doing the SOS optimization follows the same pattern across all examples.

PVTOL/Metric: 
*compute_CCM_PVTOL.m: 
- run either global optimizer (sweeps through lambda, bisection search over condition number of W), or fix lambda and condition number bound and compute and save metric
*CCM_PVTOL_Opt.m: 
- actual spotless code for computing metric
- metric parametrized as: W(x) = \sum_{i} W_i*m_i(x), where m_i is the i^th monomial and W_i are constant matrices.
- metric saved as two anonymous functions. w_poly_fnc: computes all monomials; and W_eval: computes the sum above. Also saved are any monomial gradient functions (i.e., gradients of w_poly_fnc). 

### Planning

Details provided here for PVTOL example. For FLR, Quadrotor, and TubeMPC, direct analogues can be found in similarly named files "load_<system>_config.m" and "load_solvers.m." All examples use the backend optimizer routines (using TOMLAB) for global optimization, local re-planning, NMPC, and the CCM controller. The only special case is the Quadrotor which has its own trajectory optimization and CCM controller routines, contained within the Quadrotor/Sim folder. 

Main simulation file is always called "Main.m." This calls the config and solvers files, does the simulation, and calls the plotting routines. 

PVTOL/Sim:
*Main.m: 
- load_PVTOL_config: constants, dynamics definitions, CCM, bounds, sets up environment
- load_solvers: computes global motion plan, sets up and computes first iteration of local motion plan (if enabled), and sets up and initializes CCM controller

### Simulation

PVTOL/Sim:
*Main.m
- can enable or disable local re-planning using "use_mpc"
- zoh implementation of controller at rate "dt_sim"
- the rest is simulation and plotting





