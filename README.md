# rc-ecoli
MATLAB code for modelling how synthetic gene circuits interact with the host cell's (_E. coli_) native genes and affect its growth rate, used in the manuscript 'A coarse-grained bacterial cell model for resource-aware analysis and design of synthetic gene circuits'. The folders are organised as follows:

## cell_model
Contains scripts implementing the model of an E.coli cell. These include:
- _cell_simulator.m_ - a Matlab class object enabling simulations of the host cell. The expression of synthetic circuits can be simulated by loading the 'heterologous gene' and 'external input' modules (see _het_modules_ and _ext_inputs_ folders). Note that the associated modules' parameters must be pushed into the main simulator's memory using the function _push_het()_ every time they are altered.
- _cell_params.m_ - provides default values of all parameters describing the host cell
- _cell_init_conds.m_ - provides the default initial conditions for simulating the host cell model
- _cell_formulae.m_ - contains a collection of formulae for rate and activation functions used by the cell simulator
- _get_steady.m_ - an auxiliary function that allows to determine the cell's steady state in given conditions by running _cell_simulator_ simulations.

## het_modules
Contains scripts defining class objects that allow to model the expression of different heterologous gene circuits. These include:
- _no_het.m_ - no heterologous gene being expressed. Due to being 'empty', this script can be copied and filled in to describe a synthetic gene circuit of interest
- _one_constit.m_ - one constitutive heterologous gene
- _two_constit.m_ - two constitutive heterologous genes
- _two_switches.m_ - two self-activating heterologous genes acting as bistable switches. If the genes impose a sufficiently high burden on the cell, winner-takes-all behaviour means that the activation of one switch prevent the activation of the other
- _pi_controller.m_ - a proportional-integral feedback controller for maintaining a constant extent of competition for ribosomes in the cell

## ext_modules
Contains scripts defining class objects that describe external inputs administered to the cells, such as chemical inducer concentrations or light stimuli. These include:
- _no_ext.m_ - no external signal. Due to being 'empty', this script can be copied and filled in to describe an externa input of interest
- _constant_inducer.m_ - constant concentration of a chemical inducer
- _step_inducer.m_ - a step increase in the concentration of a chemical inducer
- _pulse_inducer.m_ - time-limited step pulse above or below the baseline concentration of a chemical inducer

## figures
Running the corresponding scripts allows to reproduce the simulation results presented in the publication. Files not named _fig(...).m_ are auxiliary scripts required by some of the figure-generating scripts.

## param_fitting
Contains scripts allowing to fit the model's parameters to experimental data obtained by Scott et al. [^1] and processed as descroibed by Chure et al. [^2] using the DiffeRential Evolution Adaptive Metropolis (DREAM) algorithm [^3]. Note that these scripts require the [DREAM](https://faculty.sites.uci.edu/jasper/software/#eleven) package to work.
- _dream_fit.m_ - run the DREAM algorithm and record the outcome in a _.mat_ file
- _dream_model.m_ - defines the function for the model's log-likelihood given the experimental data for an input set of parameter values. This function is called while sampling the chains during DREAM fitting
- _dream_modelfun.m_ - defines the function for finding the sum of squared errors between the experimental measurements and model predictions for an input set of parameter values. This is required to calculate the log-likelihood
- _DREAM_fitting_outcome.mat_ - stores the sampled chains obtained as outcome of running the the algorithm for 20,000 steps in total (10 parallel chains tun for 1,000 steps each)
- _DREAM_postproc_amended.mat_ - script for postprocessing the fitting outcome (find the mode of the sampled distribution, plot the sampled chains, etc.). This script has been taken from the original DREAM package, but minor amendments have been made to it so as to improve the output plots' interpretability

## data
Experimental data used to fit the parameters and compare model predictions with real-life measurements, taken from [^1] and [^2]. The text file _annotation.txt_ explains the meaning of each dataset present.

---

REFERENCES

[^1] Matthew Scott, Carl W. Gunderson, Eduard M. Mateescu, Zhongge Zhang, and Terence Hwa. Interdependence of cell growth and gene expression: Origins and consequences. Science, 330(6007):1099–1102, 2010.

[^2] Griffin Chure and Jonas Cremer. An optimal regulation of fluxes dictates microbial growth in and out of steady-state. biorXiv, 2022.

[^3] Jasper Vrugt. Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation. Environmental Modelling & Software. 75:273-316, 2016.
