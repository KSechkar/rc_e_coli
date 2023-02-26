%% dream_fit.m
% fit parameters to experimental data using DREAM

% NOTE:
% due to the nature of the package, the simulator and simulation parameters
% are defined within the function file dream_model.m and not here

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% DEFINE starting parameter values

% All parameters are constrained to be positive. Initial values from Weisse et al. 2015 and Chure et al. 2022
params = {
    {'a_r/a_a', 1} % metabolic gene transcription rate
    {'nu_max', 6000} % max. tRNA aminoacylatio rate
    {'K_t', 80000} % MM constants for translation elongation and tRNA charging rates
    {'kcm', 0.3594/1000} % chloramphenicol binding rate constant
    };

% record all values in a row vector
theta0=zeros(1,size(params,1)); % initialise
for i=1:size(params,1)
    theta0(i)=params{i}{2};
end

% initial sampling distirbution covariances (parameters assumed independent)
% equal to (1/4 mean)^2, where the mean is the 'initial value'
covar0=zeros(size(params,1),size(params,1));
for i=1:size(params,1)
    covar0(i,i)=(0.25.*params{i}{2}).^2;
end

%% SET UP the DREAM sampler

% User-defined problem settings
DREAMPar.d = size(params,1); % dimension of the problem
DREAMPar.N = 10; % number of Markov chains sampled in parallel
DREAMPar.T = 2000; % number of generations
DREAMPar.lik = 2; % 'model' function outputs numerical steady-state predictions
DREAMPar.outlier='iqr'; % handling outliers - 

% Backing up the progress every 10 steps
DREAMPar.save='yes';
DREAMPar.steps=10;

% Initial sampling distribution
Par_info.prior = 'normal'; % Sample initial state of chains from prior distribution
Par_info.mu=theta0; % means of prior distribution
Par_info.cov=covar0; % covariances of prior distribution (parameters distributed independently)

% Boundaries of parameter value domains
Par_info.min = [params{1}{2}/100,... % lower bound of a_a/a_r
    params{2}{2}/100,... % ditto nu_max
    params{3}{2}/100,... % ditto K_t
    params{4}{2}/100]; % ditto k_cm
Par_info.max = [params{1}{2}*100,... % upper bound of a_a/a_r
    params{2}{2}*100,... % ditto nu_max
    params{3}{2}*100,... % ditto K_t
    params{4}{2}*100]; % ditto k_cm
Par_info.boundhandling = 'fold'; % handle samples out of bounds by imagaing the domain as a torus, upholding MCMC detailed balance

% do not output diagnostic
DREAMPar.diagnostics='no';

% Parallel computing to speed the simulation up
DREAMPar.parallel = 'yes'; 
DREAMPar.CPU=2;

DREAMPar.restart='no';

Meas_info.Y=[];

%% RUN the sampler
[chain,output,fx]=DREAM('dream_model',DREAMPar,Par_info,Meas_info);