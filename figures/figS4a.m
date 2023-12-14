%% figS4a.m

% WINNER-TAKES-ALL PHENOMENON
% Figure S4a

% Two bistable switches (self-activating genes) in the same cell exhibiting
% winner-takes-all behaviour, when the activation of one switch may prevent
% the activation of the other. Here, we show that the timescale of a SINGLE
% switch's activation is determined by the inducer level in the medium

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP the simulator with default circuit parameters in WTA simulations

sim=cell_simulator; % initialise the cell simulator
sim=sim.load_heterologous_and_external('two_switches','pulse_inducer'); % load the het. gene and ext. inp. modules

% inducer addition time
inducer_add_time=72;

% 'inducer levels'
sim.het.parameters('f1')=20; % inducer for switch 1
sim.het.parameters('f2')=20; % inducer for switch 2

% max transcription rates
sim.het.parameters('a_switch1')=3000;
sim.het.parameters('a_switch2')=3000;
% Switch 1:
sim.het.parameters('K_dna(switch1)-switch1f1')=5e3; % gene reg. Hill constant
sim.het.parameters('eta_dna(switch1)-switch1f1')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline1')=0.1; % baseline promoter activity
% Switch 2:
sim.het.parameters('K_dna(switch2)-switch2f2')=5e3; % gene reg. Hill constant
sim.het.parameters('eta_dna(switch2)-switch2f2')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline2')=0.1; % baseline promoter activity

% inducer protein binding
sim.het.parameters('K_switch1-f1')=1000; % dissociation constant for inducer-protein binding
sim.het.parameters('K_switch2-f2')=1000; % dissociation constant for inducer-protein binding

% DO NOT TOUCH!
sim.ext.input_func_parameters('inducer_base_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('pulse_value_prop')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('pulse_start_time')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('pulse_duration')=inducer_add_time;% disturbance flicks transcription reg. func. from 0 to 1 at t=30
   
% push amended parameter values
sim=sim.push_het();

% simulation parameters
sim.tf =  144; % simulation time frame
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % ODE simulation tolerances

%% GET TRAJECTORIES of a SINGLE bistable switch with different inducer levels

% define range of first switch inducer concentrations
f1_sweep = [4,16,20,32,40,80,160,800,8000];

% second bistable switch absent
sim.het.parameters('c_switch2')=0;
sim=sim.push_het();

% initialise trajectory storage
t_trajs={}; % time axes of simulation
pswitch1_trajs={}; % p_switch1 levels

% run simulations
for i=1:size(f1_sweep,2)
    sim.het.parameters('f1')=f1_sweep(i);
    sim=sim.push_het();
    sim = sim.simulate_model;
    x_het=sim.x(:,10:(9+2*sim.num_het));
    % record
    t_trajs{i}=sim.t;
    pswitch1_trajs{i}=x_het(:,3);
end

%% PLOT TRAJECTORIES

% create figure
Fig_timescale = figure('Position',[0 0 385 295]);
set(Fig_timescale, 'defaultAxesFontSize', 9)
set(Fig_timescale, 'defaultLineLineWidth', 1.25)

hold on

% define trajectory colours
alpha_range=linspace(0.3,1,size(f1_sweep,2)); % less transparent = higher inducer conc.
basic_traj_colour=[0.75 0.75 0]; % basic trajectory colour

% plot the trajectories
for i=size(f1_sweep,2):(-1):1 
    plot(t_trajs{i}-inducer_add_time+5,pswitch1_trajs{i},'-', ...
        'Color',[basic_traj_colour, alpha_range(i)], ...
        'DisplayName',['f_1=',num2str(f1_sweep(i)),' nM']);
end


legend('Location','northwest','FontName','Arial')

xlabel('time [h]','FontName','Arial')
ylabel('p_{s1}, conc. of switch 1 prot. [nM]','FontName','Arial')

xticks(0:5:30)
xlim([0 30])
ylim([0 5e5])

grid on
axis square
box on
hold off

