%% figS4e.m

% WINNER-TAKES-ALL PHENOMENON
% Figure S4e

% Two bistable switches (self-activating genes) in the same cell exhibiting
% winner-takes-all behaviour, when the activation of one switch may prevent
% the activation of the other. However, here we decrease the maximum 
% expression rates of the genes - eventually, resource competition exerted
% by the two genes on each other becomes insufficient to cause
% winner-takes-all behaviour.

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% DEFINE the range of inducer 1 concs to consider
f1_sweep=[4:1:40,50:50:800,1000:1000:8000];

%% DEFINE the range of promoter strengths and Hill constant values to consider
% which factors we decrease the original values by
fold_decreases=[1,1.25,2.5,5,10];

% maximum values
a_switch1_max=3000;
a_switch2_max=3000;
K1_max=5e3;
K2_max=5e3;

%% SET UP the simulator with all other parameters
sim=cell_simulator;

sim=sim.load_heterologous_and_external('two_switches','pulse_inducer'); % load the het. gene and ext. inp. modules

% inducer addition time
inducer_add_time=10;

% 'inducer levels'
sim.het.parameters('f2')=20; % inducer for switch 2

% Switch 1:
sim.het.parameters('eta_dna(switch1)-switch1f1')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline1')=0.05; % baseline promoter activity
% Switch 2:
sim.het.parameters('eta_dna(switch2)-switch2f2')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline2')=0.05; % baseline promoter activity

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

% simulate
sim.tf =  50;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);


%% FIND equilibria for different inducer 1 concs across different rescalings of promoter strengths and Hill constants

eqbs={}; % storing the system's equilibria here

for fd_index=1:size(fold_decreases,2)
    % initialise the storage array
    eqbs{fd_index}=zeros(2, size(f1_sweep,2));

    % rescale promoter strengths and hill constants
    fold_decrease=fold_decreases(fd_index);
    sim.het.parameters('a_switch1')=a_switch1_max/fold_decrease;
    sim.het.parameters('a_switch2')=a_switch1_max/fold_decrease;
    sim.het.parameters('K_dna(switch1)-switch1f1')=K1_max/fold_decrease;
    sim.het.parameters('K_dna(switch2)-switch2f2')=K2_max/fold_decrease;
    
    % run simulations, record equilibria
    for i=1:size(f1_sweep,2)
        sim.het.parameters('f1')=f1_sweep(i);
        sim=sim.push_het();
        sim = sim.simulate_model;
        x_het=sim.x(:,10:(9+2*sim.num_het));
        eqbs{fd_index}(1,i)=x_het(end,3);
        eqbs{fd_index}(2,i)=x_het(end,4);
    end 
end

%% PHASE PLANE PLOT

% create figure
Fig_pp = figure('Position',[0 0 385 295]);
set(Fig_pp, 'defaultAxesFontSize', 9)
set(Fig_pp, 'defaultLineLineWidth', 1.25)

hold on

colour_start=[1 0 0];
colour_end=[0 0 1];
for fd_index=1:size(fold_decreases,2)
    % define a colour for a given rescaling's plot
    colour=(colour_start*(size(fold_decreases,2)-fd_index)+colour_end*(fd_index-1))/(size(fold_decreases,2)-1);

    % plot equilibria for different inducer 1 concentrations
    plot(eqbs{fd_index}(1,:),eqbs{fd_index}(2,:),'-o','Color',colour,'LineWidth',1.5,'MarkerSize',4)
end

legend({'\upsilon=1','\upsilon=1.25','\upsilon=2.5','\upsilon=5','\upsilon=10'},'FontSize',9)

xlabel('p_{s1}, conc. of switch 1 prot. [nM]','FontName','Arial')
ylabel('p_{s2}, conc. of switch 2 prot. [nM]','FontName','Arial')

xlim([0 5e5])
ylim([0 5e5])
xticks(0:1e5:5e5)
yticks(0:1e5:5e5)

grid on
axis square
box on
hold off