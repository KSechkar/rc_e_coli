%% fig9.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 9

% Showcasing how the controller enforces modularity

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP the simulators for both cases

% 2 different setups - ideal AIF and realistic scenario
sim=cell_simulator;

sim.init_conditions('s')=0.5;

sim=sim.load_heterologous_and_external('pi_controller_inducible','constant_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_level')=1; % inducer level: just full expression of xtra1 gene

% integral controller parameters
sim.het.parameters('K_dna(anti)-sens')=7000; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(anti)-sens')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('K_dna(amp)-act')=700; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(amp)-act')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('kb_anti')=300; % atcuator-annihilator binding rate constant
sim.het.parameters('c_sens')=100;
sim.het.parameters('a_sens')=50; % sensor gene transcription rate
sim.het.parameters('a_anti')=800; % annigilator transcription rate
sim.het.parameters('a_act')=400; % actuator transcription rate

sim.het.parameters('a_amp')=4000; % integral controller amplifier transcription rate

% regulated gene - assessing modularity
sim.het.parameters('c_ta')=100;
sim.het.parameters('c_x')=100;
sim.het.parameters('a_ta')=50;
sim.het.parameters('a_x')=50;
sim.het.parameters('K_ta-f')=1000;
sim.het.parameters('K_dna(x)-taf')=700; % gene reg. Hill constant
sim.het.parameters('eta_dna(x)-taf')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline')=0.1; % baseline promoter activity
   
% push amended parameter values
sim=sim.push_het();

% simulation parameters
sim.tf =  72;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% DEFINE plasmid concs. to be tested

sim=sim.push_het();
inducer_concs=logspace(log10(1e-1),log10(1e6),50);

% disturbing gene parameters
cdist_value=500;
adist_value=500;

% nutrient quality change parameters
s_value=0.25;
default_s_value=sim.init_conditions('s'); % back up the default s value


%% RUN simulations - reference

% initialise the array in which the results are stored
p_xs = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim.het.parameters('f')=inducer_concs(i);
    sim = sim.push_het();

    % simulate!
    sim = sim.simulate_model;

    % record
    p_xs(i)=sim.x(end,9+sim.num_het+7);
end

%% RUN simulations - competing gene

% distrubing gene
sim.het.parameters('c_dist')=cdist_value;
sim.het.parameters('a_dist')=adist_value;
sim=sim.push_het();

% initialise the array in which the results are stored
p_xs_comp = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim.het.parameters('f')=inducer_concs(i);
    sim = sim.push_het();

    % simulate!
    sim = sim.simulate_model;

    % record
    p_xs_comp(i)=sim.x(end,9+sim.num_het+7);
end

%% RUN simulations - nutrient availability change

% backtrack distrubing gene addition
sim.het.parameters('c_dist')=0;
sim.het.parameters('c_dist')=0;
sim=sim.push_het();

% change the medium
sim.init_conditions('s')=s_value;

% initialise the array in which the results are stored
p_xs_nutr = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim.het.parameters('f')=inducer_concs(i);
    sim = sim.push_het();

    % simulate!
    sim = sim.simulate_model;

    % record
    p_xs_nutr(i)=sim.x(end,9+sim.num_het+7);
end

%% SET UP the simulations for open loop

sim_openloop=cell_simulator;

sim_openloop.init_conditions('s')=sim.init_conditions('s');

sim_openloop=sim_openloop.load_heterologous_and_external('pi_controller_inducible','constant_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim_openloop.ext.input_func_parameters('inducer_level')=sim.ext.input_func_parameters('inducer_level'); % inducer level: just full expression of disturbing gene

sim_openloop.het.parameters('a_dist')=sim.het.parameters('a_dist');
sim_openloop.het.parameters('c_dist')=sim.het.parameters('c_dist');

% sensor protein concentration - just as a burden stand-in
sim_openloop.het.parameters('c_sens')=sim.het.parameters('c_sens');
sim_openloop.het.parameters('a_sens')=sim.het.parameters('a_sens'); % sensor gene transcription rate

% no controller or output protein expression here
sim_openloop.het.parameters('c_x')=0; % gene copy number
sim_openloop.het.parameters('c_act')=0; % gene copy number
sim_openloop.het.parameters('c_anti')=0; % gene copy number
sim_openloop.het.parameters('c_amp')=0; % gene copy number

% regulated gene - assessing modularity
sim_openloop.het.parameters('c_ta')=sim.het.parameters('c_ta');
sim_openloop.het.parameters('c_x')=sim.het.parameters('c_x');
sim_openloop.het.parameters('a_ta')=sim.het.parameters('a_ta');
sim_openloop.het.parameters('a_x')=sim.het.parameters('a_x');
sim_openloop.het.parameters('K_ta-f')=sim.het.parameters('K_ta-f');
sim_openloop.het.parameters('K_dna(x)-taf')=sim.het.parameters('K_dna(x)-taf');
sim_openloop.het.parameters('eta_dna(x)-taf')=sim.het.parameters('eta_dna(x)-taf');
sim_openloop.het.parameters('baseline')=sim.het.parameters('baseline');
  
% push amended parameter values
sim_openloop=sim_openloop.push_het();

% simulation parameters
sim_openloop.tf =  sim_openloop.tf;
sim_openloop.opt = sim_openloop.opt;

%% RUN simulations - reference - open loop

sim_openloop.het.parameters('c_dist')=0;
sim_openloop.het.parameters('a_dist')=0;
sim_openloop.init_conditions('s')=default_s_value;

% initialise the array in which the results are stored
p_xs_openloop = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim_openloop.het.parameters('f')=inducer_concs(i);
    sim_openloop = sim_openloop.push_het();

    % simulate!
    sim_openloop = sim_openloop.simulate_model;

    % record
    p_xs_openloop(i)=sim_openloop.x(end,9+sim.num_het+7);
end

%% RUN simulations - competing gene - open loop

% distrubing gene
sim_openloop.het.parameters('c_dist')=cdist_value;
sim_openloop.het.parameters('a_dist')=adist_value;
sim_openloop=sim_openloop.push_het();

% initialise the array in which the results are stored
p_xs_comp_openloop = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim_openloop.het.parameters('f')=inducer_concs(i);
    sim_openloop = sim_openloop.push_het();

    % simulate!
    sim_openloop = sim_openloop.simulate_model;

    % record
    p_xs_comp_openloop(i)=sim_openloop.x(end,9+sim.num_het+7);
end

%% RUN simulations - nutrient availability change

% backtrack distrubing gene addition
sim_openloop.het.parameters('c_dist')=0;
sim_openloop.het.parameters('c_dist')=0;
sim_openloop=sim_openloop.push_het();

% change the medium
sim_openloop.init_conditions('s')=s_value;

% initialise the array in which the results are stored
p_xs_nutr_openloop = zeros(1,size(inducer_concs,2));

for i=1:size(inducer_concs,2)
    disp(['Testing f=',num2str(inducer_concs(i))])

    % set plasmid concentration
    sim_openloop.het.parameters('f')=inducer_concs(i);
    sim_openloop = sim_openloop.push_het();

    % simulate!
    sim_openloop = sim_openloop.simulate_model;

    % record
    p_xs_nutr_openloop(i)=sim_openloop.x(end,9+sim.num_het+7);
end

%% Main Figure a - open loop

Fa = figure('Position',[0 0 350 260]);
set(Fa, 'defaultAxesFontSize', 9)
set(Fa, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(inducer_concs,p_xs_openloop,'Color',[0.6350 0.0780 0.1840])
plot(inducer_concs,p_xs_comp_openloop,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')
plot(inducer_concs,p_xs_nutr_openloop,'Color',[0.6350 0.0780 0.1840],'LineStyle',':')

xlabel('f, inducer conc. [nM]','FontName','Arial');
ylabel({'p_x, output prot. conc. [nM]'},'FontName','Arial')
legend({'c_{dist}=0 nM, \sigma=0.5', 'c_{dist}=500 nM, \sigma=0.5', 'c_{dist}=0 nM, \sigma=0.25'},...
    'Location','northwest')

xticks([0.1, 10, 1000, 100000])
ylim([0 16000])
xlim([0.1 1e5])

set(gca, 'XScale', 'log')

grid 
box on
axis square
hold off

%% Main Figure b - closed loop

Fb = figure('Position',[0 0 350 260]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(inducer_concs,p_xs,'Color',[0 0.4470 0.7410])
plot(inducer_concs,p_xs_comp,'Color',[0 0.4470 0.7410],'LineStyle',':')
plot(inducer_concs,p_xs_nutr,'Color',[0 0.4470 0.7410],'LineStyle','--')

xlabel('f, inducer conc. [nM]','FontName','Arial');
ylabel({'p_x, output prot. conc. [nM]'},'FontName','Arial')
legend({'c_{dist}=0 nM, \sigma=0.5', 'c_{dist}=500 nM, \sigma=0.5', 'c_{dist}=0 nM, \sigma=0.25'},...
    'Location','northwest')

xticks([0.1, 10, 1000, 100000])
xlim([0.1 1e5])
ylim([0 10000])

set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')

grid 
box on
axis square
hold off

%% Main Figure c - open loop

Fb = figure('Position',[0 0 350 260]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(inducer_concs,p_xs_openloop./p_xs_openloop,'Color',[0.6350 0.0780 0.1840])
plot(inducer_concs,p_xs_comp_openloop./p_xs_openloop,'Color',[0.6350 0.0780 0.1840],'LineStyle',':')
plot(inducer_concs,p_xs_nutr_openloop./p_xs_openloop,'Color',[0.6350 0.0780 0.1840],'LineStyle','--')

xlabel('f, inducer conc. [nM]','FontName','Arial');
ylabel({'p_x:p_x^0, relative output prot. conc. [nM]'},'FontName','Arial')
legend({'c_{dist}=0 nM, \sigma=0.5', 'c_{dist}=500 nM, \sigma=0.5', 'c_{dist}=0 nM, \sigma=0.25'},...
    'Location','northwest')

xticks([0.1, 10, 1000, 100000])
xlim([0.1 1e5])
ylim([0.5 2])

set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')

grid 
box on
axis square
hold off

%% Main Figure d - closed loop

Fb = figure('Position',[0 0 350 260]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(inducer_concs,p_xs./p_xs,'Color',[0 0.4470 0.7410])
plot(inducer_concs,p_xs_comp./p_xs,'Color',[0 0.4470 0.7410],'LineStyle',':')
plot(inducer_concs,p_xs_nutr./p_xs,'Color',[0 0.4470 0.7410],'LineStyle','--')

xlabel('f, inducer conc. [nM]','FontName','Arial');
ylabel({'p_x:p_x^0, relative output prot. conc. [nM]'},'FontName','Arial')
legend({'c_{dist}=0 nM, \sigma=0.5', 'c_{dist}=500 nM, \sigma=0.5', 'c_{dist}=0 nM, \sigma=0.25'},...
    'Location','northwest')

xticks([0.1, 10, 1000, 100000])
xlim([0.1 1e5])
ylim([0.5 2])

set(gca, 'XScale', 'log')
set(gca,'XMinorTick','off')

grid 
box on
axis square
hold off