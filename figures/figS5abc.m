% figS5abc.m

% EMERGENT BISTABILITY OF A NON-COOPERATIVE SELF-ACTIVATOR
% Figure 5: a,b,c

% Despite promoting its own expression non-cooperatively, a T7 RNAP gene
% controller by a constitutive T7 RNAP promoter may exhibit bistability due
% to the effect of synthetic gene expression on the host cell's growth
% rate. Here, we simulate the system staring at various initial
% conditions for different T7 RNAP toxcicity values to see if bistability
% emerges

%% CLEAR all variables

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% DEFINE the gene concetration
c_t7=5;

%% INITIALISE the cell simulator, define the circuit

sim=sim.load_heterologous_and_external('t7_selfact','no_ext'); % load heterologous gene expression module
sim.het.parameters('c_t7')=0; % first, set to zero for t7-rnap-free steady state retrieval
sim.het.parameters('a_t7')=170;
sim.het.parameters('K_t7')=550;
sim.het.parameters('b_t7')=16.64;
sim.het.parameters('d_t7')=0.18;
sim.het.parameters('baseline_t7')=0.0016; % Liao et al. 2017
sim.het.parameters('t7_toxicity')=0.03; % Liao et al. 2017
sim=sim.push_het;
s=0.5;

% simulation parameters
sim.tf =  72;

%% GET the cell's steady state in a given medium without heterologous gene expression

sim=sim.simulate_model();
steady_state_no_t7=sim.x(end,:);

% SET the initial condition according to the steady state w/out t7-rnap
% mRNAs
sim.init_conditions('m_a')=steady_state_no_t7(1);
sim.init_conditions('m_r')=steady_state_no_t7(2);
% proteins
sim.init_conditions('p_a')=steady_state_no_t7(3);
sim.init_conditions('R')=steady_state_no_t7(4);
% tRNAs
sim.init_conditions('tc')=steady_state_no_t7(5);
sim.init_conditions('tu')=steady_state_no_t7(6);
% inactivated ribosomes
sim.init_conditions('Bcm')=steady_state_no_t7(7);
% nutrient quality and chloramphenicol
sim.init_conditions('s')=steady_state_no_t7(8);
sim.init_conditions('h')=steady_state_no_t7(9);

% heterologous
x_het=steady_state_no_t7(10 : (9+2*sim.num_het+sim.num_misc) );
for j=1:sim.num_het
    % mRNA
    sim.init_conditions(['m_',sim.het.names{j}])=x_het(j);
    % protein
    sim.init_conditions(['p_',sim.het.names{j}])=x_het(sim.num_het+j);
end
% miscellaneous
for j=1:sim.num_misc
    sim.init_conditions([sim.het.misc_names{j}])=x_het(2*sim.num_het+j);
end

%% DEFINE the range of initial m_t7 concentrations

points_in_range=6;
mt7_range=5*ones(1,points_in_range);
pt7_range = logspace(log10(1),log10(1e4),points_in_range);

%% DEFINE pt7 toxicities to consider

tox_range = [0, 0.005, 0.12];

%% SIMULATE the system for different toxicities and initial conditions

% initialise T7 RNAP concentration staorage
pt7_trajs={{}}; % RNAP concentrations
t_trajs={{}}; % simulation time axes

for i=1:size(tox_range,2)
    for j=1:size(mt7_range,2)
        % starting from a low-t7 initial condition
        sim.het.parameters('t7_toxicity')=tox_range(i);
        sim.het.init_conditions('m_t7')=mt7_range(j);
        sim.het.init_conditions('p_t7')=pt7_range(j);
        sim.init_conditions('p_a')=(sim.init_conditions('p_a')*sim.parameters('n_a')-...
            sim.het.init_conditions('p_t7').*sim.het.parameters('n_t7'))./sim.parameters('n_a'); % adjust initial condition to have constant cell mass
        sim.het.parameters('c_t7')=c_t7;
        sim=sim.push_het;
        sim.init_conditions('s')=s;
        sim=sim.simulate_model();
        % record simulation outcomes
        pt7_trajs{i}{j}=sim.x(:,11);
        t_trajs{i}{j}=sim.t;
    end
end

%% FIGURE a - no toxicity

Fa = figure('Position',[0 0 215 193]);
set(Fa, 'defaultAxesFontSize', 9)
set(Fa, 'defaultLineLineWidth', 1.5)

hold on

% define trajectory colours
alpha_range=linspace(0.6,1,size(mt7_range,2)); % less transparent = higher inducer conc.
basic_traj_colour=[0 0.4470 0.7410]; % basic trajectory colour

% plot the trajectories
for j=size(mt7_range,2):(-1):1 
    plot(t_trajs{1}{j},pt7_trajs{1}{j},...
        '-','Color',[basic_traj_colour, alpha_range(j)], ...
        'DisplayName',['p_{t7}(0 h)=',num2str(pt7_range(j)),' nM']);
end

title(['\gamma_{t7}=',num2str(tox_range(1)),' nM'])
xlabel('t, time [h]','FontName','Arial');
ylabel('p_{t7}, T7 RNAP conc. [nM]','FontName','Arial')

set(gca, 'YScale', 'log')
ylim([1 1e4])
xlim([0 24])

% legend

grid 
box on
axis square
hold off


%% FIGURE b - low toxicity

Fb = figure('Position',[0 0 215 193]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.5)

hold on

% define trajectory colours
alpha_range=linspace(0.6,1,size(mt7_range,2)); % less transparent = higher inducer conc.
basic_traj_colour=[0 0.4470 0.7410]; % basic trajectory colour

% plot the trajectories
for j=size(mt7_range,2):(-1):1 
    plot(t_trajs{2}{j},pt7_trajs{2}{j},...
        '-','Color',[basic_traj_colour, alpha_range(j)], ...
        'DisplayName',['p_{t7}(0 h)=',num2str(pt7_range(j)),' nM']);
end

title(['\gamma_{t7}=',num2str(tox_range(2)),' nM'])
xlabel('t, time [h]','FontName','Arial');
ylabel('p_{t7}, T7 RNAP conc. [nM]','FontName','Arial')

set(gca, 'YScale', 'log')
ylim([1 1e4])
xlim([0 72])

% legend

grid 
box on
axis square
hold off

%% FIGURE c - high toxicity

Fc = figure('Position',[0 0 215 193]);
set(Fc, 'defaultAxesFontSize', 9)
set(Fc, 'defaultLineLineWidth', 1.5)

hold on

% define trajectory colours
alpha_range=linspace(0.3,1,size(mt7_range,2)); % less transparent = higher inducer conc.
basic_traj_colour=[0.8500 0.3250 0.0980]; % basic trajectory colour

% plot the trajectories
for j=size(mt7_range,2):(-1):1 
    plot(t_trajs{3}{j},pt7_trajs{3}{j},...
        '-','Color',[basic_traj_colour, alpha_range(j)], ...
        'DisplayName',['p_{t7}(0 h)=',num2str(pt7_range(j)),' nM']);
end

title(['\gamma_{t7}=',num2str(tox_range(3)),' nM'])
xlabel('t, time [h]','FontName','Arial');
ylabel('p_{t7}, T7 RNAP conc. [nM]','FontName','Arial')

set(gca, 'YScale', 'log')
ylim([1 1e4])
xlim([0 48])

% legend

grid 
box on
axis square
hold off