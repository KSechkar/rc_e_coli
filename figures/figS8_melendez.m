%% figS8.m
% SELF-ACTIVATING EGENE OSCILLATION SIMULATION
% Figure S8

% Compare model predictions for the dynamics of nutrient upshifts with
% experimental data collected by Erickson et al. and processed by Chure et
% al.
%% CLEAR

addpath(genpath('..'))

clear
close all

%%

% OD - cel conc. conversion
OD_to_ccells=(8e8*1000)/(6.02e23)*1e9; % OD of bacteria to concnentation of cells in nM

%% EXPERIMENTAL DATA

load exp_selfact_oscill.mat

LB_of_choice=100;

for i=1:length(LB)
    if(LB(i)==LB_of_choice)
        data.ts=time;
        data.ODs=OD_Lara(:,i);
        data.pswitchs=GFP_Lara(:,i);
        ODmax_new=max(data.ODs);
    end
end

%% SET UP the simulator

sim=cell_simulator; % initialise simulator
sim.opt = odeset('reltol',1.e-12,'abstol',1.e-12); % more lenient integration tolerances for speed
sim.init_conditions('s')=0.05; % sigma for LB broth, enabling the doubling time of about 20 min
sim=sim.load_heterologous_and_external('one_switch','no_ext');

% self-activating gene parameters
sim.het.parameters('a_switch1')=6000;
% Switch 1:
sim.het.parameters('K_dna(switch1)-switch1f1')=5e3; % gene reg. Hill constant
sim.het.parameters('eta_dna(switch1)-switch1f1')=2; % gene reg. Hill coefficient
sim.het.parameters('baseline1')=0.01; % baseline promoter activity
sim.het.parameters('K_switch1-f1')=1000; % dissociation constant for inducer-protein binding
sim=sim.push_het;


%% STEP 1: get to steady state in absence of induction

% no induction
sim.het.parameters('f1')=0; % absence of induction
sim=sim.push_het;

% getting steady state ( to find sigma values for different media)
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 6; % maximum no. iterations (checking if SS reached over first 750 h)

ss_noind = get_steady(sim,Delta,Max_iter);

%% STEP 2: 16h in batch culture in 100% LB with induction

% set up
sim.parameters('is_depletion')=1;
sim.parameters('ccells_max')=1.451*OD_to_ccells;
sim.het.parameters('f1')=10; % excess of inducer
sim.init_conditions=x2init_conds(sim,ss_noind);
sim.het.init_conditions('c_cells')=0.2*OD_to_ccells; % starting with OD 1/1000
sim.init_conditions('s')=0.25;
sim.parameters('K_nutr')=5e3;
sim=sim.push_het();
sim.tf=32;

% simulate
sim=sim.simulate_model;
sim.plot_simulation('heterologous','regulation')
% sim.plot_simulation('native','concentrations')

%%
figure()
plot(sim.t,1-sim.x(:,end-1)./sim.parameters('ccells_max'))
%% STEP 3: dilute 100-fold and put in media with different nutrient concentrations

% % set up
% sim.parameters('is_depletion')=0;
% % sim.het.parameters('f1')=17;
% sim=sim.push_het;
% sim.init_conditions=x2init_conds(sim,sim.x(end,:));
% sim.init_conditions('c_cells')=data.ODs(1)*OD_to_ccells; % 100-fold dilution
% sim.parameters('ccells_max')=ODmax_new*OD_to_ccells;
% sim.init_conditions('s')=0.25;
% 
% sim.tf=9;
% sim=sim.simulate_model;
% sim.plot_simulation('heterologous','regulation')
% % sim.plot_simulation('native','concentrations')
%% FIGURE S8a - cell OD

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

% experimentaldata
plot(data.ts, data.ODs, '.')

% model predictions
plot(sim.t,sim.x(:,end-1)./OD_to_ccells, ...
    '-','LineWidth',1.5)

xlabel('Time since upshift [h]')
ylabel('OD of culture');

xlim([0 sim.tf])
grid on
% set(gca, 'YScale', 'log')
axis square
box on
hold off

%% FIGURE S8b - self-activating protein expression

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

% experimental data
plot(data.ts,data.pswitchs./data.pswitchs(1),'.')

% model predictions
plot(sim.t,sim.x(:,end-2)./sim.x(1,end-2), ...
    '-','LineWidth',1.5)

xlabel('Time since upshift [h]')
ylabel('Het. prot. conc., nM');

xlim([0 sim.tf])
grid on
% set(gca, 'YScale', 'log')
axis square
box on
hold off

%% FIGURE S8c - cell growth rate

% get cell growth rates from model
ls=zeros(size(sim.t));
for i=1:length(sim.t)
    ls(i)=get_l(sim,sim.x(i,:));
end

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

% experimental data
plot(data.ts(1:end-1),(data.ODs(2:end)-data.ODs(1:end-1))./data.ODs(1:end-1),'.')

% model predictions
plot(sim.t,ls,'-','LineWidth',1.5)

xlabel('Time since upshift [h]')
ylabel('Growth rate [1/h]');

xlim([0 sim.tf])
grid on
% set(gca, 'YScale', 'log')
axis square
box on
hold off

%% FUNCTION for getting growth rate from the system's state vector x
function l=get_l(sim,x)
    % evaluate growth
    par=sim.parameters;
    m_a = x(1); % metabolic gene mRNA
    m_r = x(2); % ribosomal mRNA
    p_a = x(3); % metabolic proteins
    R = x(4); % operational ribosomes
    tc = x(5); % charged tRNAs
    tu = x(6); % uncharged tRNAs
    Bm = x(7); % inactivated ribosomes
    s=x(8); % nutrient quality
    h=x(9); % chloramphenicol conc.

    e=sim.form.e(par,tc); % translation elongation ratio

    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    % (no heterologous genes in this scenario)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    l=sim.form.l(par,e,B); % growth rate!
end

%% FUNCTION for getting the position of the element closest to reference
function closest_i=closest_pos(val_array,ref_val)
    distance_to_ref=sum(val_array)*100; % initialise the distance to reference growth rate with unreasonably high number
    for i=1:size(val_array,2)
        if(abs(val_array(i)-ref_val)<abs(distance_to_ref)) % if the current distance is the smallest so far, this is the neww reference
            closest_i=i;
            distance_to_ref=val_array(i)-ref_val;
        end
    end
end

%% FUNCTION for converting a state vector x into the initial condition disctionary
function init_conds=x2init_conds(sim,x)
    init_conds = containers.Map('KeyType', 'char', ...
            'ValueType', 'double');
    
    % mRNA concentrations - non-zero to avoid being stuck at lambda=0
    init_conds('m_a')=x(1); % metabolic
    init_conds('m_r')=x(2); % ribosomal
    
    % protein concentrations - start with 50/50 a/R allocation as a convention
    init_conds('p_a')=x(3); % metabolic *
    init_conds('R')=x(4); % ribosomal *

    % tRNA concentrations - 3E-5 abundance units in Chure and Cremer 2022 are equivalent to 80 uM = 80000 nM
    init_conds('tc')=x(5); % charged tRNAs
    init_conds('tu')=x(6); % free tRNAs

    % concentration of ribosomes inactivated by chloramphenicol
    init_conds('Bcm')=x(7);

    % nutrient quality s and chloramphenicol concentration h
    init_conds('s')=x(8);
    init_conds('h')=x(9); % no translation inhibition assumed by default

    % heterologous genes
    x_het=x(10:end);
    for j=1:sim.num_het
        % mRNA
        init_conds(['m_',sim.het.names{j}])=x_het(j);
        % protein
        init_conds(['p_',sim.het.names{j}])=x_het(sim.num_het+j);
    end
    % miscellaneous
    for j=1:sim.num_misc
        init_conds([sim.het.misc_names{j}])=x_het(2*sim.num_het+j);
    end
end