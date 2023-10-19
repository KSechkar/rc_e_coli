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

%% LOAD experimental data

t_for_estimate=1;

% OD - cel conc. conversion
OD_to_ccells=(8e8*1000)/(6.02e23)*1e9; % OD of bacteria to concnentation of cells in nM

[dataset,captions,~] = xlsread('data/exp_meas_depletion.csv');
data.ODs=[];
data.ts=[];
c_nutr=0.22e6;
K_nutr=5e3;


% read the experimental data, splitting them by the upshift media
for i = 2:size(dataset,1)
    caption_cell=captions(i,3);
    % ! only considering 0.22 mM glucose
    if(strcmp(caption_cell{1},'0.22 mM Glucose'))
        % record the instantaneous growth rate and time of measurement
        data.ODs=[data.ODs; dataset(i,2)];
        data.ts=[data.ts; dataset(i,1)];
    end
end

% recalibrate the origin so that time axis starts at 0
data.ts=data.ts-min(data.ts);

%% ESTIMATE the original growth rate based on first hour of measurements
% estimate the original growth rate and cell conc. to find medium's nutrient quality
ts_for_init=data.ts(data.ts<t_for_estimate);
ODs_for_init=data.ODs(data.ts<t_for_estimate);
fitted=fit(ts_for_init,ODs_for_init,fittype('exp1'));
l_init=fitted.b;
ccells_init=fitted.a*OD_to_ccells;

% find the final cell conc. to find the medium's carrying capacity
ts_for_last=data.ts(data.ts>max(data.ts)-t_for_estimate);
ODs_for_last=data.ODs(data.ts>max(data.ts)-t_for_estimate);
ccells_max=mean(ODs_for_last)*OD_to_ccells; % final OD estimated as mean of values over last t_for_estimate h

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% getting steady state ( to find sigma values for different media)
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 6; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed

sim=sim.load_heterologous_and_external('one_switch','no_ext');
sim.het.parameters('c_switch1')=0; % no heterologous gene expression
sim.parameters('f1')=1e5; % excess of inducer
sim=sim.push_het;

%% FIND steady-state growth rates for a wide range of nutrient qualities

% define nutr. qual. range
nutr_quals=linspace(0.01,1,100);

% initialise growth rates array
l_map=zeros(size(nutr_quals));
for j=1:length(nutr_quals)
    disp(nutr_quals(j))
    % reset simulator
    sim.parameters=cell_params();
    sim.init_conditions=cell_init_conds(sim.parameters); % reset intial conditions
    % set nutrient quality
    sim.init_conditions('s')=nutr_quals(j)*(c_nutr/(c_nutr+K_nutr)); % assuming nutrient concentration almost unchanging within the first hour
    sim=sim.push_het();
    % Run
    ss = get_steady(sim,Delta,Max_iter);
    % get growth rate
    l = get_l(sim,ss);
    l_map(j)=l;
end

% HENCE, RETRIEVE sigma values for different nutrients in the experiment

sigma=nutr_quals(closest_pos(l_map,l_init));
disp('Sigma values for different media retrieved')

%% GET the initial condition (steady state)
% before dilution, cells growing in a largeb excess (110000 nM) of glucose
sim.init_conditions('s')=sigma;

% get steady state as the starting point
x_steady_base=get_steady(sim,Delta,Max_iter);

% simulate the model before nutrient upshift
sim.init_conditions=x2init_conds(sim,x_steady_base);
sim=sim.push_het();

%% SIMULATE nutrient depletion

% simulating nutirent depletionnow
sim.parameters('is_depletion')=1;
sim.parameters('ccells_max')=ccells_max;
sim.parameters('K_nutr')=K_nutr;
sim.het.init_conditions('c_cells')=fitted.a*OD_to_ccells;
sim.het.init_conditions('c_nutr')=c_nutr; % 0.22 mM, 1 mM=10^6 nM

sim=sim.push_het();

% simulate
sim.tf=10;
sim=sim.simulate_model();


%% FIGURE S8a

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

% experimental data
plot(data.ts,data.ODs,'.')

% model predictions
plot(sim.t,sim.x(:,end-1)./OD_to_ccells, ...
    '-','LineWidth',1.5)

xlabel('Time since upshift [h]')
ylabel('OD of culture');

xlim([0 7])
grid on
set(gca, 'YScale', 'log')
axis square
box on
hold off

%% FIGURE S8b

% Fd = figure('Position',[0 0 385 281]);
% set(Fd, 'defaultAxesFontSize', 9)
% hold on
% 
% 
% % model predictions
% plot(sim.t,sim.x(:,end), ...
%     '-','LineWidth',1.5)
% 
% xlabel('Time since upshift [h]')
% ylabel('Nutrient conc., nM');
% 
% xlim([0 sim.tf])
% grid on
% % set(gca, 'YScale', 'log')
% axis square
% box on
% hold off

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