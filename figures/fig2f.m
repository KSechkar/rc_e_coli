%% fig2bcd.m
% GROWTH PHENOMENA PREDICTION BY THE MODEL
% Figure 2: f

% Compare model predictions for the dynamics of nutrient upshifts with
% experimental data collected by Erickson et al. and processed by Chure et
% al.
%% CLEAR

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% getting steady state ( to find sigma values for different media)
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 4; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed

%% SPECIFY which upshift experiments are considered
% ribosomal mass fractions vs growth rates
[dataset,captions,~] = xlsread('data/exp_meas_upshifts.csv');
% upshift_media={'arabinose','xylose','glycerol','gluconate'};
% upshift_ls={0.7,0.76, 0.85, 0.94};
% base_medium='pyruvate';
% base_l=0.6;

upshift_media={'arabinose','xylose','glycerol','glucose'};
upshift_ls={0.55,0.73, 0.75, 0.86};
base_medium='pyruvate';
base_l=0.45;

% upshift_media={'arabinose','gluconate'};
% upshift_ls={0.55,0.85};
% base_l=0.45;

% colours for plots
upshift_colours={[0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880],...
    [0.4940 0.1840 0.5560],...
    [0.9290 0.6940 0.1250],...
    [0.8500 0.3250 0.0980],...
    [0 0.4470 0.7410]};

%% LOAD experimental data
% initialise the data storage structure
empty_data_for_one_medium.ls=[];
empty_data_for_one_medium.ts=[];
data={empty_data_for_one_medium, empty_data_for_one_medium, empty_data_for_one_medium, empty_data_for_one_medium};

% read the experimental data, splitting them by the upshift media
for i = 1:size(dataset,1)
    % ! only considering upshifts from SUCCINATE (original growth rate 0.6)
    if(dataset(i,3)==base_l)

        % find which upshift medium the current measuremnts relate to
        for j = 1:size(upshift_ls,2)
            if(dataset(i,4)==upshift_ls{j})
                % record the instantaneous growth rate and time of measurement
                data{j}.ls=[data{j}.ls; dataset(i,2)];
                data{j}.ts=[data{j}.ts; dataset(i,1)];
                break
            end
        end
    end
end

%% FIND steady-state growth rates for a wide range of nutrient qualities

% define nutr. qual. range
nutr_quals=linspace(0.01,1,1000);

% initialise growth rates array
l_map=zeros(size(nutr_quals));
for j=1:size(nutr_quals,2)
    disp(nutr_quals(j))
    % reset simulator
    sim.parameters=cell_params();
    sim.init_conditions=cell_init_conds(sim.parameters); % reset intial conditions
    % set nutrient quality
    sim.init_conditions('s')=nutr_quals(j);
    sim=sim.push_het();
    % Run
    ss = get_steady(sim,Delta,Max_iter);
    % get growth rate
    l = get_l(sim,ss);
    l_map(j)=l;
end

%% HENCE, RETRIEVE sigma values for different nutrients in the experiment

% original medium
base_sigma=nutr_quals(closest_pos(l_map,base_l));

% upshift media
upshift_sigmas={0,0,0,0}; % initialise the cell array for storing sigma values
for upshift_medium = 1:size(upshift_media,2)
    upshift_sigmas{upshift_medium}=nutr_quals(closest_pos(l_map,upshift_ls{upshift_medium}));
end

disp('Sigma values for different media retrieved')


%% SIMULATE nutrient upshifts

% specify time before and after upshift for which we run the simulation
t_before=1.5;
t_after=3;

% initialise cell array of model predictions
model={empty_data_for_one_medium, empty_data_for_one_medium, empty_data_for_one_medium, empty_data_for_one_medium};

for upshift_medium = 1:size(upshift_media,2)
    % reset simulator
    sim.parameters=cell_params();
    sim.init_conditions=cell_init_conds(sim.parameters);

    % set nutrient quality
    sim.init_conditions('s')=base_sigma; % base medium before upshift
    sim=sim.push_het();

    % get steady state as the starting point
    x_steady_base=get_steady(sim,Delta,Max_iter);

    % simulate the model before nutrient upshift
    sim.init_conditions=x2init_conds(x_steady_base);
    sim=sim.push_het();
    sim.tf=t_before;
    sim = sim.simulate_model;

    % find and record growth rates and time axis
    ts=sim.t-t_before;
    ls=zeros(size(sim.t)); % initialise
    for i = 1:size(sim.t)
        ls(i)=get_l(sim,sim.x(i,:));
    end

    % simulate the model after nutrient upshift
    sim.init_conditions=x2init_conds(sim.x(end,:));
    sim.parameters('is_upshift')=1;
    sim.parameters('s_postshift')=upshift_sigmas{upshift_medium}; % set the postshift nutr. qual. !
    sim.init_conditions('s')=sim.init_conditions('s')+...
        (sim.parameters('s_postshift')-sim.init_conditions('s' )).*...
        base_l./upshift_ls{upshift_medium};
    sim.init_conditions('s')=sim.parameters('s_postshift').*...
        base_l./upshift_ls{upshift_medium};
    sim=sim.push_het;
    sim.tf=t_after;
    sim=sim.simulate_model;

    % find and record growth rates and time axis
    ts=[ts; sim.t];
    ls_after=zeros(size(sim.t)); % initialise
    for i = 1:size(sim.t)
        ls_after(i)=get_l(sim,sim.x(i,:));
    end
    ls=[ls; ls_after];

    % record times and growth rates in permanent memory
    model{upshift_medium}.ts=ts;
    model{upshift_medium}.ls=ls;
end

disp('Nutrient upshift simulation successful')

%% FIGURE 2 F

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

for upshift_medium = 1:size(upshift_media,2)
    % experimental data
    plot(data{upshift_medium}.ts,data{upshift_medium}.ls, ...
        '.','Color',upshift_colours{upshift_medium})

    % model predictions
    plot(model{upshift_medium}.ts,model{upshift_medium}.ls, ...
        '-','Color',upshift_colours{upshift_medium},'LineWidth',1)
end

xlabel('Time since upshift [h]')
ylabel('\lambda, growth rate [1/h]');
legend('Experimental data','Model predictions',... 'Rel. rsd',
    'Location','southeast');

xlim([-t_before t_after])
ylim([0.2 1])
grid on

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
function init_conds=x2init_conds(x)
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

    % no heterologous genes in this scenario
end