%% fig5g.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 6: c

% Showcasing the controller's operation range in POOR culture media

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP the simulators for both cases

% 2 different setups - ideal AIF and realistic scenario
sim=cell_simulator;

sim.init_conditions('s')=0.25;

sim=sim.load_heterologous_and_external('pi_controller','constant_inducer'); % load the het. gene and ext. inp. modules

% output protein gene parameters
sim.het.parameters('c_x')=100;
sim.het.parameters('a_x')=100;

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
   
% push amended parameter values
sim=sim.push_het();

% simulation parameters
sim.tf =  72;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% DEFINE plasmid concs. to be tested

sim.het.parameters('a_dist')=500;
sim=sim.push_het();
plasmid_concs=linspace(0,1500,10);

%% RUN simulations - with controller

% initialise the array in which the results are stored
p_senss = zeros(1,size(plasmid_concs,2));

for i=1:size(plasmid_concs,2)
    disp(['Testing c_dist=',num2str(plasmid_concs(i))])

    % set plasmid concentration
    sim.het.parameters('c_dist')=plasmid_concs(i);
    sim = sim.push_het();

    % simulate!
    sim = sim.simulate_model;

    % record
    p_senss(i)=sim.x(end,9+sim.num_het+1);

    % display the growth rate of undisturbed system
    if(i==1)
        l_undist=get_l(sim,sim.x(end,:));
        disp(['Growth rate without disturbance = ',num2str(l_undist),' h^-1'])
    end
end

%% SET UP the simulations for open loop

sim_openloop=cell_simulator;

sim_openloop.init_conditions('s')=sim.init_conditions('s');

sim_openloop=sim_openloop.load_heterologous_and_external('pi_controller','constant_inducer'); % load the het. gene and ext. inp. modules

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
  
% push amended parameter values
sim_openloop=sim_openloop.push_het();

% simulation parameters
sim_openloop.tf =  sim_openloop.tf;
sim_openloop.opt = sim_openloop.opt;
sim_openloop.het.parameters('c_xtra2')=100;

%% RUN simulations - open loop
% initialise the array in which the results are stored
p_senss_openloop = zeros(1,size(plasmid_concs,2));

for i=1:size(plasmid_concs,2)
    disp(['Testing c_dist=',num2str(plasmid_concs(i))])

    % set plasmid concentration
    sim_openloop.het.parameters('c_dist')=plasmid_concs(i);
    sim_openloop = sim_openloop.push_het();

    % simulate!
    sim_openloop = sim_openloop.simulate_model;

    % record
    p_senss_openloop(i)=sim_openloop.x(end,9+sim_openloop.num_het+1);
end

%% MAKE PRE-CALCULATIONS for the plot
% on the x-axis, total transcription rate/growth
tot_trans=plasmid_concs*sim.het.parameters('a_dist');


%% ANALYTICALLY ESTIMATE the controller's dynamic range
% finding no-burden values of translation rate, dissociation constants,
% rib. gene transc. regulation function
sim_nb=cell_simulator;
sim_nb.init_conditions('s')=sim.init_conditions('s');
sim_nb=sim_nb.load_heterologous_and_external('pi_controller','step_inducer');
sim_nb.het.parameters('a_x')=0;
sim_nb.het.parameters('a_sens')=0;
sim_nb.het.parameters('a_anti')=0;
sim_nb.het.parameters('a_act')=0;
sim_nb.het.parameters('a_amp')=0;
sim_nb.het.parameters('a_dist')=0;
sim_nb=sim_nb.push_het();
sim_nb.tf =  10;
sim_nb.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
sim_nb = sim_nb.simulate_model;
[e_nb,Fr_nb,k_a_nb,k_r_nb,k_sens_nb,k_act_nb,k_amp_nb,k_dist_nb]=get_nb(sim_nb,sim_nb.x(end,:));

% calculate values of 'meaningful paramters'
par=sim.parameters; % IMPORTANT! parameters of the system where the controller's genes ARE expressed, not 'no burden' one
u=par('a_act')./par('a_anti'); % ideal value of F_anti
zeta=par('c_sens').*par('a_sens');

%% 
% calculating the value of lambda
lambda_est = ...
    e_nb./par('M') .* ...
    Fr_nb.*par('c_r').*par('a_r') .* ...
    par('n_sens').*k_sens_nb./(par('n_r').*k_r_nb) .* ...
    par('K_dna(anti)-sens')./(zeta) .* ...
    (1-u)./u;

% calculating the value of D
D_est = ...
    1 + ...
    (lambda_est.*zeta) ./ ...
    ((lambda_est+par('b_sens')).*k_sens_nb) ./ ...
    (par('n_sens')./par('M') .* par('K_dna(anti)-sens') .* ((1-u)./u) );

% find sum of steady-state mj/kj_nb values
maka_nb = ...
    lambda_est.*par('c_a').*par('a_a') ./ ...
    ((lambda_est+par('b_a')).*k_a_nb);
mrkr_nb = ...
    lambda_est.*Fr_nb.*par('c_r').*par('a_r') ./ ...
    ((lambda_est+par('b_r')).*k_r_nb);
msensksens_nb = ...
    lambda_est.*Fr_nb.*zeta ./ ...
    ((lambda_est+par('b_sens')).*k_sens_nb);
sum_mjkj_nb=maka_nb+mrkr_nb+msensksens_nb;

% find Omega
% Omega = ...
%     (1-par('phi_q')).*par('M').*zeta.*lambda_est ./ ...
%     (par('K_dna(anti)-sens').*par('n_sens').*k_sens_nb.*(lambda_est+par('b_sens'))) .* ...
%     (1-u)./u - ...
%     sum_mjkj_nb;
Omega=(1-par('phi_q')).*(D_est-1)-sum_mjkj_nb;

%% FIND when the controller is expected to be out of its dynamic range

% find steady-state m_dist/k_dist_nb values, assuming the same
% growth rate is maintained by the controller
mdistkdist_nbs = ...
    lambda_est.*plasmid_concs.*par('a_dist') ./ ...
    ((lambda_est+par('b_dist')).*k_dist_nb); % an array!!!
sums_x_and_dist=mdistkdist_nbs;

% find when it exceeds Omega
for i=1:size(plasmid_concs,2)
    i_exceed=i;
    if(sums_x_and_dist(i)>Omega)
        break
    end
end

%% Main Figure h

Fh = figure('Position',[0 0 250 186]);
set(Fh, 'defaultAxesFontSize', 9)
set(Fh, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(plasmid_concs,p_senss_openloop./p_senss_openloop(1),'Color',[0.6350 0.0780 0.1840])
plot(plasmid_concs,p_senss./p_senss(1),'Color',[0 0.4470 0.7410])


% mark the dynamic range of the controller
plot([plasmid_concs(i_exceed) plasmid_concs(i_exceed)],[0 2],'Color','k','LineStyle',':')

xlabel('c_{dist}, disturbing gene conc. [nM]','FontName','Arial');
ylabel({'p_{sens}:p_{sens}^0, relative', 'sensor prot. conc.'},'FontName','Arial')

ylim([0.3 1.1])
xlim([0 1500])
%xticks(0:1e5:3e5)
% yticks(0.4:0.2:1.2)

legend({'open loop','closed loop'},'FontName','Arial','FontSize',8,'Location','southwest')

grid 
box on
axis square
hold off

%% NO BURDEN VALUES FUNCTION
%

function [e,Fr,k_a,k_r,k_sens,k_act,k_amp,k_dist]=get_nb(sim,ss)   
    par=sim.parameters;
    % STATE VECTOR TO SINGLE VARIABLES
    m_a = ss(1);
    m_r = ss(2);
    p_a = ss(3);
    R = ss(4);
    tc = ss(5);
    tu = ss(6);
    Bcm = ss(7);
    s = ss(8);
    h = ss(9);
    x_het=ss(10 : (9+2*sim.num_het) );

    % USEFUL PRE-CALCULATIONS
    % translation elongation rate
    e=sim.form.e(par,tc);

    % ribosome inactivation rate due to chloramphenicol
    kcmh=par('kcm').*h;

    % ribosome dissociation constants
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
    % heterologous genes
    k_het=ones(1,sim.num_het);
    if(sim.num_het>0)
        for j=1:sim.num_het
            k_het(j)=sim.form.k(e,...
            sim.parameters(['k+_',sim.het.names{j}]),...
            sim.parameters(['k-_',sim.het.names{j}]),...
            sim.parameters(['n_',sim.het.names{j}]),...
            kcmh);
        end
    end

    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:sim.num_het)./k_het))./...
        (1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % growth rate
    l=sim.form.l(par,e,B);

    p_sens=x_het(6); % get sensor conc.
    
    % RECORD RIBOSOMAL GENE TRANSC. REGULATION FUNCTION
    Fr=sim.form.F_r(par,T);

    % RECORD DISSOCIATION CONSTANTS
    k_sens=k_het(1);
    k_act=k_het(3);
    k_amp=k_het(4);
    k_dist=k_het(5); 
end

%% GROWTH RATE CALCULATION FUNCTION
function l=get_l(sim,ss)   
    par=sim.parameters;
    % STATE VECTOR TO SINGLE VARIABLES
    m_a = ss(1);
    m_r = ss(2);
    p_a = ss(3);
    R = ss(4);
    tc = ss(5);
    tu = ss(6);
    Bcm = ss(7);
    s = ss(8);
    h = ss(9);
    x_het=ss(10 : (9+2*sim.num_het) );

    % USEFUL PRE-CALCULATIONS
    % translation elongation rate
    e=sim.form.e(par,tc);

    % ribosome inactivation rate due to chloramphenicol
    kcmh=par('kcm').*h;

    % ribosome dissociation constants
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
    % heterologous genes
    k_het=ones(1,sim.num_het);
    if(sim.num_het>0)
        for j=1:sim.num_het
            k_het(j)=sim.form.k(e,...
            sim.parameters(['k+_',sim.het.names{j}]),...
            sim.parameters(['k-_',sim.het.names{j}]),...
            sim.parameters(['n_',sim.het.names{j}]),...
            kcmh);
        end
    end

    T=tc./tu; % ratio of charged to uncharged tRNAs
    D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:sim.num_het)./k_het))./...
        (1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q

    % growth rate
    l=sim.form.l(par,e,B);
end