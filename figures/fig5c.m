%% fig5c.m
% PREDICTING HETEROLOGOUS GENE EXPRESSION NUMERICALLY AND ANALYTICALLY
% Figure 5: c

% Change in cell growth rate vs heterologous protein mass fraction

%% CLEAR all variables

addpath(genpath('..'))

clear
close all

%% LOAD experimental data
% het. prot. mass fractions vs growth rates
[dataset,captions,~] = xlsread('data/exp_meas_hetexp.csv');
data_het(:,1) = dataset(:,1); % het. prot. mass fraction
data_het(:,2) = dataset(:,3); % (growth rate w/ het. prot):(growth rate w/out)

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 40; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed


%% DEFINE gene expression paramaters

a_xtra=1000; % transcription rate of the heterologous gene PER PLASMID

sim=sim.load_heterologous_and_external('one_constit','no_ext'); % load heterologous gene expression module
sim.het.parameters('a_xtra')=a_xtra;
sim=sim.push_het;

plasmid_concs = logspace(0,5,100); % range of plasmid concentrations


%% SET UP the approximate estimator
approx=heterologous_approx;

%% INITIALISE arrays where obtained values of phys. variables will be stored

% values without het. prot. exp.
ss0=zeros((size(sim.init_conditions,1)+size(sim.het.init_conditions,1)),1); % steady states of the system
l0=0; % growth rates
e0=0;  % elongation rates
phi_het0=0; % heterologous protein fractions

% numerically obtained values
sss=zeros(size(sim.init_conditions,1),size(plasmid_concs,2)); % steady states of the system
ls = zeros(1,size(plasmid_concs,2)); % growth rates
es = zeros(1,size(plasmid_concs,2)); % elongation rates
phi_hets = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions

% approximate estimates
l0_approx=0;
phi_het0_approx=0;
ls_approx = zeros(1,size(plasmid_concs,2)); % growth rates
phi_hets_approx = zeros(1,size(plasmid_concs,2)); % heterologous protein fractions

%% RUN withput heterologous gene expression
% values w/out het. prot.
for i=1:sim.num_het
    sim.het.parameters(['c_',sim.het.names{i}])=0;
end

sim = sim.push_het(); % push initial condition to main object
ss=get_steady(sim,Delta,Max_iter);

[l,e,phi_het]=get_lephihet(sim,ss); % get values of variables

% record
ss0 = ss;
l0=l;
e0=e;
phi_het0=phi_het;

% find k_xtra^NB (mRNA-ribosome dissoc. const for het. gene - will be needed later)
par=sim.het.parameters;
kxNB=sim.form.k(e,par('k+_xtra'),par('k-_xtra'),par('n_xtra'),0);

%% GET approximate estimates
l0_approx=approx.ss_l(0,ss0,ss0,e0,sim);
ls_approx=approx.ss_l(plasmid_concs,sss,ss0,e0,sim);
phi_hets_approx=approx.ss_phi_het(plasmid_concs,sss,ss0,e0,sim);

%% OBTAIN numerical predictions

% get actual model predictions
for i=1:size(plasmid_concs,2)
    % disp(plasmid_concs(i))
    for j=1:sim.num_het
        sim.het.parameters(['c_',sim.het.names{j}])=plasmid_concs(i);
    end
    sim = sim.push_het(); % reset initial condition
    ss=get_steady(sim,Delta,Max_iter);
    [l,e,phi_het]=get_lephihet(sim,ss); % get desired values

    % record
    sss(:,i) = ss;
    ls(i)=l;
    es(i)=e;
    phi_hets(i)=phi_het;
end

%% FIGURE 3 e 
approx_colour=[0 0.8 0.8];
Fig = figure('Position',[0 0 385 290]);
set(Fig, 'defaultAxesFontSize', 9)
hold on

plot([0,phi_hets],...
    [l0,ls]/l0,['r','-'],'LineWidth',1.5) % plot model predictions
plot([0,phi_hets_approx],[l0_approx,ls_approx]/l0_approx,...
    'Color',approx_colour,'LineStyle','--','LineWidth',1.5) % plot approximate estimates

plot(data_het(:,1),data_het(:,2),'b.') % plot experimental data

plot([0.25 0.25],[0 1.05],':','Color','k','LineWidth',1)
xlabel('\phi_x, het. prot. mass fraction','FontName','Arial');
ylabel('\lambda:\lambda^{NB}, relative growth rate','FontName','Arial');
xlim([0 0.41])
ylim([0 1.05])

legend('Simulation results','Analytical estimates','Experimental data',...
    'Location','northeast','FontName','Arial')

axis square
grid on
box on
hold off

%% FUNCTION for getting growth rate, translation elongation rate and het. prot. mass fraction from the system's steady state
function [l,e,phi_het]=get_lephihet(sim,ss)
    % evaluate growth
    par=sim.parameters;
    m_a = ss(1); % metabolic gene mRNA
    m_r = ss(2); % ribosomal mRNA
    p_a = ss(3); % metabolic proteins
    R = ss(4); % operational ribosomes
    tc = ss(5); % charged tRNAs
    tu = ss(6); % uncharged tRNAs
    Bm = ss(7); % inactivated ribosomes
    s=ss(8); % nutrient quality
    h=ss(9); % chloramphenicol conc.
    ss_het=ss(10:9+2*sim.num_het); % heterologous genes

    e=sim.form.e(par,tc); % translation elongation ratio
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    l=sim.form.l(par,e,B); % growth rate!
    e=sim.form.e(par,tc); % translation elongation rate!

    % heterologous protein mass fraction
    phi_het=0;
    for i=1:sim.num_het
        phi_het=phi_het+ss_het(sim.num_het+i).*sim.parameters(['n_',sim.het.names{i}])./sim.parameters('M');
    end
end
