% fig3def_figS2.m

% PREDICTING HETEROLOGOUS GENE EXPRESSION NUMERICALLY AND ANALYTICALLY
% Figure 3: c, d, f

% Analyse the effects of heterologous gene expression: investigate growth
% rate and het. prot. mas fraction (Main Figure b,c), as well as charged/uncharged 
% tRNA conc. ratio, traslation el. rate, tRNA aminoacylation rate and  
% ribosomal genetranscription regulation as a function of 
% het. gene transcription. Determine which heterologous protein mass fraction 
% leads to the highest protein production rates in a population of bacterai
% when the culturing starts (Main Figure d). 

%% CLEAR all variables

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 20; % single integraton step timeframe
Delta = 0.1; % threshold that determines if we're in steady state
Max_iter = 100; % maximum no. iterations (checking if SS reached over first 2000 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed


%% DEFINE gene expression paramaters

a_xtra=1000; % transcription rate of the heterologous gene PER PLASMID

sim=sim.load_heterologous_and_external('one_constit','no_ext'); % load heterologous gene expression module
sim.het.parameters('a_xtra')=a_xtra;
sim=sim.push_het;

plasmid_concs = logspace(0,3.04,100); % range of plasmid concentrations

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

%% RUN without heterologous gene expression
par=sim.parameters; % the model's parameters: saved separately for ease of notation

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

%% MAIN FIGURE C - het. prot. mass fraction as a function of burden
approx_colour=[0.6350 0.0780 0.1840];

Fc = figure('Position',[0 0 385 290]);
set(Fc, 'defaultAxesFontSize', 9)
set(Fc, 'defaultLineLineWidth', 1)

hold on

plot([0,plasmid_concs*a_xtra]./kxNB,[phi_het0,phi_hets],['r','-']) % plot model predictions
plot([0,plasmid_concs*a_xtra]./kxNB,[phi_het0_approx,phi_hets_approx], ...
    'Color',approx_colour,'LineStyle','--') % plot approximation

xlabel('\xi, translational burden','FontName','Arial');
ylabel('\phi_x, het. prot. mass fraction','FontName','Arial');

legend('Simulation results','Analytical estimates','Location','northeast','FontName','Arial');

ylim([0 0.35])
yticks([0:0.1:0.3,0.35])

grid on
box on
axis square
hold off 

%% MAIN FIGURE D cell growth as a function of burden

Fb = figure('Position',[0 0 385 287]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1)
hold on

plot([0,plasmid_concs*a_xtra]./kxNB,[l0,ls],['r','-']) % plot model predictions
plot([0,plasmid_concs*a_xtra]./kxNB,[l0_approx,ls_approx],...
    'Color',approx_colour,'LineStyle','--') % plot approximation

xlabel({'\xi, translational burden'},'FontName','Arial');
ylabel('\lambda, growth rate [1/h]','FontName','Arial');

legend('Simulation results','Analytical estimates','Location','northeast','FontName','Arial')
grid on
box on
axis square
hold off

%% MAIN FIGURE F - total heterologous protein production rate at t=0 in a population of cells (mu_het)

delta=0.25;
Fd = figure('Position',[0 0 385 290]);
hold on

set(Fd, 'defaultAxesFontSize', 9)
set(Fd, 'defaultLineLineWidth', 1)

mu_het0=0; % initialise
mu_het0_approx=0; % initialise


% make delta the same dimension as growth rate
deltas=delta*ones(size(ls));

% actual value of mu_het
mu_hets=(ls-deltas).*phi_hets.*par('M');

% approximation of mu_het
mu_hets_approx=approx.ss_mu_het(plasmid_concs,sss,ss0,e0,l0,sim,delta);

plot([0,plasmid_concs*a_xtra]./kxNB,[mu_het0,mu_hets],['r','-']) % plot model predictions
plot([0,plasmid_concs*a_xtra]./kxNB,[mu_het0_approx,mu_hets_approx], ...
    'Color',approx_colour,'LineStyle','--') % plot approximation

% draw a line to mark phi_het that maximises production (NUMERICAL RESULT)
[~,index_of_max_mu]=max(mu_hets); % get index
mu_max_approx=mu_hets(index_of_max_mu); % find the corresponding mu
% get the value
plot([plasmid_concs(index_of_max_mu)*a_xtra./kxNB plasmid_concs(index_of_max_mu)*a_xtra./kxNB],[0 mu_max_approx], ...
    'Color',[0,0.8,0.8],'LineStyle','-','LineWidth',0.75) % plot

% draw a line to mark xi that maximises production  (ANALYTICAL PREDICTION)
xi_max_approx=approx.xi_max(ss,ss0,e0,l0,sim,delta); % estimate optimal burden
phi_het_max_approx=approx.phi_het_max(l0,sim,delta); % estimate optimal mass fraction of synthetic protein
mu_max_approx=par('M').*phi_het_max_approx.*...
    (l0.*(1-phi_het_max_approx./(1-par('phi_q')))-delta); % find corresponding mu
plot([xi_max_approx xi_max_approx],[0 mu_max_approx], ...
    'Color',[0,0.8,0.8],'LineStyle','--','LineWidth',0.75) % plot

xlabel('\xi, translational burden','FontName','Arial');
ylabel({'rate of total protein', 'production [aa/h/cell]'},'FontName','Arial');
legend('Simulation results','Analytical estimates','Location','northeast','FontName','Arial');

ylim([0 17.5*10^7])

yticks((0:2.5:17.5)*10^7)

grid on
box on
axis square
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
        phi_het=phi_het+ss_het(sim.num_het+i).*par(['n_',sim.het.names{i}])./par('M');
    end
end
