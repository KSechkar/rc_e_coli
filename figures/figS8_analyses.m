% fig4cdf.m

% PREDICTING HETEROLOGOUS GENE EXPRESSION NUMERICALLY AND ANALYTICALLY
% Figure 4: c, d, f

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
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 6; % maximum no. iterations (checking if SS reached over first 2000 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed


%% INITIALISE the cell simulator, define the circuit

a_xtra=1000; % transcription rate of the heterologous gene PER PLASMID

sim=sim.load_heterologous_and_external('t7_selfact','no_ext'); % load heterologous gene expression module
sim.het.parameters('a_t7')=3000;
sim.het.parameters('K_t7')=10000;
sim.het.parameters('baseline_t7')=0.01;
sim=sim.push_het;
sim.init_conditions('s')=01;

%% SET UP the approximate estimator
approx=heterologous_approx;

%% RUN without heterologous gene expression
par=sim.parameters; % the model's parameters: saved separately for ease of notation

% values w/out het. prot.
gene_concs=containers.Map('KeyType', 'char','ValueType', 'double');
for i=1:sim.num_het
    gene_concs(['c_',sim.het.names{i}])=sim.het.parameters(['c_',sim.het.names{i}]);
    sim.het.parameters(['c_',sim.het.names{i}])=0;
end

sim = sim.push_het(); % push initial condition to main object
ss=get_steady(sim,Delta,Max_iter);

[l,e,phi_het,F_r]=get_lephihetFr(sim,ss); % get values of variables

% record
ss0 = ss;
l0=l;
e0=e;
phi_het0=phi_het;
Fr0=F_r;

% find k_x^NB (mRNA-ribosome dissoc. consts for different genes - will be needed later)
kaNB=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),0);
krNB=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),0);
kt7NB=sim.form.k(e,par('k+_t7'),par('k-_t7'),par('n_t7'),0);

% restore gene concs
for i=1:sim.num_het
    sim.het.parameters(['c_',sim.het.names{i}])=gene_concs(['c_',sim.het.names{i}]);
end
sim=sim.push_het();

%% GET approximate estimates

% define range of t7 levels
p_t7_max=(1-par('phi_q'))*par('M')/par('n_t7');
p_t7s=linspace(0,p_t7_max,100);

% estimate cell growth rates
ls_approx=l0.*(1-p_t7s*par('n_t7')./(par('M')*(1-par('phi_q'))));

% estimate gene transcription regulation functions
Fs_approx=p_t7s./(p_t7s+par('K_t7'));
Fs_approx=par('baseline_t7')+(1-par('baseline_t7'))*Fs_approx;

% estimate t7 RNAP synthesis rates
synths_approx=ls_approx.*...
    (1-par('phi_q')).*par('M')./par('n_t7').*...
    (Fs_approx.*par('c_t7').*par('a_t7')./(kt7NB.*(ls_approx+par('b_t7'))))./...
    (Fs_approx.*par('c_t7').*par('a_t7')./(kt7NB.*(ls_approx+par('b_t7')))+...
    par('c_a').*par('a_a').*ls_approx./(kaNB.*(ls_approx+par('b_a')))+Fr0.*par('c_r').*par('a_r').*ls_approx./(krNB.*(ls_approx+par('b_r'))));

%% PLOT
figure()
hold on

% plot synthesis rate
plot(p_t7s,synths_approx)
% plot degradation rate
plot(p_t7s,(ls_approx+0.18).*p_t7s)

xlabel('p_{t7} [nM]')
ylabel('dp_{t7}/dt [nM/h]')

legend('Synthesis rate','Degradation rate')

%% FUNCTION for getting growth rate, translation elongation rate and het. prot. mass fraction from the system's steady state
function [l,e,phi_het,F_r]=get_lephihetFr(sim,ss)
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

    % ribosomal gene transcription regulation function
    T=tc./tu;
    F_r=sim.form.F_r(par,T);
end
