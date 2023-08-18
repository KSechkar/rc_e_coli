%% fig2e.m
% GROWTH PHENOMENA PREDICTION
% Figure 2: e

% Find growth rates for different fixed values of the ribosomal gene
% transcription function, compare to those produced by Flux-Parity
% Regulation. Plot ratios of ribosome mass fractions

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 4; % maximum no. iterations (checking if SS reached over first 48 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% DEFINE nutrient condiitons to explore

% vector of nurtient qualities
nutrients=flip(logspace(log10(0.05),log10(1),16));

% range of fixed ribosome transcription regulation function values
fixed_F_rs=0.005:0.005:0.99;
% corresponding reciprocals of ppGpp levels
fixed_Ts=fixed_F_rs./(sim.parameters('tau').*ones(size(fixed_F_rs))-fixed_F_rs);

%% INITIALISE arrays that will store differen growth rates
% storing growth rates - optimal and flux-parity
l_map = zeros(2,size(nutrients,2));

% storing ribosome transcription regulation function values - flux-parity
F_r_map = zeros(2,size(nutrients,2));
phi_r_map = zeros(2,size(nutrients,2));

% storing growth rates - for flux-parity and different fixed regulation function values
ls=zeros(size(nutrients,2),size(fixed_Ts,2));

% auxiliary array required to find optimal ribosomal mass fraction
phi_rs_with_fixed_Ts=zeros(size(fixed_Ts));

%% RUN simulations
for j=1:size(nutrients,2)
    % RESET
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality

    % Run with regulated ribosome transcription
    sim.parameters('is_fixed_T')=0; % F_r regulated now

    ss=get_steady(sim,Delta,Max_iter);
    [l,F_r,phi_r]=get_lFrphir(sim,ss);
    l_map(1,j)=l;
    F_r_map(1,j)=F_r;
    phi_r_map(1,j)=phi_r;

    % RESET
    sim=sim.set_default_parameters(); % reset initial consitions and parameters
    sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
    sim.init_conditions('s')=nutrients(j); % set nutrient quality 

    % Run with a range of fixed F_r values
    for i=1:size(fixed_Ts,2)
        sim=sim.set_default_parameters(); % reset initial consitions and parameters
        sim=sim.set_default_init_conditions(); % reset initial consitions and parameters
        sim.init_conditions('s')=nutrients(j); % set nutrient quality
        sim.parameters('is_fixed_T')=1; % F_r fixed now

        disp(fixed_Ts(i))
        sim.parameters('fixed_T')=fixed_Ts(i);
        ss=get_steady(sim,Delta,Max_iter);
        [l,F_r,phi_r]=get_lFrphir(sim,ss);
        ls(j,i)=l;
        phi_rs_with_fixed_Ts(j,i)=phi_r;
    end

    % finding optimal ribosome content
    [maxl,maxlpos]=max(ls(j,:));
    l_map(2,j)=maxl;
    F_r_map(2,j)=fixed_F_rs(maxlpos);
    phi_r_map(2,j)=phi_rs_with_fixed_Ts(maxlpos);
    disp(['Nutrient quality ',num2str(j),'/',num2str(size(nutrients,2)),' considered'])
end


%% FIGURE 2 E
Fe = figure('Position',[0 0 385 280]);
set(Fe, 'defaultAxesFontSize', 9)
set(Fe, 'defaultLineLineWidth', 1.5)

plot(nutrients,l_map(1,:)./l_map(2,:),'-')

xlabel('\sigma, nutrient quality','FontName','Arial');
ylabel({'\lambda:\lambda^{opt}, ratio of predicted', 'to optimal growth rate'},'FontName','Arial');

xlim([0 1])
ylim([0.85 1.15])
xticks(0:0.25:1)
yticks(0.85:0.05:1.15)

axis square
grid on
box on
hold off

%% FUNCTION for getting growth rate and F_r from the system's steady state
function [l,F_r,phi_r]=get_lFrphir(sim,ss)
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

    e=sim.form.e(par,tc); % translation elongation ratio
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    l=sim.form.l(par,e,B); % growth rate!
    
    T = tc./tu; % ratio of charged to uncharged tRNAs
    F_r = sim.form.F_r(par,T) ; % ribosome regulation
    phi_r=R.*par('n_r')./par('M');
end