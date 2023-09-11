%% fig3bcd.m
% GROWTH PHENOMENA PREDICTION BY THE MODEL
% Figure 3: b,c,d

% Compare model predictions of the cell's growth rate, ribosomal mass
% fraction, translation elongation rate and ppGpp levels with experimental data from the 
% last 55 years of measurements, compiled by Chure et al. 2022

%% CLEAR

addpath(genpath('..'))

clear
close all

%% LOAD experimental data
% ribosomal mass fractions vs growth rates
[dataset,captions,~] = xlsread('data/exp_meas_ribosomes.csv');
data_rib(:,1) = dataset(:,1); % het. prot. mass fraction
data_rib(:,2) = dataset(:,2); % (growth rate w/ het. prot):(growth rate w/out)

% translation elongation rates vs growth rates
[dataset,captions,~] = xlsread('data/exp_meas_elongation.csv');
data_el(:,1) = dataset(:,1); % growth rate (1/h)
data_el(:,2) = dataset(:,2); % translation elongation rate

% ppGpp levels vs growth rates
[dataset,captions,~] = xlsread('data/exp_meas_ppGpp.csv');
data_ppGpp(:,1) = dataset(:,1); % growth rate (1/h)
data_ppGpp(:,2) = dataset(:,2); % ppGpp conc. relative to reference (l=0.9/h)

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 6; % maximum no. iterations (checking if SS reached over first 60 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% SPECIFY range of nutrient qualit for which we run the simulation

nutr_quals=logspace(log10(0.01),log10(1),32);

%% INITIALISE the arrays where model predictions will be stored

l_map=zeros(size(nutr_quals));
phir_map=zeros(size(nutr_quals));
el_map=zeros(size(nutr_quals));
ppGpp_map=zeros(size(nutr_quals));

%% GET model predictions

for j=1:size(nutr_quals,2)
    disp(nutr_quals(j))

    % reset simulator
    sim.parameters=cell_params();
    theta=[1.03184, 4046.89, 1239.69, 0.000356139];
    sim.parameters('a_a') = 3.89e5; % metabolic prot. transcription rate (/h)
    sim.parameters('a_r') = sim.parameters('a_a').*theta(1); % ribosome transcription rate (/h) - rescaled!
    sim.parameters('nu_max') = theta(2); % max metabolic rate (/h)
    sim.parameters('K_e') = theta(3); % elongation rate Hill constant (nM)
    sim.parameters('K_nut') = theta(3); % tRNA charging rate Hill constant (nM)
    sim.parameters('kcm') = theta(4); % tRNA charging rate Hill constant (nM)
    
    sim.init_conditions=cell_init_conds(sim.parameters); % reset intial conditions
    
    % set nutrient quality
    sim.init_conditions('s')=nutr_quals(j);
    sim=sim.push_het();

    % Run
    sim.parameters('is_fixed_F_r')=0; % F_r regulated now
    ss = get_steady(sim,Delta,Max_iter);
    [l, e, phi_r,T] = get_lephirT(sim,ss);
    l_map(j)=l;
    phir_map(j)=phi_r;
    el_map(j)=e./3600; % record in aa/s, NOT aa/h
    ppGpp_map(j)=1./T;
end

%% FIND the reference ppGpp concentration (for l=1/h)
distance_to_ref=100; % initialise the distance to reference growth rate with unreasonably high number
for i=1:size(l_map,2)
    if(abs(l_map(i)-1)<abs(distance_to_ref)) % if the current distance is the smallest so far, this is the neww reference
        distance_to_ref=l_map(i)-1;
        closest_i=i;
    end
end

ppGpp_ref=ppGpp_map(closest_i); % reference ppGpp level

%% FIGURE 2 B

Fb = figure('Position',[0 0 385 280]);
set(Fb, 'defaultAxesFontSize', 9)

hold on

plot(data_rib(:,1),data_rib(:,2),'.','Color','b')
plot(l_map,phir_map,'-','Color','r','LineWidth',1.5)
plot([0.3 0.3],[0 0.3],':','Color','k')

ylabel('\phi_r, ribosomal mass fraction','FontName','Arial');
xlabel('\lambda, growth rate [1/h]','FontName','Arial');
legend('Experimental data','Model predictions',...
    'Location','southeast','FontName','Arial')
grid on
axis square
box on
hold off

%% FIGURE 2 C

Fc = figure('Position',[0 0 385 280]);
set(Fc, 'defaultAxesFontSize', 9)
hold on

plot(data_el(:,1),data_el(:,2),'.','Color','b')
plot(l_map,el_map,'-','Color','r','LineWidth',1.5)
%plot([0.3 0.3],[6 20],':','Color','k')

ylabel('\epsilon, translation elong. rate [aa/s]','FontName','Arial');
xlabel('\lambda, growth rate [1/h]','FontName','Arial');
legend('Experimental data','Model predictions',...
    'Location','southeast','FontName','Arial')
grid on
axis square
box on
hold off

%% FIGURE 2 D

Fd = figure('Position',[0 0 385 281]);
set(Fd, 'defaultAxesFontSize', 9)
hold on

plot(data_ppGpp(:,1),data_ppGpp(:,2),'.','Color','b')
plot(l_map,ppGpp_map./ppGpp_ref,'-','Color','r','LineWidth',1.5)
%plot([0.9 0.9],[10^(-1) 10^(1.5)],'--','Color','m','LineWidth',1)

xlabel('\lambda, growth rate [1/h]');
ylabel('Normalized ppGpp conc.');
legend('Experimental data','Model predictions',... 'Rel. rsd',
    'Location','northeast');

xlim([0 2.5])
ylim([10^(-1) 10^(1.5)])
set(gca, 'YScale', 'log')
grid on

axis square
box on
hold off

%% FUNCTION for getting growth rate, translation elongation rate and rib. mass fraction from the system's steady state

function [l, e, phi_r,T]=get_lephirT(sim,ss)
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

    % rib. mass fraction
    phi_r=R.*par('n_r')./par('M');

    % inverse of ppGpp level
    T=tc./tu;
end
