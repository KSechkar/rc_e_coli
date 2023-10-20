% fig4de.m

% EMERGENT BISTABILITY OF A NON-COOPERATIVE SELF-ACTIVATOR
% Figure 4: d, e

% Despite rpomoting its own expression non-cooperatively, a T7 RNAP gene
% controller by a constitutive T7 RNAP promoter may exhibit bistability due
% to the effect of synthetic gene expression on the host cell's growth
% rate

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

%% SIMULATE the system startung from a low-t7 initial condition
sim.het.init_conditions('m_t7')=5;
sim.het.init_conditions('p_t7')=0;
sim.init_conditions('p_a')=(sim.init_conditions('p_a')*sim.parameters('n_a')-...
    sim.het.init_conditions('p_t7').*sim.het.parameters('n_t7'))./sim.parameters('n_a'); % adjust initial condition to have constant cell mass
sim.het.parameters('c_t7')=c_t7;
sim=sim.push_het;
sim.init_conditions('s')=s;
sim=sim.simulate_model();
% record simulation outcomes
t_fromlow=sim.t;
x_fromlow=sim.x;
% calculate growth rates
l_fromlow=get_ls(sim);

%% SIMULATE the system startung from a HIGH-t7 initial condition
sim.het.init_conditions('m_t7')=50;
sim.het.init_conditions('p_t7')=0;
sim=sim.push_het();
sim.init_conditions('p_a')=(sim.init_conditions('p_a')*sim.parameters('n_a')-...
    sim.het.init_conditions('p_t7').*sim.het.parameters('n_t7'))./sim.parameters('n_a'); % adjust initial condition to have constant cell mass
sim.init_conditions('s')=s;

% simulate
sim=sim.simulate_model();
% record simulation outcomes
t_fromhigh=sim.t;
x_fromhigh=sim.x;
% calculate growth rates
l_fromhigh=get_ls(sim);

%% FIGURE a - protein levels

Fa = figure('Position',[0 0 215 193]);
set(Fa, 'defaultAxesFontSize', 9)
set(Fa, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(t_fromlow,x_fromlow(:,11))
plot(t_fromhigh,x_fromhigh(:,11))

xlabel('time [h]','FontName','Arial');
ylabel('m_{t7}, T7-RNAP conc. [nM]','FontName','Arial')

legend({'m_{t7}(0 h)=25 nM','p_{t7}(0 h)=50 nM'})

set(gca, 'YScale', 'log')
ylim([1 1e4])
xlim([0 sim.tf])

grid 
box on
axis square
hold off

%% FIGURE b - growth rates

Fb = figure('Position',[0 0 215 193]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

hold on

% plot p_x values with and without the controller
plot(t_fromlow,l_fromlow)
plot(t_fromhigh,l_fromhigh)

xlabel('time [h]','FontName','Arial');
ylabel('\lambda, growth rate [1/h]','FontName','Arial')

legend({'p_{t7}(0 h)=25 nM','p_{t7}(0 h)=50 nM'})

xlim([0 sim.tf])

grid 
box on
axis square
hold off

%% FUNCTION for getting growth rates over time

function ls=get_ls(sim)
    ls=zeros(size(sim.t));
    par=sim.parameters;
    for i=1:length(sim.t)
        % STATE VECTOR TO SINGLE VARIABLES
        m_a = sim.x(i,1);
        m_r = sim.x(i,2);
        p_a = sim.x(i,3);
        R = sim.x(i,4);
        tc = sim.x(i,5);
        tu = sim.x(i,6);
        Bcm = sim.x(i,7);
        s = sim.x(i,8);
        h = sim.x(i,9);
        x_het=sim.x(i,10 : (9+2*sim.num_het) );
    
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

        % calculate RC denbominator with protein degradation
        prot_deg_rate_div_eR=par('d_t7').*x_het(2).*par('n_t7')./(e.*R);
        sum_mk_not_q=m_a./k_a+m_r./k_r+sum(x_het(1:sim.num_het)./k_het);
        mq_div_kq=(par('phi_q').*(1-prot_deg_rate_div_eR).*sum_mk_not_q-par('phi_q').*prot_deg_rate_div_eR)./...
            (1-par('phi_q').*(1-prot_deg_rate_div_eR));
        D=1+sum_mk_not_q+mq_div_kq;
    
        T=tc./tu; % ratio of charged to uncharged tRNAs
        B=R.*(1-1./D); % actively translating ribosomes - INCLUDING Q
    
        % FIND GROWTH RATE
        ls(i)=sim.form.l(par,e,B)-par('d_t7').*x_het(2).*par('n_t7')./par('M');
    end
end