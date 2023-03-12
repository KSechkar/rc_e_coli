%% fig5bcdef.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 5: b,c,d,e,f

% Showcasing how our proportional-integral controller mitigates cisturbances 
% caused by extra heterologous mRNA expression and keeps the burden almost constant 

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP and RUN the simulator

% 2 different setups - ideal AIF and realistic scenario
sim=cell_simulator;

sim=sim.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_base_level')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('inducer_final_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('step_time')=30; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('slope_duration')=0.1;% disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.het.parameters('c_dist')=100; % gene copy number
sim.het.parameters('a_dist')=300; % max. gene transcription rate

% no output protein expression here
sim.het.parameters('c_x')=0; % gene copy number
sim.het.parameters('a_x')=0; % max. gene transcription rate

% integral controller parameters
sim.het.parameters('K_dna(anti)-sens')=4000; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(anti)-sens')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('K_dna(amp)-act')=4000; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(amp)-act')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('kb_anti')=300; % atcuator-annihilator binding rate constant
sim.het.parameters('a_sens')=22.5; % sensor gene transcription rate
sim.het.parameters('a_anti')=800; % annigilator transcription rate
sim.het.parameters('a_act')=400; % actuator transcription rate

sim.het.parameters('a_amp')=4000; % integral controller amplifier transcription rate
   
% push amended parameter values
sim=sim.push_het();

% simulate
sim.tf =  60;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
sim = sim.simulate_model;


%% GET relevant time frame (start 'inter_rad' h before disturbance, end 'inter_rad' h after)
% note: plus one point either side of the interval to have prettier plots

inter_rad=15;

dist_time=sim.ext.input_func_parameters('step_time'); % time of disturbance

% find first point in the frame
for i=1:size(sim.t,1)
    if(sim.t(i)>=dist_time-inter_rad)
        first_pt=i-1;
        
        % correction if first point already relevant
        if(first_pt==0)
            first_pt=1;
        end

        break
    end
end

% find last point in the frame
for i=size(sim.t,1):(-1):1
    if(sim.t(i)<=dist_time+inter_rad)
        last_pt=i+1;

        % correction if last point still relevant
        if(last_pt>size(sim.t,1))
            last_pt=i;
        end

        break
    end
end

% record relevant time points
rel_t=sim.t(first_pt:last_pt);

% record states of the cell at these times
rel_dist=sim.x(first_pt:last_pt,:);


%% ANALYTICALLY ESTIMATE lambda and D

% finding no-burden values of translation rate, dissociation constants,
% rib. gene transc. regulation function
sim_nb=cell_simulator;
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
[e_nb,Fr_nb,k_a_nb,k_r_nb,k_sens_nb,k_act_nb,k_amp_nb]=get_nb(sim_nb,sim_nb.x(end,:));


% calculate values of 'meaningful paramters'
par=sim.parameters; % IMPORTANT! parameters of the system where the controller's genes ARE expressed, not 'no burden' one
u=par('a_act')./par('a_anti'); % ideal value of F_anti
zeta=par('c_sens').*par('a_sens');

% calculating the value of lambda
lambda_estimated = ...
    e_nb./par('M') .* ...
    Fr_nb.*par('c_r').*par('a_r') .* ...
    par('n_sens').*k_sens_nb./(par('n_r').*k_r_nb) .* ...
    par('K_dna(anti)-sens')./(zeta) .* ...
    (1-u)./u;

% calculating the value of D
D_estimated = ...
    1 + ...
    (lambda_estimated.*zeta) ./ ...
    ((lambda_estimated+par('b_sens')).*k_sens_nb) ./ ...
    (par('n_sens')./par('M') .* par('K_dna(anti)-sens') .* ((1-u)./u) );

    

%% FIGURE 5 b

Fb = figure('Position',[0 0 250 186]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

% colours for plots
colours={[1 0.586 0 1],... % sens
        [0 0.781 0.781 1],... % anti
        [0.781 0 0.781 1],... % act
        [1 0 0.586 1],... % amp
        [0 0 0 0.3] % x
        };

hold on
for i=4:5
    plot(rel_t-dist_time,rel_dist(:,9+i),'Color',colours{i});
end


legend('m_{amp}','m_{dist}',...
    'Location','west','FontName','Arial')

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel('m_i, mRNA concentration [nM]','FontName','Arial')

xlim([(dist_time-inter_rad) (dist_time+inter_rad)])
xticks(-inter_rad:inter_rad/2:inter_rad)
xlim([-inter_rad inter_rad])

grid on
box on
axis square
hold off

%% FIGURE 5 c

Fc = figure('Position',[0 0 250 185]);
set(Fc, 'defaultAxesFontSize', 9)
set(Fc, 'defaultLineLineWidth', 1.25)

% colours for plots
colours={[1 0.586 0 1],... % sens
        [0 0.781 0.781 1],... % anti
        [0.781 0 0.781 1],... % act
        [1 0 0.586 1],... % amp
        [0.3 0.3 0.3 1] % x
        };

hold on
for i=1:3
    plot(rel_t-dist_time,rel_dist(:,9+i),'Color',colours{i});
end


legend('m_{sens}','m_{anti}','m_{act}',...
    'Location','northeast','FontName','Arial')

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel('m_i, mRNA concentration [nM]','FontName','Arial')

xlim([(dist_time-inter_rad) (dist_time+inter_rad)])
xticks(-inter_rad:inter_rad/2:inter_rad)
xlim([-inter_rad inter_rad])
ylim([0 750])
yticks(0:250:750)

grid on
box on
axis square
hold off


%% CALCULATE annihilator transcription regulation func., growth rate, RC denominator

calculated=calc(sim,rel_dist,rel_t);
Fs=calculated.Fs;
ls=calculated.ls;
Ds=calculated.Ds;

%% FIGURE 5 d - control error

Fd = figure('Position',[0 0 250 186]);
set(Fd, 'defaultAxesFontSize', 9)
set(Fd, 'defaultLineLineWidth', 1.25)
hold on

% plot model predictions
plot(rel_t-30,u*ones(size(Fs))-Fs,'Color',[0 0.4470 0.7410])

% plot ideal value
plot([-inter_rad inter_rad],[0 0],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel(', control error','FontName','Arial')

ylim([-0.2 0.2])
xlim([-inter_rad inter_rad])
yticks(-0.2:0.1:0.2)
xticks(-inter_rad:inter_rad/2:inter_rad)

grid on
box on
axis square
hold off

%% FIGURE 5 e - growth rate

Fe = figure('Position',[0 0 250 186]);
set(Fe, 'defaultAxesFontSize', 9)
set(Fe, 'defaultLineLineWidth', 1.25)

hold on

% plot model predictions
plot(rel_t-30,ls,'Color',[0 0.4470 0.7410])

% plot analytically calculated target value
plot([-inter_rad inter_rad],[lambda_estimated lambda_estimated],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel('\lambda, growth rate [1/h]','FontName','Arial')

ylim([0.6 1.8])
yticks(0.6:0.3:1.8)
xlim([-inter_rad inter_rad])
xticks(-inter_rad:inter_rad/2:inter_rad)

grid on
box on
axis square
hold off

%% FIGURE 5 f - resource competition denominator

Ff = figure('Position',[0 0 250 189]);
set(Ff, 'defaultAxesFontSize', 9)
set(Ff, 'defaultLineLineWidth', 1.25)

hold on

% plot model predictions
plot(rel_t-30,Ds,'Color',[0 0.4470 0.7410])

% plot analytically calculated target value
plot([-inter_rad inter_rad],[D_estimated D_estimated],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel('D, RC denominator','FontName','Arial')

ylim([6e4 10e4])
xlim([-inter_rad inter_rad])
xticks(-inter_rad:inter_rad/2:inter_rad)

grid on
box on
axis square
hold off

%% CALCULATOR FUNCTION

% finding annihilator transcription regulation func., growth rate, RC
% denominator on the relevant time interval
function calculated=calc(sim,rel_dist,rel_t)
    calculated.ls=zeros(size(rel_t));
    calculated.Fs=zeros(size(rel_t));   
    calculated.Ds=zeros(size(rel_t));
    
    par=sim.parameters;
    for i=1:size(rel_t)
        % STATE VECTOR TO SINGLE VARIABLES
        m_a = rel_dist(i,1);
        m_r = rel_dist(i,2);
        p_a = rel_dist(i,3);
        R = rel_dist(i,4);
        tc = rel_dist(i,5);
        tu = rel_dist(i,6);
        Bcm = rel_dist(i,7);
        s = rel_dist(i,8);
        h = rel_dist(i,9);
        x_het=rel_dist(i,10 : (9+2*sim.num_het) );
    
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
        
        % RECORD ANNIHILATOR REG.
        calculated.Fs(i)=sim.het.regulation(sim.het.names{2},rel_dist(i,:),0); 
    
        % RECORD GROWTH RATE
        calculated.ls(i)=l;
        
        % RECORD RC DENOMINATOR
        calculated.Ds(i)=D;
    end
end

%% NO BURDEN VALUES FUNCTION
%

function [e,Fr,k_a,k_r,k_sens,k_act,k_amp]=get_nb(sim,ss)   
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
    
end