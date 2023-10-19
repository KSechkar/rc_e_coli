%% fig5bcdef.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% FIGURE S6: b,c,d,e,f

% Showcasing how our proportional-integral controller mitigates cisturbances 
% caused by extra heterologous mRNA expression and keeps the burden almost
% constant - with su

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP and RUN the simulator

sim=cell_simulator;

sim.init_conditions('s')=0.5;

sim=sim.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules
% sim=sim.load_heterologous_and_external('pi_controller','oscillating_inducer');

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_base_level')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('inducer_final_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('step_time')=72; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('slope_duration')=0.1;% disturbance flicks transcription reg. func. from 0 to 1 at t=72
% sim.ext.input_func_parameters('wave_period')=2.5; % period of the sine wave (h)
% sim.ext.input_func_parameters('wave_amplitude')=1; % amplitude of the sine wave
% sim.ext.input_func_parameters('oscillation_start_time')=72; % start time of oscillations

sim.het.parameters('c_dist')=1000; % gene copy number
sim.het.parameters('a_dist')=500; % max. gene transcription rate

% no output protein expression here
sim.het.parameters('c_x')=0; % gene copy number
sim.het.parameters('a_x')=0; % max. gene transcription rate

% integral controller parameters
sim.het.parameters('K_dna(anti)-sens')=7000; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(anti)-sens')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('K_dna(amp)-act')=700; % sensor prot.-DNA binding Hill constant
sim.het.parameters('eta_dna(amp)-act')=1; % sensor prot.-DNA binding Hill coefficient

sim.het.parameters('kb_anti')=3000; % atcuator-annihilator binding rate constant
sim.het.parameters('c_sens')=100;
sim.het.parameters('a_sens')=50; % sensor gene transcription rate
sim.het.parameters('a_anti')=800; % annigilator transcription rate
sim.het.parameters('a_act')=400; % actuator transcription rate

sim.het.parameters('a_amp')=100; % integral controller amplifier gene copy number
sim.het.parameters('a_amp')=4000; % integral controller amplifier transcription rate
   
% push amended parameter values
sim=sim.push_het();

% simulate
sim.tf =  200;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
sim = sim.simulate_model;


%% GET relevant time frame (start 'inter_rad' h before disturbance, end 'inter_rad' h after)
% note: plus one point either side of the interval to have prettier plots

inter_rad=15;

dist_time=sim.ext.input_func_parameters('step_time'); % time of disturbance
% dist_time=sim.ext.input_func_parameters('oscillation_start_time');

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

%% GET l and D before disturbance in order to plot relative values

% find last point before disturbance
for i=1:size(sim.t,1)
    if(sim.t(i)>=dist_time)
        predist_pt=i-1;
        break
    end
end
predist_x=sim.x(predist_pt,:); % x just before the disturbance
predist_psens=predist_x(9+sim.num_het+1); % sens is the first gene in the list
[predist_l, predist_D] = get_predist(sim,predist_x); % get the pre-disturbance l and D values

%% SIMULATE THE SYSTEM WITHOUT THE CONTROLLER

sim_openloop=cell_simulator;

sim_openloop.init_conditions('s')=sim.init_conditions('s');

sim_openloop=sim_openloop.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules
% sim_openloop=sim_openloop.load_heterologous_and_external('pi_controller','oscillating_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim_openloop.ext.input_func_parameters('inducer_base_level')=sim.ext.input_func_parameters('inducer_base_level'); 
sim_openloop.ext.input_func_parameters('inducer_final_level')=sim.ext.input_func_parameters('inducer_final_level'); 
sim_openloop.ext.input_func_parameters('step_time')=sim.ext.input_func_parameters('step_time'); 
sim_openloop.ext.input_func_parameters('slope_duration')=sim.ext.input_func_parameters('slope_duration');

% sim_openloop.ext.input_func_parameters('wave_period')=sim.ext.input_func_parameters('wave_period');
% sim_openloop.ext.input_func_parameters('wave_amplitude')=sim.ext.input_func_parameters('wave_amplitude');
% sim_openloop.ext.input_func_parameters('oscillation_start_time')=sim.ext.input_func_parameters('oscillation_start_time');

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

% simulate
sim_openloop.tf=sim.tf;
sim_openloop = sim_openloop.simulate_model;

% find first point in the frame
for i=1:size(sim_openloop.t,1)
    if(sim_openloop.t(i)>=dist_time-inter_rad)
        first_pt_openloop=i-1;
        
        % correction if first point already relevant
        if(first_pt_openloop==0)
            first_pt_openloop=1;
        end

        break
    end
end

% find last point in the frame
for i=size(sim_openloop.t,1):(-1):1
    if(sim_openloop.t(i)<=dist_time+inter_rad)
        last_pt_openloop=i+1;

        % correction if last point still relevant
        if(last_pt_openloop>size(sim_openloop.t,1))
            last_pt_openloop=i;
        end

        break
    end
end

% record relevant time points
rel_t_openloop=sim_openloop.t(first_pt_openloop:last_pt_openloop);

% record states of the cell at these times
rel_dist_openloop=sim_openloop.x(first_pt_openloop:last_pt_openloop,:);

% find last point before disturbance and the respective l and D values
for i=1:size(sim_openloop.t,1)
    if(sim_openloop.t(i)>=dist_time)
        predist_pt_openloop=i-1;
        break
    end
end
predist_x_openloop=sim_openloop.x(predist_pt_openloop,:); % x just before the disturbance
predist_psens_openloop=predist_x_openloop(9+sim_openloop.num_het+1); % sens is the first gene in the list
[predist_l_openloop, predist_D_openloop] = get_predist(sim_openloop,predist_x_openloop); % get the pre-disturbance l and D values

%% SIMULATE THE SYSTEM WITHOUT PROPORTIONAL CONTROL

sim_noprop=cell_simulator_no_proportional;

sim_noprop.init_conditions('s')=sim.init_conditions('s');

sim_noprop=sim_noprop.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules
% sim_noprop=sim_noprop.load_heterologous_and_external('pi_controller','oscillating_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim_noprop.ext.input_func_parameters('inducer_base_level')=sim.ext.input_func_parameters('inducer_base_level'); 
sim_noprop.ext.input_func_parameters('inducer_final_level')=sim.ext.input_func_parameters('inducer_final_level'); 
sim_noprop.ext.input_func_parameters('step_time')=sim.ext.input_func_parameters('step_time'); 
sim_noprop.ext.input_func_parameters('slope_duration')=sim.ext.input_func_parameters('slope_duration');

% sim_noprop.ext.input_func_parameters('wave_period')=sim.ext.input_func_parameters('wave_period');
% sim_noprop.ext.input_func_parameters('wave_amplitude')=sim.ext.input_func_parameters('wave_amplitude');
% sim_noprop.ext.input_func_parameters('oscillation_start_time')=sim.ext.input_func_parameters('oscillation_start_time');

sim_noprop.het.parameters('c_dist')=sim.het.parameters('c_dist'); % gene copy number
sim_noprop.het.parameters('a_dist')=sim.het.parameters('a_dist'); % max. gene transcription rate

% no output protein expression here
sim_noprop.het.parameters('c_x')=sim.het.parameters('c_x'); % gene copy number
sim_noprop.het.parameters('a_x')=sim.het.parameters('a_x'); % max. gene transcription rate

% integral controller parameters
sim_noprop.het.parameters('K_dna(anti)-sens')=sim.het.parameters('K_dna(anti)-sens'); % sensor prot.-DNA binding Hill constant
sim_noprop.het.parameters('eta_dna(anti)-sens')=sim.het.parameters('eta_dna(anti)-sens'); % sensor prot.-DNA binding Hill coefficient

sim_noprop.het.parameters('K_dna(amp)-act')=sim.het.parameters('K_dna(amp)-act'); % sensor prot.-DNA binding Hill constant
sim_noprop.het.parameters('eta_dna(amp)-act')=sim.het.parameters('eta_dna(amp)-act'); % sensor prot.-DNA binding Hill coefficient

sim_noprop.het.parameters('kb_anti')=sim.het.parameters('kb_anti'); % atcuator-annihilator binding rate constant
sim_noprop.het.parameters('c_sens')=sim.het.parameters('c_sens');
sim_noprop.het.parameters('a_sens')=sim.het.parameters('a_sens'); % sensor gene transcription rate
sim_noprop.het.parameters('a_anti')=sim.het.parameters('a_anti'); % annigilator transcription rate
sim_noprop.het.parameters('a_act')=sim.het.parameters('a_act'); % actuator transcription rate

sim_noprop.het.parameters('c_amp')=sim.het.parameters('c_amp'); % integral controller amplifier gene copy number
sim_noprop.het.parameters('a_amp')=sim.het.parameters('a_amp'); % integral controller amplifier transcription rate
  
% push amended parameter values
sim_noprop=sim_noprop.push_het();

% simulate
sim_noprop.tf=sim.tf;
sim_noprop = sim_noprop.simulate_model;

% find first point in the frame
for i=1:size(sim_noprop.t,1)
    if(sim_noprop.t(i)>=dist_time-inter_rad)
        first_pt_noprop=i-1;
        
        % correction if first point already relevant
        if(first_pt_noprop==0)
            first_pt_noprop=1;
        end

        break
    end
end

% find last point in the frame
for i=size(sim_noprop.t,1):(-1):1
    if(sim_noprop.t(i)<=dist_time+inter_rad)
        last_pt_noprop=i+1;

        % correction if last point still relevant
        if(last_pt_noprop>size(sim_noprop.t,1))
            last_pt_noprop=i;
        end

        break
    end
end

% record relevant time points
rel_t_noprop=sim_noprop.t(first_pt_noprop:last_pt_noprop);

% record states of the cell at these times
rel_dist_noprop=sim_noprop.x(first_pt_noprop:last_pt_noprop,:);

% find last point before disturbance and the respective l and D values
for i=1:size(sim_noprop.t,1)
    if(sim_noprop.t(i)>=dist_time)
        predist_pt_noprop=i-1;
        break
    end
end
predist_x_noprop=sim_noprop.x(predist_pt_noprop,:); % x just before the disturbance
predist_psens_noprop=predist_x_noprop(9+sim_noprop.num_het+1); % sens is the first gene in the list
[predist_l_noprop, predist_D_noprop] = get_predist(sim_noprop,predist_x_noprop); % get the pre-disturbance l and D values

%% ANALYTICALLY ESTIMATE lambda and D

% finding no-burden values of translation rate, dissociation constants,
% rib. gene transc. regulation function
sim_nb=cell_simulator;
sim_nb.init_conditions('s')=0.07;
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

    

%% FIGURE 6 b

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
ylim([0 6e4])

grid on
box on
axis square
hold off

%% FIGURE 6 c

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
ylim([0 800])
yticks(0:200:800)

grid on
box on
axis square
hold off


%% CALCULATE annihilator transcription regulation func., growth rate, RC denominator

calculated=calc(sim,rel_dist,rel_t);
Fs=calculated.Fs;
ls=calculated.ls;
Ds=calculated.Ds;

%% CALCULATE, for open loop, growth rate, RC denominator

calculated_openloop=calc(sim_openloop,rel_dist_openloop,rel_t_openloop);
ls_openloop=calculated_openloop.ls;
Ds_openloop=calculated_openloop.Ds;

%% CALCULATE, for no proprtional feedback, growth rate, RC denominator

calculated_noprop=calc(sim_noprop,rel_dist_noprop,rel_t_noprop);
ls_noprop=calculated_noprop.ls;
Ds_noprop=calculated_noprop.Ds;


%% FIGURE 6 d - sensor protein conc.

Fd = figure('Position',[0 0 250 186]);
set(Fd, 'defaultAxesFontSize', 9)
set(Fd, 'defaultLineLineWidth', 1.25)
hold on

% plot model predictions
plot(rel_t-dist_time,rel_dist(:,9+sim.num_het+1)/predist_psens,'Color',[0 0.4470 0.7410])
%plot(rel_t_openloop-dist_time,rel_dist_openloop(:,9+sim.num_het+1)/predist_psens_openloop,'Color',[0.6350 0.0780 0.1840])
plot(rel_t_noprop-dist_time,rel_dist_noprop(:,9+sim.num_het+1)/predist_psens_noprop,'Color',[0.8500 0.3250 0.0980])

% plot ideal value
plot([-inter_rad inter_rad],[1 1],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'p_{sens}:p_{sens}^0, relative', 'sensor prot. conc.'},'FontName','Arial')
legend({'w/ prop. fbck','no prop. fbck'},'FontName','Arial','FontSize',8,'Location','northwest')

ylim([0.95 1.01])
xlim([-inter_rad inter_rad])
xticks(-inter_rad:inter_rad/2:inter_rad)

grid on
box on
axis square
hold off

%% FIGURE 6 e - growth rate

Fe = figure('Position',[0 0 250 186]);
set(Fe, 'defaultAxesFontSize', 9)
set(Fe, 'defaultLineLineWidth', 1.25)

hold on

% plot model predictions
plot(rel_t-dist_time,ls/predist_l,'Color',[0 0.4470 0.7410])

% results with no proportional feedback
plot(rel_t_noprop-dist_time,ls_noprop/predist_l_noprop,'Color',[0.8500 0.3250 0.0980])

% plot analytically calculated target value
plot([-inter_rad inter_rad],[lambda_estimated/predist_l lambda_estimated/predist_l],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'\lambda:\lambda^0, relative growth rate'},'FontName','Arial')

ylim([0.95 1.01])
xlim([-inter_rad inter_rad])
xticks(-inter_rad:inter_rad/2:inter_rad)

legend({'w/ prop. fbck','no prop. fbck'},'FontName','Arial','FontSize',8,'Location','northwest')

grid on
box on
axis square
hold off

%% FIGURE 6 f - resource competition denominator

Ff = figure('Position',[0 0 250 189]);
set(Ff, 'defaultAxesFontSize', 9)
set(Ff, 'defaultLineLineWidth', 1.25)

hold on

% plot model predictions
plot(rel_t-dist_time,Ds/predist_D,'Color',[0 0.4470 0.7410])

% results with no proportional feedback
plot(rel_t_noprop-dist_time,Ds_noprop/predist_D_noprop,'Color',[0.8500 0.3250 0.0980])

% plot analytically calculated target value
plot([-inter_rad inter_rad],[D_estimated/predist_D D_estimated/predist_D],'k:') 

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'D:D^0, relative RC','denominator value'},'FontName','Arial')

xlim([-inter_rad inter_rad])
ylim([0.98 1.04])
xticks(-inter_rad:inter_rad/2:inter_rad)

legend({'w/ prop. fbck','no prop. fbck'},'FontName','Arial','FontSize',8,'Location','northwest')

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

%% PRE-DISTURBANCE GROWTH RATE AND D
%

% calculated from the state of the system right before disturbance
function [l,D]=get_predist(sim,predist_x) 
    par=sim.parameters;
    % STATE VECTOR TO SINGLE VARIABLES
    m_a = predist_x(1);
    m_r = predist_x(2);
    p_a = predist_x(3);
    R = predist_x(4);
    tc = predist_x(5);
    tu = predist_x(6);
    Bcm = predist_x(7);
    s = predist_x(8);
    h = predist_x(9);
    x_het = predist_x(10 : (9+2*sim.num_het) );

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