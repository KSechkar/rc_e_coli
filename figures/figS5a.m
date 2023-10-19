%% figS5a.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% FIGURE S5a

% Performance of the controller with MCMC parameter samples

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% LOAD the MCMC chains

% define number of samples
num_samples=100;

load DREAM_fitting_outcome.mat
ParSet = genparset(chain); DREAMPar.N = size(chain,3);
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : DREAMPar.d ); % take the last 25% of the posterior samples
N_Pars = size(Pars,1); % get number of posterior samples

%% SET UP the simulator

sim=cell_simulator;

sim.init_conditions('s')=0.5;

sim=sim.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_base_level')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('inducer_final_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('step_time')=72; % disturbance flicks transcription reg. func. from 0 to 1 at t=72
sim.ext.input_func_parameters('slope_duration')=0.1;% disturbance flicks transcription reg. func. from 0 to 1 at t=72

sim.het.parameters('c_dist')=100; % gene copy number
sim.het.parameters('a_dist')=500; % max. gene transcription rate

% no output protein expression here
sim.het.parameters('c_x')=0; % gene copy number
sim.het.parameters('a_x')=0; % max. gene transcription rate

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

sim.het.parameters('c_amp')=100; % integral controller amplifier gene copy number
sim.het.parameters('a_amp')=4000; % integral controller amplifier transcription rate
   
% push amended parameter values
sim=sim.push_het();

% simulation parameters
sim.tf =  200;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% SIMULATE with different parameter values

% initialise time and trajectory vector arrays - closed loop
rel_ts={};
rel_psens_trajs={};
rel_l_trajs={};
rel_D_trajs={};
% initialise time and trajectory vector arrays - open loop
rel_ts_openloop={};
rel_psens_trajs_openloop={};
rel_l_trajs_openloop={};
rel_D_trajs_openloop={};

% tried parameter values
tried_thetas=zeros(num_samples,4);

for sample_cntr=1:num_samples
    % randomly draw parameter values using CDFs
    random_draw=randi(N_Pars);
    tried_thetas(sample_cntr,:)=Pars(random_draw,:);
    disp(['Sample ',num2str(sample_cntr),': Theta=',num2str(tried_thetas(sample_cntr,:))]) % print resultant sum of squared errors
    
    % set the parametr values accordingly
    sim.parameters('a_r') = sim.parameters('a_a').*tried_thetas(sample_cntr,1); % ribosome transcription rate (/h) - rescaled!
    sim.parameters('nu_max') = tried_thetas(sample_cntr,2); % max metabolic rate (/h)
    sim.parameters('K_e') = tried_thetas(sample_cntr,3); % elongation rate Hill constant (nM)
    sim.parameters('K_nut') = tried_thetas(sample_cntr,3); % tRNA charging rate Hill constant (nM)
    sim.parameters('kcm') = tried_thetas(sample_cntr,4); % chloramphenical binding rate constant (/h/nM)

    sim = sim.simulate_model;
    
    % GET relevant time frame (start 'inter_rad' h before disturbance, end 'inter_rad' h after)
    % note: plus one point either side of the interval to have prettier
    % plots
    inter_rad=7.5;
    
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
    
    
    % CALCULATE annihilator transcription regulation func., growth rate, RC denominator
    calculated=calc(sim,rel_dist,rel_t);

    % RECORD the times and trajectories
    rel_ts{sample_cntr}=rel_t;
    rel_psens_trajs{sample_cntr}=rel_dist(:,9+sim.num_het+1)./predist_psens;
    rel_l_trajs{sample_cntr}=calculated.ls./predist_l;
    rel_D_trajs{sample_cntr}=calculated.Ds/predist_D;

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
    sim_openloop.het.parameters('a_sens')=sim.het.parameters('a_sens');
    
    % no controller or output protein expression here
    sim_openloop.het.parameters('c_x')=0; % gene copy number
    sim_openloop.het.parameters('c_act')=0; % gene copy number
    sim_openloop.het.parameters('c_anti')=0; % gene copy number
    sim_openloop.het.parameters('c_amp')=0; % gene copy number
      
    % push amended parameter values
    sim_openloop=sim_openloop.push_het();

    % set the sampled parameter values
    sim_openloop.parameters('a_r') = sim_openloop.parameters('a_a').*tried_thetas(sample_cntr,1); % ribosome transcription rate (/h) - rescaled!
    sim_openloop.parameters('nu_max') = tried_thetas(sample_cntr,2); % max metabolic rate (/h)
    sim_openloop.parameters('K_e') = tried_thetas(sample_cntr,3); % elongation rate Hill constant (nM)
    sim_openloop.parameters('K_nut') = tried_thetas(sample_cntr,3); % tRNA charging rate Hill constant (nM)
    sim_openloop.parameters('kcm') = tried_thetas(sample_cntr,4); % chloramphenical binding rate constant (/h/nM)
    
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

    % CALCULATE annihilator transcription regulation func., growth rate, RC denominator
    calculated_openloop=calc(sim_openloop,rel_dist_openloop,rel_t_openloop);

    % RECORD the times and trajectories
    rel_ts_openloop{sample_cntr}=rel_t_openloop;
    rel_psens_trajs_openloop{sample_cntr}=rel_dist_openloop(:,9+sim.num_het+1)./predist_psens_openloop;
    rel_l_trajs_openloop{sample_cntr}=calculated_openloop.ls./predist_l_openloop;
    rel_D_trajs_openloop{sample_cntr}=calculated_openloop.Ds/predist_D_openloop;
end

%% FIGURE 6 d - sensor protein conc.

Fd = figure('Position',[0 0 250 195],'Renderer','painters');
set(Fd, 'defaultAxesFontSize', 9)
set(Fd, 'defaultLineLineWidth', 1.25)
hold on

% open loop
for sample_cntr=1:num_samples
    plot(rel_ts_openloop{sample_cntr}-dist_time-2.5,rel_psens_trajs_openloop{sample_cntr},'Color',[0.6350 0.0780 0.1840 0.02])
end

% closed loop
for sample_cntr=1:num_samples
    plot(rel_ts{sample_cntr}-dist_time-2.5,rel_psens_trajs{sample_cntr},'Color',[0 0.4470 0.7410 0.02])
end

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'p_{sens}:p_{sens}^0, relative', 'sensor prot. conc.'},'FontName','Arial')

ylim([0.85 1.1])
xlim([-5 5])
xticks(-15:2.5:15)

grid on
box on
axis square
hold off

%% FIGURE 6 e - growth rate

Fe = figure('Position',[0 0 250 195],'Renderer','painters');
set(Fe, 'defaultAxesFontSize', 9)
set(Fe, 'defaultLineLineWidth', 1.25)

hold on

% open loop
for sample_cntr=1:num_samples
    plot(rel_ts_openloop{sample_cntr}-dist_time-2.5,rel_l_trajs_openloop{sample_cntr},'Color',[0.6350 0.0780 0.1840 0.02])
end

% closed loop
for sample_cntr=1:num_samples
    plot(rel_ts{sample_cntr}-dist_time-2.5,rel_l_trajs{sample_cntr},'Color',[0 0.4470 0.7410 0.02])
end

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'\lambda:\lambda^0, relative growth rate'},'FontName','Arial')

ylim([0.85 1.1])
xlim([-5 5])
xticks(-15:2.5:15)

% legend({'open loop','closed loop'},'FontName','Arial','FontSize',8,'Location','northwest')

grid on
box on
axis square
hold off

%% FIGURE 6 f - resource competition denominator

Ff = figure('Position',[0 0 250 195],'Renderer','painters');
set(Ff, 'defaultAxesFontSize', 9)
set(Ff, 'defaultLineLineWidth', 1.25)

hold on

% open loop
for sample_cntr=1:num_samples
    plot(rel_ts_openloop{sample_cntr}-dist_time-2.5,rel_D_trajs_openloop{sample_cntr},'Color',[0.6350 0.0780 0.1840 0.02])
end

% closed loop
for sample_cntr=1:num_samples
    plot(rel_ts{sample_cntr}-dist_time-2.5,rel_D_trajs{sample_cntr},'Color',[0 0.4470 0.7410 0.02])
end

xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'D:D^0, relative RC','denominator value'},'FontName','Arial')

xlim([-5 5])
ylim([0.85 1.1])
xticks(-15:2.5:15)

% legend({'open loop','closed loop'},'FontName','Arial','FontSize',8,'Location','southwest')

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