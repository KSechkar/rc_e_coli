%% figS5a.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% FIGURE S5a

% Performance of the controller with MCMC parameter samples

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% DEFINE number of trajectories and how we simulate the trajectories

num_trajs=16;

%% LOAD prematurely terminated stochastic trajectory arrays and simulators

load('figS5b_controller_16-10.mat')
load('figS5b_openloop_16-10.mat')
load('sims.mat')

%% RUN parallel simulations
tic
for traj_cntr=1:num_trajs
    % adjust disturbance time based on how much of the trajectory has been simulated
    sim{traj_cntr}.ext.input_func_parameters('step_time')=sim{traj_cntr}.ext.input_func_parameters('step_time')-rel_ts{traj_cntr}(end);
    sim_openloop{traj_cntr}.ext.input_func_parameters('step_time')=sim_openloop{traj_cntr}.ext.input_func_parameters('step_time')-rel_ts_openloop{traj_cntr}(end);
    sim{traj_cntr}=sim{traj_cntr}.push_het;

    %% CLOSED LOOP
    % STOCHASTIC simulation parameters
    sim{traj_cntr}.tf=15-rel_ts{traj_cntr}(end); % simulation time frame

    
    % SET the initial condition based on tht trajector's last point
    traj_endpoint=sim{traj_cntr}.x(end,:);
    
    % mRNAs
    sim{traj_cntr}.init_conditions('m_a')=traj_endpoint(1);
    sim{traj_cntr}.init_conditions('m_r')=traj_endpoint(2);
    % proteins
    sim{traj_cntr}.init_conditions('p_a')=traj_endpoint(3);
    sim{traj_cntr}.init_conditions('R')=traj_endpoint(4);
    % tRNAs
    sim{traj_cntr}.init_conditions('tc')=traj_endpoint(5);
    sim{traj_cntr}.init_conditions('tu')=traj_endpoint(6);
    % inactivated ribosomes
    sim{traj_cntr}.init_conditions('Bcm')=traj_endpoint(7);
    % nutrient quality and chloramphenicol
    sim{traj_cntr}.init_conditions('s')=traj_endpoint(8);
    sim{traj_cntr}.init_conditions('h')=traj_endpoint(9);
    
    % heterologous
    x_het=traj_endpoint(10 : (9+2*sim{traj_cntr}.num_het+sim{traj_cntr}.num_misc) );
    for j=1:sim{traj_cntr}.num_het
        % mRNA
        sim{traj_cntr}.init_conditions(['m_',sim{traj_cntr}.het.names{j}])=x_het(j);
        % protein
        sim{traj_cntr}.init_conditions(['p_',sim{traj_cntr}.het.names{j}])=x_het(sim{traj_cntr}.num_het+j);
    end
    % miscellaneous
    for j=1:sim{traj_cntr}.num_misc
        sim{traj_cntr}.init_conditions([sim{traj_cntr}.het.misc_names{j}])=x_het(2*sim{traj_cntr}.num_het+j);
    end

    % find deterministic steady-statte p_sens, l and D
    steady_psens{traj_cntr}=sim{traj_cntr}.x(1,9+sim{traj_cntr}.num_het+1);
    calculated_detss=calc(sim{traj_cntr},sim{traj_cntr}.x,sim{traj_cntr}.t);
    steady_l{traj_cntr}=calculated_detss.ls(1);
    steady_D{traj_cntr}=calculated_detss.Ds(1);
    
    % SIMULATE stochastic trajectories
    sim{traj_cntr}=sim{traj_cntr}.generate_het_stoichiometry_matrix();
    sim{traj_cntr} = sim{traj_cntr}.simulate_model_hybrid_tauleap;
    
    % calculate l and D
    calculated=calc(sim{traj_cntr},sim{traj_cntr}.x,sim{traj_cntr}.t);

    % record trajectories
    rel_ts{traj_cntr}=[rel_ts{traj_cntr}; sim{traj_cntr}.t(2:end)+rel_ts{traj_cntr}(end)];
    rel_psens_trajs{traj_cntr}=[rel_psens_trajs{traj_cntr}; sim{traj_cntr}.x(2:end,9+sim{traj_cntr}.num_het+1)./steady_psens{traj_cntr}];
    rel_l_trajs{traj_cntr}=[rel_l_trajs{traj_cntr}; (calculated.ls(2:end))./steady_l{traj_cntr}];
    rel_D_trajs{traj_cntr}=[rel_D_trajs{traj_cntr}; (calculated.Ds(2:end))./steady_D{traj_cntr}];

    %% OPEN LOOP
    sim_openloop{traj_cntr}.ext.input_func_parameters('step_time')=sim_openloop{traj_cntr}.ext.input_func_parameters('step_time')-rel_ts_openloop{traj_cntr}(end);
    sim_openloop{traj_cntr}.tf=15-rel_ts_openloop{traj_cntr}(end); % simulation time frame
    sim_openloop{traj_cntr}=sim_openloop{traj_cntr}.push_het;
    
    % SET the initial condition according to the deterministic steady state
    % get steady state
    traj_endpoint_openloop = sim_openloop{traj_cntr}.x(end,:);

    % mRNAs
    sim_openloop{traj_cntr}.init_conditions('m_a')=traj_endpoint_openloop(1);
    sim_openloop{traj_cntr}.init_conditions('m_r')=traj_endpoint_openloop(2);
    % proteins
    sim_openloop{traj_cntr}.init_conditions('p_a')=traj_endpoint_openloop(3);
    sim_openloop{traj_cntr}.init_conditions('R')=traj_endpoint_openloop(4);
    % tRNAs
    sim_openloop{traj_cntr}.init_conditions('tc')=traj_endpoint_openloop(5);
    sim_openloop{traj_cntr}.init_conditions('tu')=traj_endpoint_openloop(6);
    % inactivated ribosomes
    sim_openloop{traj_cntr}.init_conditions('Bcm')=traj_endpoint_openloop(7);
    % nutrient quality and chloramphenicol
    sim_openloop{traj_cntr}.init_conditions('s')=traj_endpoint_openloop(8);
    sim_openloop{traj_cntr}.init_conditions('h')=traj_endpoint_openloop(9);
    
    % heterologous
    x_het=traj_endpoint_openloop(10 : (9+2*sim_openloop{traj_cntr}.num_het+sim_openloop{traj_cntr}.num_misc) );
    for j=1:sim_openloop{traj_cntr}.num_het
        % mRNA
        sim_openloop{traj_cntr}.init_conditions(['m_',sim_openloop{traj_cntr}.het.names{j}])=x_het(j);
        % protein
        sim_openloop{traj_cntr}.init_conditions(['p_',sim_openloop{traj_cntr}.het.names{j}])=x_het(sim_openloop{traj_cntr}.num_het+j);
    end
    % miscellaneous
    for j=1:sim_openloop{traj_cntr}.num_misc
        sim_openloop{traj_cntr}.init_conditions([sim_openloop{traj_cntr}.het.misc_names{j}])=x_het(2*sim_openloop{traj_cntr}.num_het+j);
    end

    % find deterministic steady-statte l and D
    steady_psens_openloop{traj_cntr}=sim_openloop{traj_cntr}.x(1,9+sim{traj_cntr}.num_het+1);
    calculated_detss_openloop=calc(sim_openloop{traj_cntr},sim_openloop{traj_cntr}.x,sim_openloop{traj_cntr}.t);
    steady_l_openloop{traj_cntr}=calculated_detss_openloop.ls(1);
    steady_D_openloop{traj_cntr}=calculated_detss_openloop.Ds(1);
    
    % SIMULATE stochastic trajectories
    sim_openloop{traj_cntr} = sim_openloop{traj_cntr}.simulate_model_hybrid_tauleap;
    
    dist_time=sim_openloop{traj_cntr}.ext.input_func_parameters('step_time'); % time of disturbance
    
    % calculate l and D
    calculated_openloop=calc(sim_openloop{traj_cntr},sim_openloop{traj_cntr}.x,sim_openloop{traj_cntr}.t);

    % record trajectories
    rel_ts_openloop{traj_cntr}=[rel_ts_openloop{traj_cntr}; sim_openloop{traj_cntr}.t(2:end)+rel_ts_openloop{traj_cntr}(end)];
    rel_psens_trajs_openloop{traj_cntr}=[rel_psens_trajs_openloop{traj_cntr}; sim_openloop{traj_cntr}.x(2:end,9+sim_openloop{traj_cntr}.num_het+1)./steady_psens_openloop{traj_cntr}];
    rel_l_trajs_openloop{traj_cntr}=[rel_l_trajs_openloop{traj_cntr}; (calculated_openloop.ls(2:end))./steady_l_openloop{traj_cntr}];
    rel_D_trajs_openloop{traj_cntr}=[rel_D_trajs_openloop{traj_cntr}; (calculated_openloop.Ds(2:end))./steady_D_openloop{traj_cntr}];
end
toc
% save all trjectories for future use
save('figS5b_controller.mat', 'rel_ts', 'rel_psens_trajs', 'rel_l_trajs', 'rel_D_trajs')
save('figS5b_openloop.mat', 'rel_ts_openloop', 'rel_psens_trajs_openloop', 'rel_l_trajs_openloop', 'rel_D_trajs_openloop')



%% FIGURE 6 d - sensor protein conc.

Fd = figure('Position',[0 0 350 260]);
set(Fd, 'defaultAxesFontSize', 9)
set(Fd, 'defaultLineLineWidth', 1.25)
hold on

% open loop plots
for traj_cntr=1:num_trajs
    plot(rel_ts_openloop{traj_cntr},rel_psens_trajs_openloop{traj_cntr},'Color',[0.6350 0.0780 0.1840 0.05])
end
% average trajectory
plot(rel_ts{1},mean(cat(2,rel_psens_trajs_openloop{1}),2),'Color',[0.6350 0.0780 0.1840 1]);

% closed loop plots
psens_ends=zeros(1,num_trajs);
for traj_cntr=1:num_trajs
    plot(rel_ts{traj_cntr},rel_psens_trajs{traj_cntr},'Color',[0 0.4470 0.7410 0.05])
end
% average trajectory
plot(rel_ts{1},mean(cat(2,rel_psens_trajs{1}),2),'Color',[0 0.4470 0.7410 1]);

% xlabel('t, time since disturbance [h]','FontName','Arial')
ylabel({'p_{sens}:p_{sens}^0, relative', 'sensor prot. conc.'},'FontName','Arial')

% ylim([0.95 1.05])
% xlim([-inter_rad inter_rad])
% xticks(-inter_rad:inter_rad/2:inter_rad)

grid on
box on
axis square
hold off

%% FIGURE 6 e - growth rate
% 
% Fe = figure('Position',[0 0 250 186]);
% set(Fe, 'defaultAxesFontSize', 9)
% set(Fe, 'defaultLineLineWidth', 1.25)
% 
% hold on
% 
% % open loop plots
% for traj_cntr=1:num_trajs
%     plot(rel_ts_openloop{traj_cntr}-dist_time,rel_l_trajs_openloop{traj_cntr},'Color',[0.6350 0.0780 0.1840 0.5])
% end
% 
% % closed loop plots
% for traj_cntr=1:num_trajs
%     plot(rel_ts{traj_cntr}-dist_time,rel_l_trajs{traj_cntr},'Color',[0 0.4470 0.7410 0.05])
% end
% 
% xlabel('t, time since disturbance [h]','FontName','Arial')
% ylabel({'\lambda:\lambda^0, relative growth rate'},'FontName','Arial')
% 
% % ylim([0.95 1.05])
% % xlim([-inter_rad inter_rad])
% % xticks(-inter_rad:inter_rad/2:inter_rad)
% % 
% % legend({'open loop','closed loop'},'FontName','Arial','FontSize',8,'Location','northwest')
% 
% grid on
% box on
% axis square
% hold off
% 
% %% FIGURE 6 f - resource competition denominator
% 
% Ff = figure('Position',[0 0 250 189]);
% set(Ff, 'defaultAxesFontSize', 9)
% set(Ff, 'defaultLineLineWidth', 1.25)
% 
% hold on
% 
% % open loop plots
% for traj_cntr=1:num_trajs
%     plot(rel_ts_openloop{traj_cntr}-dist_time,rel_D_trajs_openloop{traj_cntr},'Color',[0.6350 0.0780 0.1840 0.05])
% end
% 
% % closed loop plots
% for traj_cntr=1:num_trajs
%     plot(rel_ts{traj_cntr}-dist_time,rel_D_trajs{traj_cntr},'Color',[0 0.4470 0.7410 0.05])
% end
% 
% xlabel('t, time since disturbance [h]','FontName','Arial')
% ylabel({'D:D^0, relative RC','denominator value'},'FontName','Arial')
% 
% % xlim([-inter_rad inter_rad])
% % ylim([0.95 1.05])
% % xticks(-inter_rad:inter_rad/2:inter_rad)
% 
% % legend({'no controller','w/ controller'},'FontName','Arial','FontSize',8,'Location','southwest')
% 
% grid on
% box on
% axis square
% hold off

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