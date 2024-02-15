% figS5defg.m

% EMERGENT BISTABILITY OF A NON-COOPERATIVE SELF-ACTIVATOR
% Figure S5: d,e,f,g

% Despite promoting its own expression non-cooperatively, a T7 RNAP gene
% controller by a constitutive T7 RNAP promoter may exhibit bistability due
% to the effect of synthetic gene expression on the host cell's growth
% rate. Here, we obtain the bifrucation diagram using CBC

%% CLEAR all variables

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

%% INITIALISE the cell simulator, define default circuit parameters

sim=sim.load_heterologous_and_external('t7_selfact','no_ext'); % load heterologous gene expression module
sim.het.parameters('c_t7')=5;
sim.het.parameters('a_t7')=170;
sim.het.parameters('K_t7')=550;
sim.het.parameters('b_t7')=16.64;
sim.het.parameters('d_t7')=0.18;
sim.het.parameters('baseline_t7')=0.0016; % Liao et al. 2017
sim.het.parameters('t7_toxicity')=0.13; % Liao et al. 2017
sim.het.init_conditions('m_t7')=50;
sim=sim.push_het;

% simulation parameters
sim.tf =  240;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% GET the cell's steady state in a given medium without heterologous gene expression

c_t7=sim.het.parameters('c_t7'); % first, set to zero for t7-rnap-free steady state retrieval
sim.het.parameters('c_t7')=0;
sim=sim.push_het();
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

%% GET the initial state of the system for default t7_toxicity=0
sim.het.parameters('c_t7')=c_t7; % restore c_t7 value
sim.het.parameters('t7_toxicity')=0.13; % Liao et al. 2017
sim.het.init_conditions('m_t7')=1e4;
sim.het.init_conditions('p_t7')=0;
sim=sim.push_het;
sim=sim.simulate_model();
x_init=sim.x(end,:); % initial steady state
pt7_init=x_init(11); % initial steady-state p_t7 concentration
u_init = sim.het.parameters('t7_toxicity'); % initial control input (toxicity value)

%% DEFINE the CBC parameters

% sampling time
cbc_par.sampling_time = 0.05; % (h)
cbc_par.end_time = 240; % (h)

% position of pt7 in the state vector
cbc_par.ref_in_state = 11;

% threshold for declaring x to be in steady state
cbc_par.dxdt_threshold=0.001;

% PID controller gains - just initialised for now
cbc_par.ctrl.Kp=0;
cbc_par.ctrl.Ki=0;
cbc_par.ctrl.Kd=0;

% maximum and minimum  input (toxicity)
cbc_par.ctrl.min_u=0; % (1/nM)
cbc_par.ctrl.max_u=0.15; % (1/nM)

% integral clamping
cbc_par.ctrl.clamp=0.15;%/cbc_par.ctrl.Ki;

% dipslay diagnostics
cbc_par.verbose=true;

% information about error between state and reference
errinfo.errint=0;%u_init/cbc_par.ctrl.Ki; % integral of error
errinfo.preverr=0; % error calculated in the past step (initialised with initial state error)

%% DEFINE the range of p_t7 concentrations to retrieve the equilibrium values for

% pt7_range=(pt7_init+u_init/cbc_par.ctrl.Kp):(-0.01/cbc_par.ctrl.Kp):(pt7_init+0.01/cbc_par.ctrl.Kp);
pt7_range=[...
    3450:(-100):2050,...
    2050:(-100):250,...
    160:(-20):60,...
    35:(-2.5):25,...
    15:(-0.25):7.5,...
    5:(-1):2,...
    1.9:(-0.1):1.6,...
    1.5:(-0.03):1.41];

%% DEFINE feedback gains

% initialise
Kp_vals=zeros(size(pt7_range));

% set gains
for j=1:size(pt7_range,2)
    if(pt7_range(j)>=250)
        Kp_vals(j)=1e-4;
    elseif(pt7_range(j)>=100)
        Kp_vals(j)=5e-4;
    elseif(pt7_range(j)>=60)
        Kp_vals(j)=1e-3;
    elseif(pt7_range(j)>=25)
        Kp_vals(j)=5e-3;
    elseif(pt7_range(j)>=7)
        Kp_vals(j)=2.5e-2;
    else
        Kp_vals(j)=5e-2;
    end
end

%% INITIALISE cbc results storage
% retrieved equilibirium points
cbc_results.pt7=zeros(size(pt7_range)); % p_t7 concentration
cbc_results.tox=zeros(size(pt7_range)); % toxicity (bifurcation parameter)

% overall controlled system trajectories, recorded for every sampling step
x_traj=[x_init];
pt7_traj=[pt7_init];
u_traj=[u_init];
    
%% RUN CBC
% tic
for j = 1:size(pt7_range,2)
    disp(pt7_range(j))

    % set gains
    cbc_par.ctrl.Kp=Kp_vals(j);
    cbc_par.ctrl.Ki=0;
    cbc_par.ctrl.Kd=0;
    
    % for each reference steady state p_t7, call the inner CBC loop
    [STATES, CONTROL, endtime, errinfo]=rc_cbc_inner_loop_own(sim, cbc_par, pt7_range(j), x_traj(end,:), u_traj, errinfo);

    % store the recorded equilibrium point
    cbc_results.pt7(j)=STATES.pt7_end;
    cbc_results.tox(j)=CONTROL.u_end;

    % store the trajectories
    x_traj=[x_traj; STATES.x_traj(1:endtime,:)];
    pt7_traj=[pt7_traj STATES.pt7_traj(1:endtime)];
    u_traj=CONTROL.u_traj;

end
% toc

%% SAVE the outcome
save('figS5defg.mat')

%% FIGURE D - bifurcation diagram

% create figure
Fig_bd = figure('Position',[0 0 320 360]);
set(Fig_bd, 'defaultAxesFontSize', 9)
set(Fig_bd, 'defaultLineLineWidth', 1.25)

hold on

% plot retireived steady states
plot(cbc_results.tox,cbc_results.pt7,'-o',...
    'LineWidth',1,'Color','r','MarkerSize',4,...
    'DisplayName','Retrieved equilbiria')

% plot(u_traj,pt7_traj)

%legend('Location','east')

xlabel('\gamma_{tox}, T7 RNAP toxicity [1/nM]','FontName','Arial')
ylabel('p_{t7}, T7 RNAP concentration [nM]','FontName','Arial')

xlim([0 0.13])
ylim([1 5e3])
set(gca, 'YScale', 'log')

grid 
box on

%% FIGURE E - reference and retrieved p_t7 values

% time dimension of a single CBC loop iteration
time_dim = round(cbc_par.end_time/cbc_par.sampling_time)-1;

% get the sequence of reference values over time
pt7_refs_of_time=zeros(1,size(pt7_traj,2));
for j = 1:size(pt7_range,2)
    pt7_refs_of_time((j-1)*time_dim+1:j*time_dim)=pt7_range(j);
end

% overall time axis
time_axis=cbc_par.sampling_time*((1:size(pt7_traj,2))-1);

% PLOT
% create figure
Fig_pt7 = figure('Position',[0 0 300 160]);
set(Fig_pt7, 'defaultAxesFontSize', 9)
set(Fig_pt7, 'defaultLineLineWidth', 1.25)

hold on

% plot reference values
plot(time_axis,pt7_refs_of_time,...
    '-','LineWidth',1,'Color','k',...
    'DisplayName','Ref. values')

% plot CBC trajectory
plot(time_axis,pt7_traj,...
    '-','LineWidth',1,'Color','b',...
    'DisplayName','CBC trajectory')

legend('Location','northeast')

xlabel('t, overall simulation time [h]','FontName','Arial')
ylabel('p_{t7}, T7 RNAP conc. [nM]','FontName','Arial')

xlim([0 time_axis(end)])
ylim([1 5e3])
set(gca, 'YScale', 'log')

grid 
box on

%% FIGURE F - \gamma_tox values

% PLOT
% create figure
Fig_tox = figure('Position',[0 0 304 160]);
set(Fig_tox, 'defaultAxesFontSize', 9)
set(Fig_tox, 'defaultLineLineWidth', 1.25)

hold on

% plot CBC trajectory
plot(time_axis,u_traj,...
    '-','LineWidth',1,'Color','b',...
    'DisplayName','CBC trajectory')

legend('Location','northeast')

xlabel('t, overall simulation time [h]','FontName','Arial')
ylabel('\gamma_{tox}, T7 RNAP tox. [1/nM]','FontName','Arial')

xlim([0 time_axis(end)])
ylim([0 0.15])

grid 
box on

%% FIGURE G - feedback gains

% get the sequence of Kp values over time
Kp_vals_of_time=zeros(1,size(pt7_traj,2));
for j = 1:size(Kp_vals,2)
    Kp_vals_of_time((j-1)*time_dim+1:j*time_dim)=Kp_vals(j);
end

% PLOT
% create figure
Fig_kp = figure('Position',[0 0 303 160]);
set(Fig_kp, 'defaultAxesFontSize', 9)
set(Fig_kp, 'defaultLineLineWidth', 1.25)

hold on

% plot reference values
plot(time_axis, Kp_vals_of_time,...
    '-','LineWidth',1,'Color','m')

xlabel('t, overall simulation time [h]','FontName','Arial')
ylabel('K_{p}, feedback gain','FontName','Arial')

xlim([0 time_axis(end)])
ylim([5e-5 1e-1])
yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
set(gca, 'YScale', 'log')

grid 
box on

%% CBC LOOP FUNCTION
% input = cell smulator, cbc parameters, pt7 reference, starting state, input array over time, error information
function [STATES, CONTROL, endtime, errinfo] = rc_cbc_inner_loop_own(sim, cbc_par, ref_pt7, x0, u_traj, errinfo) 
    % sampling and end times
    Ts=cbc_par.sampling_time; % get sampling time
    Tend=cbc_par.end_time; % get end time
    time_dim=round(Tend/Ts); % input time dimension

    % position of pt7 in the state vector
    ref_in_state=cbc_par.ref_in_state;

    
    %% INITIALISE the experiment
    reference=ref_pt7*ones(1,time_dim);

    j=1; % sample counter

    % overall state
    x_traj=zeros(time_dim,size(x0,2));
    x_traj(j,:) = x0;

    % p_t7 trajectory storage
    pt7_traj=zeros(size(reference));
    pt7_traj(j)=x_traj(j,ref_in_state);

    keeprunning=true; % indication whether thje CBC loop should keep running
    
    if(cbc_par.verbose)
        disp('Initialisation Concluded. Reference for the Control Loop evaluated.');
    end

    %% RUN the CBC loop
    while(keeprunning)        
        %% CALCULATE (and display) controller input
        [u, errinfo] = p_control(reference(j), pt7_traj(j), cbc_par.ctrl, Ts, errinfo, (j==1));
        u_traj = [u_traj u];% we need to save it as many times as TC

        %% DISPLAY diagnostics
        if(cbc_par.verbose==true)
            clc;
            disp(['Control iteration number: ',num2str(j)]);
            disp('STATE: EVALUATION OF CONTROL ACTION u(k).');
            disp(['reference value: ',num2str(reference(j))])
            disp(['state: ',num2str(pt7_traj(j))])
            disp(['Control input: ',num2str(u)])
        end


        %% SIMULATE the plant
        sim.het.parameters('t7_toxicity')=u_traj(end);
        sim=sim.push_het();
        [~,xout]=ode15s(@sim.ss_model, [0, Ts], [x_traj(j,:)], sim.opt);
        
        % record
        x_traj(j+1,:)=xout(end,:);
        pt7_traj(j+1)=x_traj(j,ref_in_state);

        % CHECK for steady state
%         dxdt=sim.ss_model(j*Ts,x_traj(j+1,:));
%         if(abs(dxdt)<0.001)
%             keeprunning=false;
%         end
    
        %% UPDATE counter and see if we should keep running the loop
        endtime=j;
        j=j+1;
        if(j==time_dim)
            keeprunning=false;
        end
    end

    %% OUTPUT variables of interest
    STATES.x_traj = x_traj(1:j,:);
    STATES.x_end=x_traj(j,:);
    STATES.pt7_traj=pt7_traj(1:j);
    STATES.pt7_end=pt7_traj(j);

    CONTROL.u_end=u_traj(end);
    CONTROL.u_traj=u_traj;
end

%% CONTROLLER
function [u, errinfo] = p_control(ref, y, ctrl, Ts, errinfo, ref_just_changed)
    % Compute error
    err = ref-y ;

    % Compute integral of error
    errint = errinfo.errint+Ts*(errinfo.preverr+err)/2; % include the latest error in the integral
    % BUT clamp the integral - restore if saturation limit reached
    if(abs(errint)>ctrl.clamp)
        errint=errinfo.errint;
    end

    % Compute derivative of error - only if reference not just changed to avoid derivatuive jumps
    if(ref_just_changed)
        errder = 0;
    else
        errder=(err-errinfo.preverr)/2;
    end
    

    % Compute the PID control signal
    u = ctrl.Kp*err + ctrl.Ki*errint + ctrl.Kd*errder;

    % saturation
    if u <= ctrl.min_u
        u = ctrl.min_u ;
    end
    if u >= ctrl.max_u
        u = ctrl.max_u ;
    end

    % Update error information
    errinfo.preverr=err; % the current error will be 'previous' in the next step
    errinfo.errint=errint; % record error integral
end