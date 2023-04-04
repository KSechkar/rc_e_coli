%% fig3b.m

% WINNER-TAKES-ALL PHENOMENON
% Figure 3: b

% Two bistable switches (self-activating genes) in the same cell exhibiting
% winner-takes-all behaviour, when the activation of one switch may prevent
% the activation of the other

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP and RUN the simulators to get three trajectories

% switch 2 only, both switches on, switch 1 only
f1s={25 125 1000};


sims={cell_simulator cell_simulator cell_simulator};

for scenario=1:size(f1s,2)  
    sims{scenario}=sims{scenario}.load_heterologous_and_external('two_switches','pulse_inducer'); % load the het. gene and ext. inp. modules
    
    % inducer addition time
    inducer_add_time=10;
    
    % 'inducer levels'
    sims{scenario}.het.parameters('f1')=f1s{scenario}; % inducer for switch 1
    sims{scenario}.het.parameters('f2')=125; % inducer for switch 2
    
    % max transcription rates
    sims{scenario}.het.parameters('a_switch1')=3000;
    sims{scenario}.het.parameters('a_switch2')=3000;
    % Switch 1:
    sims{scenario}.het.parameters('K_dna(switch1)-switch1f1')=5e4; % gene reg. Hill constant
    sims{scenario}.het.parameters('eta_dna(switch1)-switch1f1')=2; % gene reg. Hill coefficient
    sims{scenario}.het.parameters('baseline1')=0.1; % baseline promoter activity
    % Switch 2:
    sims{scenario}.het.parameters('K_dna(switch2)-switch2f2')=5e4; % gene reg. Hill constant
    sims{scenario}.het.parameters('eta_dna(switch2)-switch2f2')=2; % gene reg. Hill coefficient
    sims{scenario}.het.parameters('baseline2')=0.1; % baseline promoter activity
    
    % DO NOT TOUCH!
    sims{scenario}.ext.input_func_parameters('inducer_base_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
    sims{scenario}.ext.input_func_parameters('pulse_value_prop')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
    sims{scenario}.ext.input_func_parameters('pulse_start_time')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
    sims{scenario}.ext.input_func_parameters('pulse_duration')=inducer_add_time;% disturbance flicks transcription reg. func. from 0 to 1 at t=30
       
    % push amended parameter values
    sims{scenario}=sims{scenario}.push_het();
    
    % simulate
    sims{scenario}.tf =  50;
    sims{scenario}.opt = odeset('reltol',1.e-6,'abstol',1.e-9);
    sims{scenario} = sims{scenario}.simulate_model;
end

%% FIND equilibria across different inducer 1 concentrations

f1_sweep=f1s{1}:12.5:f1s{3};

eqbs=zeros(2, size(f1_sweep,2));

% setting up the simulator
sim=cell_simulator; % initialise simulator
sim=sim.load_heterologous_and_external('two_switches','pulse_inducer'); % load the het. gene and ext. inp. modules
sim.het.parameters=sims{1}.het.parameters; % same parameters as before
sim.ext.input_func_parameters=sims{1}.ext.input_func_parameters; % same parameters as before
sim=sim.push_het();

sim.tf =  50;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

% run simulations, record equilibria
for i=1:size(f1_sweep,2)
    sim.het.parameters('f1')=f1_sweep(i);
    sim=sim.push_het();
    sim = sim.simulate_model;
    x_het=sim.x(:,10:(9+2*sim.num_het));
    eqbs(1,i)=x_het(end,3);
    eqbs(2,i)=x_het(end,4);
end

%% PHASE PLANE PLOT

% create figure
Fig_pp = figure('Position',[0 0 385 295]);
set(Fig_pp, 'defaultAxesFontSize', 9)
set(Fig_pp, 'defaultLineLineWidth', 1.25)

hold on

% plot the three trajectories
traj_colours={[0 0.75 0.75],[0.375 0.75 0.375],[0.75 0.75 0]}; % colours of trajectories
for scenario=1:size(f1s,2) 
    % find time point at which inducer was added
    for i=1:size(sims{scenario}.t,1)
        if(sims{scenario}.t(i)>inducer_add_time)
            i_start_obs=i-1;
            break
        end
    end
    
    x_het=sims{scenario}.x(:,10:(9+2*sims{scenario}.num_het));
    
    plot(x_het(i_start_obs:end,3),x_het(i_start_obs:end,4),'-','Color',traj_colours{scenario},'LineWidth',1);
end

% plot equilibria for different inducer 1 concentrations
plot(eqbs(1,:),eqbs(2,:),'Color',[1 0 0],'LineWidth',2)

legend('f_1=25nM','f_1=125uM','f_1=1000uM',...
    'Location','northeast','FontName','Arial')

xlabel('p_{s1}, conc. of switch 1 prot. [nM]','FontName','Arial')
ylabel('p_{s2}, conc. of switch 2 prot. [nM]','FontName','Arial')

xlim([0 5e5])
ylim([0 5e5])
xticks(0:1e5:5e5)
yticks(0:1e5:5e5)

grid on
axis square
box on
hold off