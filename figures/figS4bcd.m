%% figS4ab.m

% WINNER-TAKES-ALL PHENOMENON
% Figure S4bcd

% Two bistable switches (self-activating genes) in the same cell exhibiting
% winner-takes-all behaviour, when the activation of one switch may prevent
% the activation of the other. Here, we show that the WTA effect can be
% reproduced by changing the sequence in which we activate the bistable
% switches

%% CLEAR all variables

addpath(genpath('..'))

close all
clear
%% DEFINE the addition times and the EQUAL inducer concentrations added

added_first={[1 0],... % switch1 activated first
    [1 1],... % both switches activated simultaneously
    [0 1]... % switch2 activated first
    };
add_time_gap=15; % time gap between addition of different inducers
fi=20;

%% SET UP the simulators with default circuit parameters in WTA simulations

sim={}; % initialise the cell simulator array

for scenario=1:size(added_first,2)
    sim{scenario}=cell_simulator;
    sim{scenario}=sim{scenario}.load_heterologous_and_external('two_switches','constant_inducer'); % load the het. gene and ext. inp. modules
    
    % 'inducer levels' - originally zero
    sim{scenario}.het.parameters('f1')=0; % inducer for switch 1
    sim{scenario}.het.parameters('f2')=0; % inducer for switch 2
    
    % max transcription rates
    sim{scenario}.het.parameters('a_switch1')=3000;
    sim{scenario}.het.parameters('a_switch2')=3000;
    % Switch 1:
    sim{scenario}.het.parameters('K_dna(switch1)-switch1f1')=5e3; % gene reg. Hill constant
    sim{scenario}.het.parameters('eta_dna(switch1)-switch1f1')=2; % gene reg. Hill coefficient
    sim{scenario}.het.parameters('baseline1')=0.05; % baseline promoter activity
    % Switch 2:
    sim{scenario}.het.parameters('K_dna(switch2)-switch2f2')=5e3; % gene reg. Hill constant
    sim{scenario}.het.parameters('eta_dna(switch2)-switch2f2')=2; % gene reg. Hill coefficient
    sim{scenario}.het.parameters('baseline2')=0.05; % baseline promoter activity
    
    % inducer protein binding
    sim{scenario}.het.parameters('K_switch1-f1')=1000; % dissociation constant for inducer-protein binding
    sim{scenario}.het.parameters('K_switch2-f2')=1000; % dissociation constant for inducer-protein binding
    
    % DO NOT TOUCH!
    sim{scenario}.ext.input_func_parameters('inducer_level')=1;
       
    % push amended parameter values
    sim{scenario}=sim{scenario}.push_het();
    
    % simulation parameters
    sim{scenario}.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % ODE simulation tolerances
end

%% GET TRAJECTORIES of a SINGLE bistable switch with different inducer levels

% initialise trajectory storage
t_trajs={}; % time axes of simulation
pswitch1_trajs={}; % p_switch1 levels
pswitch2_trajs={}; % p_switch2 levels

for scenario=1:size(added_first,2)
    %% SIMULATE the system prior to induction
    sim{scenario}.tf =  72; % simulation time frame
    sim{scenario} = sim{scenario}.simulate_model;
    % record trajectories
    t_trajs{scenario}=sim{scenario}.t;
    x_het=sim{scenario}.x(:,10:(9+2*sim{scenario}.num_het));
    pswitch1_trajs{scenario}=x_het(:,3);
    pswitch2_trajs{scenario}=x_het(:,4);

    %% SET the initial condition according to the uninduced state
    state_before = sim{scenario}.x(end,:);
     % mRNAs
    sim{scenario}.init_conditions('m_a')=state_before(1);
    sim{scenario}.init_conditions('m_r')=state_before(2);
    % proteins
    sim{scenario}.init_conditions('p_a')=state_before(3);
    sim{scenario}.init_conditions('R')=state_before(4);
    % tRNAs
    sim{scenario}.init_conditions('tc')=state_before(5);
    sim{scenario}.init_conditions('tu')=state_before(6);
    % inactivated ribosomes
    sim{scenario}.init_conditions('Bcm')=state_before(7);
    % nutrient quality and chloramphenicol
    sim{scenario}.init_conditions('s')=state_before(8);
    sim{scenario}.init_conditions('h')=state_before(9);
    % heterologous
    x_het=state_before(10 : (9+2*sim{scenario}.num_het+sim{scenario}.num_misc) );
    for j=1:sim{scenario}.num_het
        % mRNA
        sim{scenario}.het.init_conditions(['m_',sim{scenario}.het.names{j}])=x_het(j);
        % protein
        sim{scenario}.het.init_conditions(['p_',sim{scenario}.het.names{j}])=x_het(sim{scenario}.num_het+j);
    end
    % miscellaneous
    for j=1:sim{scenario}.num_misc
        sim{scenario}.het.init_conditions([sim{scenario}.het.misc_names{j}])=x_het(2*sim{scenario}.num_het+j);
    end
    sim{scenario}=sim{scenario}.push_het;
    
    %% INDUCE first switch to be activated and SIMULATE the system
    % add inducer
    sim{scenario}.het.parameters('f1')=fi*added_first{scenario}(1);
    sim{scenario}.het.parameters('f2')=fi*added_first{scenario}(2);
    sim{scenario}=sim{scenario}.push_het;

    % simulate
    sim{scenario}.tf=add_time_gap;
    sim{scenario}=sim{scenario}.simulate_model;

    % record trajectories
    t_trajs{scenario}=[t_trajs{scenario};t_trajs{scenario}(end)+sim{scenario}.t];
    x_het=sim{scenario}.x(:,10:(9+2*sim{scenario}.num_het));
    pswitch1_trajs{scenario}=[pswitch1_trajs{scenario};x_het(:,3)];
    pswitch2_trajs{scenario}=[pswitch2_trajs{scenario};x_het(:,4)];

    %% SET the initial condition ccording to the system's state after first induction
    state_before = sim{scenario}.x(end,:);
    % mRNAs
    sim{scenario}.init_conditions('m_a')=state_before(1);
    sim{scenario}.init_conditions('m_r')=state_before(2);
    % proteins
    sim{scenario}.init_conditions('p_a')=state_before(3);
    sim{scenario}.init_conditions('R')=state_before(4);
    % tRNAs
    sim{scenario}.init_conditions('tc')=state_before(5);
    sim{scenario}.init_conditions('tu')=state_before(6);
    % inactivated ribosomes
    sim{scenario}.init_conditions('Bcm')=state_before(7);
    % nutrient quality and chloramphenicol
    sim{scenario}.init_conditions('s')=state_before(8);
    sim{scenario}.init_conditions('h')=state_before(9);
    % heterologous
    x_het=state_before(10 : (9+2*sim{scenario}.num_het+sim{scenario}.num_misc) );
    for j=1:sim{scenario}.num_het
        % mRNA
        sim{scenario}.het.init_conditions(['m_',sim{scenario}.het.names{j}])=x_het(j);
        % protein
        sim{scenario}.het.init_conditions(['p_',sim{scenario}.het.names{j}])=x_het(sim{scenario}.num_het+j);
    end
    % miscellaneous
    for j=1:sim{scenario}.num_misc
        sim{scenario}.het.init_conditions([sim{scenario}.het.misc_names{j}])=x_het(2*sim{scenario}.num_het+j);
    end
    sim{scenario}=sim{scenario}.push_het;

    %% INDUCE the rest of switches and SIMULATE the system
    % add inducer
    sim{scenario}.het.parameters('f1')=fi;
    sim{scenario}.het.parameters('f2')=fi;
    sim{scenario}=sim{scenario}.push_het;

    % simulate
    sim{scenario}.tf=72;
    sim{scenario}=sim{scenario}.simulate_model;

    % record trajectories
    t_trajs{scenario}=[t_trajs{scenario};t_trajs{scenario}(end)+sim{scenario}.t];
    x_het=sim{scenario}.x(:,10:(9+2*sim{scenario}.num_het));
    pswitch1_trajs{scenario}=[pswitch1_trajs{scenario};x_het(:,3)];
    pswitch2_trajs{scenario}=[pswitch2_trajs{scenario};x_het(:,4)];

end

%% PLOT TRAJECTORIES

% create figure
Figs={};
for scenario=1:size(added_first,2)
    Fig{scenario} = figure('Position',[0 0 250 186]);
    set(Fig{scenario}, 'defaultAxesFontSize', 9)
    set(Fig{scenario}, 'defaultLineLineWidth', 1.25)
    
    hold on
    
    switch1_colour=[0.75 0.75 0];
    switch2_colour=[0 0.75 0.75];
    
    % plot the three trajectories
    plot(t_trajs{scenario}-67, ...
        pswitch1_trajs{scenario},'-', ...
        'Color',switch1_colour, 'DisplayName','p_{s1}');
    plot(t_trajs{scenario}-67, ...
        pswitch2_trajs{scenario},'--', ...
        'Color',switch2_colour, 'DisplayName','p_{s2}');
    
    
    legend('Location','northwest','FontName','Arial')
    
    xlabel('time [h]','FontName','Arial')
    ylabel('p_{si}, conc. of switch prot. [nM]','FontName','Arial')
    
    xticks(0:5:30)
    xlim([0 30])
    ylim([0 5e5])
    
    grid on
    axis square
    box on
    hold off
end

