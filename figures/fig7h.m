%% fig7h.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 7: h

% Showcasing how increasing our controller's amplifier gain xi reduces the 
% adaptation error caused by the expression of extra synthetic mRNA 

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP the simulators for both cases

% 2 different setups - ideal AIF and realistic scenario
sim=cell_simulator;

sim=sim.load_heterologous_and_external('pi_controller','step_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_base_level')=0; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('inducer_final_level')=1; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('step_time')=30; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('slope_duration')=0.1;% disturbance flicks transcription reg. func. from 0 to 1 at t=30

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

sim.het.parameters('a_amp')=100; % integral controller amplifier gene copy number
sim.het.parameters('a_amp')=4000; % integral controller amplifier transcription rate
   
% push amended parameter values
sim=sim.push_het();

% simulation parameters
sim.tf =  72;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% DEFINE plasmid concs. to be tested

sim.het.parameters('a_dist')=500;
sim=sim.push_het();
cdists=linspace(0,200,11);
chis=linspace(0,8e5,11);

% initialise array of adaptation errors
rel_psens=zeros(size(cdists,2),size(chis,2));

%% RUN simulations - with controller

% as we scale the integral controller gain, ratio of (c_act*a_act):(c_anti*a_anti) stays the same
u=(sim.het.parameters('c_act').*sim.het.parameters('a_act'))./...
    (sim.het.parameters('c_anti').*sim.het.parameters('a_anti'));

thetachi=sim.het.parameters('kb_anti')./(sim.het.parameters('c_anti').*sim.het.parameters('a_anti'));

% initialise the array in which the results are stored
phi_dists_with = zeros(1,size(cdists,2));

for i=1:size(cdists,2)
    for j=1:size(chis,2)
        disp(['Testing c_dist=',num2str(cdists(i))...
            ' chi=',num2str(chis(j))])
        
        % set plasmid concentration
        sim.het.parameters('c_dist')=cdists(i);
        % set amplifier gain
        sim.het.parameters('a_amp')=chis(j)./sim.het.parameters('c_amp');
        sim = sim.push_het();
        
        % simulate!
        sim = sim.simulate_model;
        
        % find p_sens value before disturbance
        for t_counter=1:size(sim.t,1)
            if(sim.t(t_counter)>sim.ext.input_func_parameters('step_time'))
                psens_before=sim.x(t_counter-1,9+sim.num_het+1);
                break
            end
        end

        % find steady-state p_sens value after disturbance
        psens_after=sim.x(end,9+sim.num_het+1);

        % record adaptation error
        rel_psens(i,j)=psens_after./psens_before;

    end
end


%% FIGURE 7 h

Fg = figure('Position',[0 0 239 220]);
set(Fg, 'defaultAxesFontSize', 9)
set(Fg, 'defaultLineLineWidth', 1.25)

% with aif controller
hmap=heatmap(flip(rel_psens,1),'ColorLimits',[0.85,1.15],'Colormap',parula);

% display labels
x_text_labels=string(chis);
for i=1:size(chis,2)
    midpoint=round(size(chis,2)/2,0);
    if (i~=1 && i~=size(chis,2) && i~=midpoint)
        x_text_labels(i)='';
    else
        if(chis(i)~=0)
            x_text_labels(i)=[num2str(chis(i)/1e5),'e5'];
        else
            x_text_labels(i)='0';
        end
    end
end

cdists_flip=flip(cdists);
y_text_labels=string(cdists_flip);
for i=1:size(cdists_flip,2)
    midpoint=round(size(cdists_flip,2)/2,0);
    if (i~=1 && i~=size(cdists_flip,2) && i~=midpoint)
        y_text_labels(i)=' ';
    else
        if(cdists_flip(i)~=0)
            y_text_labels(i)=num2str(cdists_flip(i));
        else
            y_text_labels(i)='0';
        end
    end
end
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

ylabel('c_{dist}, dist. gene conc. [nM]');
xlabel('\chi, amplifier gain');

hmap.GridVisible = 'off';