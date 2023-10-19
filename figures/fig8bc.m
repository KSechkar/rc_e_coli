%% fig8bc.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 8: b,c

% Showcasing how increasing our controller's amplifier gain xi reduces the 
% adaptation error caused by the expression of extra synthetic mRNA 

%% CLEAR all variables

addpath(genpath('..'))

close all
clear

%% SET UP the simulators for closed loop

% 2 different setups - ideal AIF and realistic scenario
sim=cell_simulator;

sim.init_conditions('s')=0.5;

sim=sim.load_heterologous_and_external('pi_controller','constant_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim.ext.input_func_parameters('inducer_level')=1; % inducer level: just full expression of disturbing gene
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

%% DEFINE parameter ranges to be tested

sim.het.parameters('a_dist')=500;
sim=sim.push_het();
cdists=linspace(0,1500,11);
delta_sigmas=linspace(-0.43,0.43,11);

% record the default nutrient quality
default_sigma=sim.init_conditions('s');

%% RUN simulations - closed loop

% initialise array of steady-state values
ss_values=zeros(size(cdists,2),size(delta_sigmas,2));

for i=1:size(cdists,2)
    for j=1:size(delta_sigmas,2)
        disp(['Testing c_dist=',num2str(cdists(i))...
            ' delta(sigma)=',num2str(delta_sigmas(j))])
        
        % set plasmid concentration
        sim.het.parameters('c_dist')=cdists(i);
        % set amplifier gain
        sim.init_conditions('s')=default_sigma+delta_sigmas(j);
        sim = sim.push_het();
        
        % simulate!
        sim = sim.simulate_model;
        
        % find steady-state D value in given conditions
        ss_values(i,j)=sim.x(end,9+sim.num_het+1); % sens is the first gene in the list
    end
end

%% FIND reference p_sens value for 0 disturbance and 0 change in sigma - closed loop

for i=1:size(cdists,2)
    for j=1:size(delta_sigmas,2)
        if(cdists(i)==0 && delta_sigmas(j)==0)
            ref_value=ss_values(i,j);
        end
    end
end


%% SET UP the simulations for open loop

sim_openloop=cell_simulator;

sim_openloop.init_conditions('s')=sim.init_conditions('s');

sim_openloop=sim_openloop.load_heterologous_and_external('pi_controller','constant_inducer'); % load the het. gene and ext. inp. modules

% disturbance signal parameters
sim_openloop.ext.input_func_parameters('inducer_level')=sim.ext.input_func_parameters('inducer_level'); % inducer level: just full expression of disturbing gene

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

% simulation parameters
sim_openloop.tf =  sim_openloop.tf;
sim_openloop.opt = sim_openloop.opt;

%% RUN simulations and find reference value of p_sens - closed loop

% initialise array of steady-state values
ss_values_openloop=zeros(size(cdists,2),size(delta_sigmas,2));

disp('OPEN LOOP');
for i=1:size(cdists,2)
    for j=1:size(delta_sigmas,2)
        disp(['Testing c_dist=',num2str(cdists(i))...
            ' delta(sigma)=',num2str(delta_sigmas(j))])
        
        % set plasmid concentration
        sim_openloop.het.parameters('c_dist')=cdists(i);
        % set amplifier gain
        sim_openloop.init_conditions('s')=default_sigma+delta_sigmas(j);
        sim_openloop = sim_openloop.push_het();
        
        % simulate!
        sim_openloop = sim_openloop.simulate_model;
        
        % find steady-state D value in given conditions
        ss_values_openloop(i,j)=sim_openloop.x(end,9+sim_openloop.num_het+1); % sens is the first gene in the list
    end
end

% find reference
for i=1:size(cdists,2)
    for j=1:size(delta_sigmas,2)
        if(cdists(i)==0 && delta_sigmas(j)==0)
            ref_value_openloop=ss_values_openloop(i,j);
        end
    end
end

%% FIGURE 8 b - closed loop

Fa = figure('Position',[0 0 239 220]);
set(Fa, 'defaultAxesFontSize', 9)
set(Fa, 'defaultLineLineWidth', 1.25)

% with aif controller
hmap=heatmap(flip(ss_values./ref_value,1),'ColorLimits',[0.4,1.6],'Colormap',parula);

% display labels
x_text_labels=string(delta_sigmas);
for i=1:size(delta_sigmas,2)
    midpoint=round(size(delta_sigmas,2)/2,0);
    if (i~=1 && i~=size(delta_sigmas,2) && i~=midpoint)
        x_text_labels(i)='';
    else
        if(delta_sigmas(i)~=0)
            x_text_labels(i)=[num2str(delta_sigmas(i))];
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
xlabel('\Delta \sigma, nutr. qual. change');

hmap.GridVisible = 'off';

%% FIGURE 8 c - open loop

Fb = figure('Position',[0 0 239 220]);
set(Fb, 'defaultAxesFontSize', 9)
set(Fb, 'defaultLineLineWidth', 1.25)

% with aif controller
hmap=heatmap(flip(ss_values_openloop./ref_value_openloop,1),'ColorLimits',[0.4,1.6],'Colormap',parula);

% display labels
x_text_labels=string(delta_sigmas);
for i=1:size(delta_sigmas,2)
    midpoint=round(size(delta_sigmas,2)/2,0);
    if (i~=1 && i~=size(delta_sigmas,2) && i~=midpoint)
        x_text_labels(i)='';
    else
        if(delta_sigmas(i)~=0)
            x_text_labels(i)=[num2str(delta_sigmas(i))];
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
xlabel('\Delta \sigma, nutr. qual. change');
% title('p_{sens}:p_{sens}^{0, 0.5}, rel. sensor prot. conc.')

hmap.GridVisible = 'off';

%% FUNCTION for getting translation elongation rate
function D=get_D(sim,ss,t)
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

    kcmh=par('kcm').*h;
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    % ribosome dissociation constants
            k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

            % heterologous genes rib. dissoc. constants
            k_het=ones(sim.num_het,1); % initialise with default value 1
            if(sim.num_het>0)
                for i=1:sim.num_het
                    k_het(i)=sim.form.k(e,...
                    sim.parameters(['k+_',sim.het.names{i}]),...
                    sim.parameters(['k-_',sim.het.names{i}]),...
                    sim.parameters(['n_',sim.het.names{i}]),...
                    kcmh);
                end
            end

    D=1+(m_a./k_a+m_r./k_r+sum(ss_het(1:sim.num_het)./(k_het.')))./...
            (1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    ext_inp=sim.ext.input(ss,t);
end