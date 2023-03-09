%% fig5f.m

% PROPORTIONAL-INTEGRAL CONTROLLER
% Figure 5: f

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
sim.ext.input_func_parameters('step_time')=100; % disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.ext.input_func_parameters('slope_duration')=0.1;% disturbance flicks transcription reg. func. from 0 to 1 at t=30
sim.het.parameters('c_dist')=100; % gene copy number
sim.het.parameters('a_dist')=300; % max. gene transcription rate

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

% simulation parameters
sim.tf =  200;
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9);

%% DEFINE plasmid concs. to be tested

sim.het.parameters('a_dist')=1000;
sim=sim.push_het();
cdists=linspace(0,150,3);
chis=linspace(0,7500.*sim.het.parameters('c_amp'),3);

% initialise array of adaptation errors
adaptation_errors=zeros(size(cdists,2),size(chis,2));

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
        
        % find D value before disturbance
        for t_counter=1:size(sim.t,1)
            if(sim.t(t_counter)>sim.ext.input_func_parameters('step_time'))
                t_before=sim.t(t_counter-1);
                x_before=sim.x(t_counter-1,:);
                D_before=get_D(sim,x_before,t_before);
                break
            end
        end

        % find steady-state D value after disturbance
        D_after=get_D(sim,sim.x(end,:),sim.t(end));

        % record adaptation error
        adaptation_errors(i,j)=abs((D_after-D_before)./D_before)*100;

    end
end


%% FIGURE 5 f

Fg = figure('Position',[0 0 239 220]);
set(Fg, 'defaultAxesFontSize', 9)
set(Fg, 'defaultLineLineWidth', 1.25)

% with aif controller
hmap=heatmap(flip(adaptation_errors,1));

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