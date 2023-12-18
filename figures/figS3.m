% figS3.m

% HETEROLOGOUS GENE EXPRESSION AND ITS EFFECT ON CELLULAR VARIABLES
% Figure S3: all subfigures

% Show how heterologous gene expression affects cellular variables in
% culture media of different qualities

%% CLEAR all variables

addpath(genpath('..'))

clear
close all

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

sim.init_conditions('s')=0.25;

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 40; % maximum no. iterations (checking if SS reached over first 48 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

par=sim.parameters; % the model's parameters: saved separately for ease of notation

% save simulation options to pass on to parallel computing instances
tf=sim.tf;
opt=sim.opt;
%% DEFINE culture media

sigmas=linspace(0.05,1,25); %logspace(log(0.05)/log(10),0,25); % range of culture media's nutrient qualities

%% DEFINE gene expression paramaters

a_xtra=1000; % transcription rate of the heterologous gene PER PLASMID

sim=sim.load_heterologous_and_external('one_constit','no_ext'); % load heterologous gene expression module
sim.het.parameters('a_xtra')=a_xtra;
sim=sim.push_het;

%plasmid_concs = [0,linspace(0,log(par('a_a')*2/a_xtra)/log(10),24)]; % range of plasmid concentrations
plasmid_concs = linspace(0,par('a_a')*2/a_xtra,25);

%% INITIALISE arrays where obtained values of phys. variables will be stored

% values without het. prot. exp.
e0=zeros(size(sigmas,2),1);  % elongation rates
F_r0=zeros(size(sigmas,2),1); % ribosomal gene transcription regulation function
T0=zeros(size(sigmas,2),1); % proportional to 1/ppGpp concentration
nu0=zeros(size(sigmas,2),1); % tRNA charging rate

% numerically obtained values
es = zeros(size(sigmas,2),size(plasmid_concs,2)); % elongation rates
F_rs = zeros(size(sigmas,2),size(plasmid_concs,2)); % ribosomal gene transcription regulation function
Ts = zeros(size(sigmas,2),size(plasmid_concs,2)); % proportional to 1/ppGpp concentration
nus = zeros(size(sigmas,2),size(plasmid_concs,2)); % tRNA charging rate

%% RUN with different burden levels

% get actual model predictions
for i=1:size(sigmas,2)
    for j=1:size(plasmid_concs,2)
        parsim=cell_simulator(); % create cell simulator for parallel computing
        parsim=parsim.load_heterologous_and_external('one_constit','no_ext'); % load heterologous gene expression module
        parsim.het.parameters('a_xtra')=a_xtra;
        parsim.push_het();
        parsim.opt=opt;
        parsim.tf=tf;

        parsim.init_conditions('s')=sigmas(i);
        for het_cntr=1:parsim.num_het
            parsim.het.parameters(['c_',parsim.het.names{het_cntr}])=plasmid_concs(j);
        end
        parsim = parsim.push_het(); % push initial condition and parameters to main object
    
        ss=get_steady(parsim,Delta,Max_iter);
        
        [e,nu,F_r,T]=get_enuFrT(parsim,ss); % get values of variables
        
        % record
        es(i,j)=e;%(e-e0(i))/e0(i);
        F_rs(i,j)=F_r;%(F_r-F_r0(i))/F_r0(i);
        nus(i,j)=nu;
        Ts(i,j)=T;
    end
    disp([num2str(i),' out of ',num2str(size(sigmas,2)),' nutirent qualities considered: ',num2str(sigmas(i))])
end

%% FIND RELATIVE CHANGES

e0=es(:,1);
F_r0=F_rs(:,1);
nu0=nus(:,1);
T0=Ts(:,1);
for i=1:size(sigmas,2)
    for j=1:size(plasmid_concs,2)
        es(i,j)=(es(i,j)-e0(i))/e0(i);
        F_rs(i,j)=(F_rs(i,j)-F_r0(i))/F_r0(i);
        nus(i,j)=(nus(i,j)-nu0(i))/nu0(i);
        Ts(i,j)=(Ts(i,j)-T0(i))/T0(i);
    end
end

%% MAKE axis labels
x_text_labels=string(sigmas);
for i=1:size(sigmas,2)
    midpoint=round(size(sigmas,2)/2,0);
    if (i~=1 && i~=size(sigmas,2) && i~=midpoint)
        x_text_labels(i)='';
    else
        if(sigmas(i)~=0)
            x_text_labels(i)=num2str(round(sigmas(i),2));
        else
            x_text_labels(i)='0';
        end
    end
end

plasmid_concs_flip=flip(plasmid_concs);
y_text_labels=string(plasmid_concs_flip);
for i=1:size(plasmid_concs,2)
    midpoint=round(size(plasmid_concs_flip,2)/2,0);
    if (i~=1 && i~=size(plasmid_concs_flip,2) && i~=midpoint)
        y_text_labels(i)=' ';
    else
        if(plasmid_concs_flip(i)~=0)
            y_text_labels(i)=num2str(round(plasmid_concs_flip(i),2));
        else
            y_text_labels(i)='0';
        end
    end
end

%% CREATE FIGURE

F = figure('Position',[0 0 720 560]);
set(F, 'defaultAxesFontSize', 9)
set(F, 'defaultLineLineWidth', 1.25)

%% SUBFIGURE A - TRANSLATION ELONGATION RATES
subplot(2,2,1)
hmap=heatmap(100*abs(flip(es.',1)),'ColorMap', parula);

% axis labels and ticks
xlabel('\sigma, nutr. qual.');
ylabel('c_x, gene DNA conc.');
title('Change in \epsilon due to burden, %');
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

%appaearance of the heatmap
caxis(hmap,[0 1.5e-4]); % colour map
hmap.GridVisible = 'off'; % don't show grid lines

%% SUBFIGURE B - RIBOSOMAL GENE TRANSCRIPTION REGULATION

subplot(2,2,2)
hmap=heatmap(100*abs(flip(F_rs.',1)),'ColorMap', parula);

% axis labels and ticks
xlabel('\sigma, nutr. qual.');
ylabel('c_x, gene DNA conc.');
title('Change in F_r due to burden, %');
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

%appaearance of the heatmap
caxis(hmap,[0 1.5e-4]); % colour map
hmap.GridVisible = 'off'; % don't show grid lines

%% SUBFIGURE C - tRNA AMINOACYLATION RATES
subplot(2,2,3)
hmap=heatmap(100*abs(flip(nus.',1)),'ColorMap', parula);

% axis labels and ticks
xlabel('\sigma, nutr. qual.');
ylabel('c_x, gene DNA conc.');
title('Change in \nu due to burden, %');
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

%appaearance of the heatmap
caxis(hmap,[0 1.5e-4]); % colour map
hmap.GridVisible = 'off'; % don't show grid lines

%% SUBFIGURE D - ppGpp levels
subplot(2,2,4)
hmap=heatmap(100*abs(flip(Ts.',1)),'ColorMap', parula);

% axis labels and ticks
xlabel('\sigma, nutr. qual.');
ylabel('c_x, gene DNA conc.');
title('Change in T due to burden, %');
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

%appaearance of the heatmap
caxis(hmap,[0 1.5e-4]); % colour map
hmap.GridVisible = 'off'; % don't show grid lines

%% FUNCTION for getting growth rate, translation elongation rate and het. prot. mass fraction from the system's steady state
function [e,nu,F_r,T]=get_enuFrT(sim,ss)
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
    k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),par('kcm').*h); % ribosome-mRNA dissociation constant (metab. genes)
    k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),par('kcm').*h); % ribosome-mRNA dissociation constant (rib. genes)
    D=1+(m_a./k_a+m_r./k_r)./(1-par('phi_q')); % denominator in ribosome competition calculations
    B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

    e=sim.form.e(par,tc); % translation elongation rate
    nu=sim.form.nu(par,tu,s); % tRNA aminoacylation rate
    T=tc./tu; % inversely proportional to ppGpp conc.
    F_r=sim.form.F_r(par,T); % ribosomal gene transcription regulation function
end