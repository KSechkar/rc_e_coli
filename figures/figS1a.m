%% figS1a.m

% MCMC FITTING
% Figure S1: a

% Sample parameter combinations from the MCMC chains and display the corresponding model
% predictions

% Warning: requires the DREAM package

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% DEFINE number of samples

num_samples=100;

%% LOAD the MCMC chains
load DREAM_fitting_outcome.mat
ParSet = genparset(chain); DREAMPar.N = size(chain,3);
Pars = ParSet ( floor ( 0.75 * size(ParSet,1) ) : size(ParSet,1), 1 : DREAMPar.d ); % take the last 25% of the posterior samples
N_Pars = size(Pars,1); % get number of posterior samples

%% DEFINE parameter names

params = {'a_r/a_a',... % metabolic gene transcription rate
    'nu_max',...        % max. tRNA aminoacylatio rate
    'K_t',...           % MM constants for translation elongation and tRNA charging rates
    'kcm'               % chloramphenicol binding rate constant
    };

%% LOAD Experimental data (to compare with the fit)

% read the experimental dataset (eq2 strain of Scott 2010)
dataset = readmatrix('data/growth_rib_fit_notext.csv');

% nutrient qualities are equally log-spaced points
nutr_quals=logspace(log10(0.08),log10(0.5),6);

% get inputs: nutrient quality and h; get outputs: l and rib mass frac
data.xdata=[]; % initialise inputs array
data.ydata=[]; % intialise outputs array
for i = 1:size(dataset,1)
    if(dataset(i,1)>0.3)
        % inputs
        nutr_qual = nutr_quals(fix((i-1)/5)+1); % records start from worst nutrient quality
        h = dataset(i,4)*1000; % all h values for same nutr quality same go one after another. Convert to nM from uM!
        data.xdata=[data.xdata; [nutr_qual,h]];
    
        % outputs
        l = dataset(i,1); % growth rate (1/h)
        phi_r = dataset(i,3); % ribosome mass fraction
        data.ydata=[data.ydata; [l,phi_r]];
    end
end

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 4; % maximum no. iterations (checking if SS reached over first 48 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed

model.ssfun = @(theta,data) rc_ecoli_sos(theta,data,sim,Delta,Max_iter);

%% GET model predictions with fitted parameters
% Find Maximum A Posteriori estimate (mode of the distribution) - these are the values we ultimately use in our model
[~,max_index] = max(ParSet(:,end)); max_index = max_index(1);
fitted_theta = ParSet(max_index,1:DREAMPar.d);

disp(['Fitted parameters: Theta=',num2str(fitted_theta)]) % print resultant sum of squared errors

ymodel_fitted=dream_modelfun(fitted_theta,data.xdata,sim,Delta,Max_iter); % find model predictions
disp(['Fitted parameters: SOS=',num2str(sum((ymodel_fitted-data.ydata).^2))]) % print resultant sum of squared errors

%% GET model predictions with sampled parameters - randomly draw from the chains

ymodels={};

tried_thetas=zeros(num_samples,size(params,2));
for sample_cntr=1:num_samples
    % randomly draw parameter values using CDFs
    random_draw=randi(N_Pars);
    tried_thetas(sample_cntr,:)=Pars(random_draw,:);

    disp(['Sample ',num2str(sample_cntr),': Theta=',num2str(tried_thetas(sample_cntr,:))]) % print resultant sum of squared errors

    ymodels{sample_cntr}=dream_modelfun(tried_thetas(sample_cntr,:),data.xdata,sim,Delta,Max_iter); % find model predictions
    disp(['Sample ',num2str(sample_cntr),': SOS=',num2str(sum((ymodels{sample_cntr}-data.ydata).^2))]) % print resultant sum of squared errors
end

%% PREPARE THE PLOT

colours={[0.6350 0.0780 0.1840 0.25],...
    [0.4660 0.6740 0.1880 0.25],...
    [0.4940 0.1840 0.5560 0.25],...
    [0.9290 0.6940 0.1250 0.25],...
    [0.8500 0.3250 0.0980 0.25],...
    [0 0.4470 0.7410 0.25]};

colours_dark={}; % colours for different media - fitted and original data (darker than sampled)
for colourind=1:size(colours,2)
    colours_dark{colourind}=colours{colourind}.*[0.66 0.66 0.66 0]+[0 0 0 0.75];
end

Fdata = figure('Position',[0 0 385 265]);
set(Fdata, 'defaultAxesFontSize', 9)
hold on


%% PLOT SAMPLED PREDICTIONS
for sample_cntr=1:num_samples
    % group data by nutrient qualities
    x_grouped={[]};
    y_grouped={[]};
    nutrind=1;
    chlorind=1;
    last_nutr_qual=data.xdata(1,1);
    for i=1:size(data.xdata,1)
        if(data.xdata(i,1)~=last_nutr_qual)
            nutrind=nutrind+1;
            chlorind=1;
            last_nutr_qual=data.xdata(i,1);
            x_grouped{nutrind}=[];
            y_grouped{nutrind}=[];
        end
        x_grouped{nutrind}(chlorind)=ymodels{sample_cntr}(i,1);
        y_grouped{nutrind}(chlorind)=ymodels{sample_cntr}(i,2);
        chlorind=chlorind+1;
    end

    % plot
    for nutrind=1:size(x_grouped,2)
        plot(x_grouped{nutrind},y_grouped{nutrind},'Color',colours{nutrind},'LineWidth',0.5)
    end
end

%% PLOT ORIGINAL DATAPOINTS
colourind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        colourind=colourind+1;
        last_nutr_qual=data.xdata(i,1);
    end
    plot(data.ydata(i,1),data.ydata(i,2),'o','Color',colours_dark{colourind},'LineWidth',2) % real data
end

%% PLOT FITTED PREDICTIONS
% group data by nutrient qualities
    x_grouped={[]};
    y_grouped={[]};
    nutrind=1;
    chlorind=1;
    last_nutr_qual=data.xdata(1,1);
    for i=1:size(data.xdata,1)
        if(data.xdata(i,1)~=last_nutr_qual)
            nutrind=nutrind+1;
            chlorind=1;
            last_nutr_qual=data.xdata(i,1);
            x_grouped{nutrind}=[];
            y_grouped{nutrind}=[];
        end
        x_grouped{nutrind}(chlorind)=ymodel_fitted(i,1);
        y_grouped{nutrind}(chlorind)=ymodel_fitted(i,2);
        chlorind=chlorind+1;
    end

    % plot
    for nutrind=1:size(x_grouped,2)
        plot(x_grouped{nutrind},y_grouped{nutrind},'-x','Color',colours_dark{nutrind},'LineWidth',2)
    end

%% ADD PLOT LABELS
ylabel('\phi_r, ribosomal mass fraction','FontName','Arial');
xlabel('\lambda, growth rate [1/h]','FontName','Arial')
xlim([0.2 1.8])
ylim([0.05 0.35])
axis square
grid on
box on
hold off