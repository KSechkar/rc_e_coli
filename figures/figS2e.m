%% figS1b.m

% MCMC FITTING
% Figure S1: e

% Analyse the model's parameter sensitivity by going through different
% pairwise combinations of the values of the four fitted parameters

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all


%% LOAD experimental data

% read the experimental dataset (eq2 strain of Scott 2010)
dataset = readmatrix('data/growth_rib_fit_notext.csv');

% nutrient qualities are equally log-spaced points
nutr_quals=logspace(log10(0.08),log10(0.5),6);

% get inputs: nutrient quality and h; get outputs: l and rib mass frac
data.xdata=[]; % initialise inputs array
data.ydata=[]; % intialise outputs array
for i = 1:size(dataset,1)
    if dataset(i,1)>0.3
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

% standard measurement errors
l_stdev=0.04467;
phir_stdev=0.018976;

%% DEFINE parameter values to check
% names of parameters to consider
parnames={'\alpha_r : \alpha_a', ...
    '\nu_{max}', ...
    'K_\epsilon = K_\nu', ...
    'k_{cm}'};

% define default parameter values
default_par=cell_params();
default_theta=[default_par('a_r')/default_par('a_a'), ...
    default_par('nu_max'), ...
    default_par('K_e'), ...
    default_par('kcm')];

% ranges of parameter values to consider
points_in_range=75;
ratio_range=logspace(log10(1/10),log10(10),points_in_range)*default_theta(1);
numax_range=logspace(log10(1/10),log10(10),points_in_range)*default_theta(2);
K_range=logspace(log10(1/10),log10(10),points_in_range)*default_theta(3);
kcm_range=logspace(log10(1/10),log10(10),points_in_range)*default_theta(4);
ranges = [ratio_range; numax_range; K_range; kcm_range];

% define parameter combinations to consider
parcombs=nchoosek([1,2,3,4],2);

% initialise array of log-likelihoods
loglikes=zeros(points_in_range,points_in_range);
for i=2:size(parcombs,1)
    loglikes(:,:,i)=zeros(points_in_range,points_in_range);
end

% backup save
bckup.parnames=parnames;
bckup.ranges=ranges;
bckup.parcombs=parcombs;
bckup.loglikes=loglikes;
save('S1e.mat','bckup')

%% RUN the simulations and calculate log-likelihoods
for i_comb=1:size(parcombs,1)
    % set up the simulator
    sim=cell_simulator; % initialise simulator

    % settings
    sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed
    
    % parameters for getting steady state
    sim.tf = 12; % single integraton step timeframe
    Delta = 0.001; % threshold that determines if we're in steady state
    Max_iter = 4; % maximum no. iterations (checking if SS reached over first 48 h)

    % which parameter combination is being considered - x and y histogram ranges
    parcomb=parcombs(i_comb,:);
    x_par=parcomb(1);
    y_par=parcomb(2);

    for i_x=1:points_in_range
        parfor i_y=1:points_in_range
            % set the parameter vactor to considered values
            theta=default_theta;
            theta(x_par)=ranges(x_par,i_x);
            theta(y_par)=ranges(y_par,i_y);

            % get model predictions
            ymodel=dream_modelfun(theta,data.xdata,sim,Delta,Max_iter,sim.parameters('a_a'));

            % get sum of squared errors
            sos = sum((ymodel - data.ydata).^2);
        
            % for testing - print the results
            %disp(['SOS=',num2str(sos)])
        
            % get logarithm of mutually independent normal distributions
            loglikes(i_x,i_y,i_comb)=-log(l_stdev.*sqrt(2.*pi))-0.5.*sos(1)./(l_stdev.^2)... % growth rate
                -log(phir_stdev.*sqrt(2.*pi))-0.5.*sos(2)./(phir_stdev.^2); % ribosomal mass fraction
        end
        disp([parnames{x_par},' vs ' ,parnames{y_par},' : ',num2str(i_x/points_in_range*100),' % completed'])
    end
end
% save results
bckup.loglikes=loglikes;
save('S1e.mat','bckup')

%% LOAD the saved simulation results (optional)

% load('S2e.mat');
% parnames=bckup.parnames;
% ranges=bckup.ranges;
% parcombs=bckup.parcombs;
% loglikes=bckup.loglikes;
% points_in_range=size(ranges,2);


%% PLOT
% heatmap colour axis range
hmap_crange=[-3000,0];

Fs1e = figure('Position',[0 0 900 700]);
set(Fs1e, 'defaultAxesFontSize', 9)
set(Fs1e, 'defaultLineLineWidth', 1.25)

% number of plots in one row of subplot layout
plots_in_row=size(parnames,2)-1;


for i_comb=1:size(parcombs,1)
    % which parameter combination is being considered - x and y histogram ranges
    parcomb=parcombs(i_comb,:);
    x_par=parcomb(1);
    y_par=parcomb(2);

    % select subplot
    subplot(plots_in_row,plots_in_row,10-(plots_in_row*(y_par-2)+(plots_in_row+1-x_par)))

    % plot heatmap
    hmap=heatmap(flip(loglikes(:,:,i_comb).',1),'ColorMap', parula);
%     hmap=heatmap(flip(loglikes{i_comb},1),'ColorMap', jet(100));

    % display labels
    x_text_labels=string(ranges(x_par,:));
    for i=1:points_in_range
        midpoint=round(points_in_range/2,0);
        if (i~=1 && i~=points_in_range && i~=midpoint)
            x_text_labels(i)='';
        else
            if(ranges(x_par,i)~=0)
                % find which power of 10 to use in the engineering notation
                for power_of_ten=-100:100
                    if(10^power_of_ten>ranges(x_par,i))
                        break;
                    end
                end
                power_of_ten=power_of_ten-1;
    
                % make axis label
                x_text_labels(i)=[num2str(round(ranges(x_par,i)/10^power_of_ten,2)),'e',num2str(power_of_ten)];
            else
                x_text_labels(i)='0';
            end
        end
    end
    
    ys_flip=flip(ranges(y_par,:));
    y_text_labels=string(ys_flip);
    for i=1:size(ys_flip,2)
        midpoint=round(size(ys_flip,2)/2,0);
        if (i~=1 && i~=size(ys_flip,2) && i~=midpoint)
            y_text_labels(i)=' ';
        else
            if(ys_flip(i)~=0)
                % find which power of 10 to use in the engineering notation
                for power_of_ten=-100:100
                    if(10^power_of_ten>ys_flip(i))
                        break;
                    end
                end
                power_of_ten=power_of_ten-1;
    
                % make axis label
                y_text_labels(i)=[num2str(round(ys_flip(i)/10^power_of_ten,2)),'e',num2str(power_of_ten)];
            else
                y_text_labels(i)='0';
            end
        end
    end
    hmap.XDisplayLabels = x_text_labels;
    hmap.YDisplayLabels = y_text_labels;
    
    xlabel(parnames{x_par});
    ylabel(parnames{y_par});
    
    caxis(hmap,hmap_crange);
    %set(gca,'ColorScale','log')
    hmap.GridVisible = 'off';

end
