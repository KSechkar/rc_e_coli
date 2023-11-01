%% fig3a.m
% GROWTH PHENOMENA PREDICTION BY THE MODEL
% Figure 3: a

% Parameter fitting results: model predictions for growth rates and ribosomal
% mass fractions in different media and chloramphenicol concs, compared to
% real-life measurement to which the parameters were fitted. Includes
% linear fits illustrating the bacterial growth laws

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% VECTOR of fitted parameter values

%theta=[0.9998, 4.0140e+03, 1.1233e+03, 3.5527e-04];
theta=[1.03184, 4046.89, 1239.69, 0.000356139];

%% DEFINE starting parameter values (to compare with the fit)

params = {
    {'a_r/a_a', 1,  0} % metabolic gene transcription rate
    {'nu_max', 6000,  0} % max. tRNA aminoacylatio rate
    {'K_t', 80000, 0} % MM constants for translation elongation and tRNA charging rates
    {'kcm', 0.3594/1000, 0} % chloramphenicol binding rate constant
    };

% record original params into a single vector
theta_origin=zeros([size(params,1) 1]);
for i=1:size(theta_origin,1)
    theta_origin(i) = params{i}{2};
end

%% LOAD Experimental data (to compare with the fit)

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

%% SET UP the simulator

sim=cell_simulator; % initialise simulator

% parameters for getting steady state
sim.tf = 12; % single integraton step timeframe
Delta = 0.001; % threshold that determines if we're in steady state
Max_iter = 4; % maximum no. iterations (checking if SS reached over first 750 h)

sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed


%% GET model predictions with fitted parameters
ymodel=dream_modelfun(theta,data.xdata,sim,Delta,Max_iter,sim.parameters('a_a'));


l_stdev=0.04467;
phir_stdev=0.018976;
sos=sum((ymodel-data.ydata).^2);
loglike=-log(l_stdev.*sqrt(2.*pi))-0.5.*sos(1)./(l_stdev.^2)... % growth rate
        -log(phir_stdev.*sqrt(2.*pi))-0.5.*sos(2)./(phir_stdev.^2); % ribosomal mass fraction
disp(['SOS=',num2str(sos)]) % print resultant sum of squared errors
disp(['loglike=',num2str(loglike)])
% get points for original/strarting parameter values (optional)
ymodel_origin=dream_modelfun(theta_origin,data.xdata,sim,Delta,Max_iter, ...
    sim.parameters('a_a')); % log-likelihoods found with the UPDATED a_a value

%% COLOURS FOR THE PLOT 

colours={[0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880],...
    [0.4940 0.1840 0.5560],...
    [0.9290 0.6940 0.1250],...
    [0.8500 0.3250 0.0980],...
    [0 0.4470 0.7410]};

%% FIGURE 2 A

Fa = figure('Position',[0 0 385 280]);
set(Fa, 'defaultAxesFontSize', 9)
set(Fa, 'defaultLineLineWidth', 1)

hold on
colourind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        colourind=colourind+1;
        last_nutr_qual=data.xdata(i,1);
    end
    plot(data.ydata(i,1),data.ydata(i,2),'o','Color',colours{colourind},'LineWidth',1) % real data
    plot(ymodel(i,1),ymodel(i,2),'+','Color',colours{colourind},'MarkerSize',8,'LineWidth',1.25) % model predictions

end
ylabel('\phi_r, ribosomal mass fraction','FontName','Arial');
xlabel('\lambda, growth rate [1/h]','FontName','Arial')
 
%% ADD lines for 1ST GROWTH LAW FITS

% group data by nutrient quality
xs_1={[],[],[],[]};
ys_1={[],[],[],[]};
nutrind=1;
chlorind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        nutrind=nutrind+1;
        chlorind=1;
        last_nutr_qual=data.xdata(i,1);
    end
    xs_1{chlorind}(end+1)=ymodel(i,1);
    ys_1{chlorind}(end+1)=ymodel(i,2);
    chlorind=chlorind+1;
end

% make linear fits
fit_coeffs=zeros([nutrind 2]);
for chlorind=1:size(xs_1,2)
    if(size(xs_1{chlorind},2)>=2)
        linfit=polyfit(xs_1{chlorind},ys_1{chlorind},1);
        fit_coeffs(chlorind,1)=linfit(1);
        fit_coeffs(chlorind,2)=linfit(2);
    end
end

% plot
dashings={'--','-.',':','-'};
for chlorind=4:(-1):1
    if(fit_coeffs(chlorind,1)~=0)
        xpoints=linspace(0,xs_1{chlorind}(end)*1.1,100); % points for which we plot the linear fit
        ypoints=polyval(fit_coeffs(chlorind,:),xpoints);
        plot(xpoints,ypoints,'Color','k','LineStyle',dashings{chlorind},'LineWidth',0.5)
    end
end

%% ADD lines for 2ND GROWTH LAW FITS

% group data by nutrient quality
xs_2={[]};
ys_2={[]};
nutrind=1;
chlorind=1;
last_nutr_qual=data.xdata(1,1);
for i=1:size(data.xdata,1)
    if(data.xdata(i,1)~=last_nutr_qual)
        nutrind=nutrind+1;
        chlorind=1;
        last_nutr_qual=data.xdata(i,1);
        xs_2{nutrind}=[];
        ys_2{nutrind}=[];
    end
    xs_2{nutrind}(chlorind)=ymodel(i,1);
    ys_2{nutrind}(chlorind)=ymodel(i,2);
    chlorind=chlorind+1;
end

% make linear fits - only based on points with up to 4 uM chloramphenicol
fit_coeffs=zeros([nutrind 2]);
for nutrind=1:size(xs_2,2)
    if(size(xs_2{nutrind},2)==2 || size(xs_2{nutrind},2)==3)
        linfit=polyfit(xs_2{nutrind},ys_2{nutrind},1);
        fit_coeffs(nutrind,1)=linfit(1);
        fit_coeffs(nutrind,2)=linfit(2);
    elseif(size(xs_2{nutrind},2)>3)
        linfit=polyfit(xs_2{nutrind}(1:3),ys_2{nutrind}(1:3),1);
        fit_coeffs(nutrind,1)=linfit(1);
        fit_coeffs(nutrind,2)=linfit(2);
    end
end

% plot
for nutrind=flip(1:size(xs_2,2))
    if(fit_coeffs(nutrind,1)~=0)
        xpoints=linspace(0,xs_2{nutrind}(1)*1.1,100); % points for which we plot the linear fit
        ypoints=polyval(fit_coeffs(nutrind,:),xpoints);
        plot(xpoints,ypoints,'Color',colours{nutrind},'LineWidth',0.5)
    end
end

%% SETTINGS of the plot
xlim([0 1.8])
ylim([0 0.46])

axis square
grid on
box on
hold off