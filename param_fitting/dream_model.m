%% dream_model.m
% Log-likelihood for our model with a given set of parameter values being
% fitted - for parameter fitting using DREAM

%%
function loglike=dream_model(theta)
    %% LOAD experimental data
    persistent data % experimental data retained after first function call
    
    if(isempty(data))
        dataset = readmatrix('data/growth_rib_fit_notext.csv');

        % nutrient qualities are equally log-spaced points
        nutr_quals=logspace(log10(0.08),log10(0.5),6);
        
        % get inputs: nutrient quality and h; get outputs: l and rib mass frac
        data.xdata=[]; % initialise inputs array
        data.ydata=[]; % initialise outputs array

        for i = 1:size(dataset,1)
            if(dataset(i,1)>0.3)
                % inputs
                nutr_qual = nutr_quals(fix((i-1)/5)+1); % records start from worst nutrient quality
                h = dataset(i,4)*1000; % all h values for same nutr quality same go one after another. Convert to nM from uM!
                data.xdata=[data.xdata; [nutr_qual,h]];
                
                l = dataset(i,1); % growth rate (1/h)
                phi_r = dataset(i,3); % ribosome mass fraction
                data.ydata=[data.ydata; [l,phi_r]];
            end
        end
    end

    % standard measurement errors
    l_stdev=0.04467;
    phir_stdev=0.018976;
    
    %% DEFINE the simulator
    persistent sim % enough to initialise the simulator just once
    if(isempty(sim))
        sim=cell_simulator; % initialise simulator
        
        % settings
        sim.tf = 12; % single integraton step timeframe
        sim.opt = odeset('reltol',1.e-6,'abstol',1.e-6); % more lenient integration tolerances for speed
    end

    % parameters for getting steady state
    Delta = 0.001; % threshold that determines if we're in steady state
    Max_iter = 4; % maximum no. iterations (checking if SS reached over first 750 h)
    
    %% GET model predictions
    ymodel=dream_modelfun(theta,data.xdata,sim,Delta,Max_iter,3.89e5); % using crude estimate 3.89e5 for a_a

    %% FIND log likelihood

    % get sum of squares
    sos = sum((ymodel - data.ydata).^2);

    % for testing - print the results
    disp(['SOS=',num2str(sos)])

    % get logarithm of mutually independent normal distributions
    loglike=-log(l_stdev.*sqrt(2.*pi))-0.5.*sos(1)./(l_stdev.^2)... % growth rate
        -log(phir_stdev.*sqrt(2.*pi))-0.5.*sos(2)./(phir_stdev.^2); % ribosomal mass fraction

end