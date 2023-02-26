%% figS1b.m

% MCMC FITTING
% Figure S1: b

% Demonstrate that the fit's accuracy is relatively unaffected over a wide
% range of ribosomal and metabolic gene promoter strengths as long as their
% ratio stays the same

%% CLEAR parameters, add paths of all files

addpath(genpath('..'))

clear
close all

%% LOAD experimental data
 
dataset = readmatrix('data/growth_rib_fit_notext.csv');

% nutrient qualities are equally log-spaced points
nutr_quals=logspace(log10(0.08),log10(0.5),6);

% get inputs: nutrient quality and h; get outputs: l and rib mass frac
xdata=[]; % initialise inputs array
ydata=[]; % initialise outputs array
for i = 1:size(dataset,1)
    if(dataset(i,1)>0.3)
        % inputs
        nutr_qual = nutr_quals(fix((i-1)/5)+1); % records start from worst nutrient quality
        h = dataset(i,4)*1000; % all h values for same nutr quality same go one after another. Convert to nM from uM!
        xdata=[xdata; [nutr_qual,h]];
        
        l = dataset(i,1); % growth rate (1/h)
        phi_r = dataset(i,3); % ribosome mass fraction
        ydata=[ydata; [l,phi_r]];
    end
end

% standard measurement errors
l_stdev=0.04467;
phir_stdev=0.018976;

%% DEFINE the simulator
sim=cell_simulator; % initialise simulator

% settings
sim.tf = 10; % single integraton step timeframe
sim.opt = odeset('reltol',1.e-6,'abstol',1.e-9); % more lenient integration tolerances for speed

% parameters for getting steady state
Delta = 0.1; % threshold that determines if we're in steady state
Max_iter = 50; % maximum no. iterations (checking if SS reached over first 500 h)

%% DEFINE a_a and a_r value to check
a_as=logspace(log10(sim.parameters('a_a')/10),log10(sim.parameters('a_a')*30),9);
a_rs=logspace(log10(sim.parameters('a_r')/10),log10(sim.parameters('a_r')*30),9);

% initialise array of adaptation errors
loglikes=zeros(size(a_as,2),size(a_rs,2));

% save the explored ranges
mcmcsave.a_as=a_as;
mcmcsave.a_rs=a_rs;
save('mcmcsave.mat','mcmcsave');

%% SET the starting point for simulation (not 1 if we're restarting a saved simulation batch)

mcmcsave.i_a_first=1;

%% RUN the simulations and calculate log-likelihoods

for i_a=mcmcsave.i_a_first:size(a_as,2)
    for i_r=1:size(a_rs,2)
        % reset parameters and initial condition
        sim=sim.set_default_parameters();
        sim=sim.set_default_init_conditions();
        
        % change fitted parameters to current values
        sim.parameters('a_a') = a_as(i_a); % metabolic prot. transcription rate (/h)
        sim.parameters('a_r') = a_rs(i_r); % ribosome transcription rate (/h) - rescaled!
        
        % find steady state values for every input
        ymodel=zeros([size(xdata,1) 2]); % initialise array of avlues predicted by the model
        for i=1:size(xdata,1)
            %disp(i)
            % change parameter values to relevant inputs
            sim.init_conditions('s')=xdata(i,1);
            sim.init_conditions('h')=xdata(i,2);
            
            % evaluate steady state
            ss=get_steady(sim,Delta,Max_iter); % evaluate steady state value
        
            % get growth rate and ribosome mass fraction
            par=sim.parameters;
            m_a = ss(1);
            m_r = ss(2);
            p_a = ss(3);
            R = ss(4);
            tc = ss(5);
            tu = ss(6);
            Bcm = ss(7);
            s = ss(8);
            h = ss(9);
            ss_het=ss(10 : (9+2*sim.num_het) ).';
        
            e=sim.form.e(par,tc); % translation elongation rate
            kcmh=par('kcm').*h; % ribosome inactivation rate due to chloramphenicol
        
            % ribosome dissociation constants
            k_a=sim.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=sim.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
            k_het=ones(sim.num_het,1);
        
            % ribosome dissociation constants for heterologous genes
            if(sim.num_het>0)
                for j=1:sim.num_het
                    k_het(i)=sim.form.k(e,...
                    sim.parameters(['k+_',sim.het.names{j}]),...
                    sim.parameters(['k-_',sim.het.names{j}]),...
                    sim.parameters(['n_',sim.het.names{j}]),...
                    kcmh);
                end
            end
        
            D=1+(m_a./k_a+m_r./k_r+sum(ss_het(1:sim.num_het)./k_het))./...
                (1-par('phi_q')); % denominator in ribosome competition calculations
            B=R.*(1-1./D); % actively translating ribosomes - INCLUDING HOUSEKEEPING GENES
        
            % growth rate
            l=sim.form.l(par,e,B);
            
            % RECORD!
            ymodel(i,1)=l; % record growth rate
            ymodel(i,2)=(R+Bcm).*sim.parameters('n_r')./sim.parameters('M'); % record ribosome mass fraction
        end

        sos=sum((ymodel - ydata).^2); % find Sum Of Squared errors

        loglikes(i_a,i_r)=-log(l_stdev.*sqrt(2.*pi))-0.5.*sos(1)./(l_stdev.^2)... % growth rate
        -log(phir_stdev.*sqrt(2.*pi))-0.5.*sos(2)./(phir_stdev.^2); % ribosomal mass fraction
    end

    % checkpoint save
    mcmcsave.i_a_first=i_a+1;
    mcmcsave.loglikes=loglikes;
    save('mcmcsave.mat','mcmcsave')
    disp([num2str(i_a/size(a_as,2)*100),' % completed and saved'])
end

%% PLOT

Fg = figure('Position',[0 0 300 300]);
set(Fg, 'defaultAxesFontSize', 9)
set(Fg, 'defaultLineLineWidth', 1.25)

% with aif controller
hmap=heatmap(flip(loglikes,1),'ColorMap', jet(100));
% display labels
x_text_labels=string(a_as);
for i=1:size(a_as,2)
    midpoint=round(size(a_as,2)/2,0);
    if (i~=1 && i~=size(a_as,2) && i~=midpoint)
        x_text_labels(i)='';
    else
        if(a_as(i)~=0)
            % find which power of 10 to use in the engineering notation
            for power_of_ten=0:100
                if(10^power_of_ten>a_as(i))
                    break;
                end
            end
            power_of_ten=power_of_ten-1;

            % make axis label
            x_text_labels(i)=[num2str(round(a_as(i)/10^power_of_ten,2)),'e',num2str(power_of_ten)];
        else
            x_text_labels(i)='0';
        end
    end
end

a_rs_flip=flip(a_rs);
y_text_labels=string(a_rs_flip);
for i=1:size(a_rs_flip,2)
    midpoint=round(size(a_rs_flip,2)/2,0);
    if (i~=1 && i~=size(a_rs_flip,2) && i~=midpoint)
        y_text_labels(i)=' ';
    else
        if(a_rs_flip(i)~=0)
            % find which power of 10 to use in the engineering notation
            for power_of_ten=0:100
                if(10^power_of_ten>a_rs_flip(i))
                    break;
                end
            end
            power_of_ten=power_of_ten-1;

            % make axis label
            y_text_labels(i)=[num2str(round(a_rs_flip(i)/10^power_of_ten,2)),'e',num2str(power_of_ten)];
        else
            y_text_labels(i)='0';
        end
    end
end
hmap.XDisplayLabels = x_text_labels;
hmap.YDisplayLabels = y_text_labels;

ylabel('\alpha_r');
xlabel('\alpha_a');

caxis(hmap,[-3000 0]);
hmap.GridVisible = 'off';