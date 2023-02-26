%% rc_ecoli_modelfun.m
% retrieve steady-state growth rate and ribosomal mass fraction - for
% parameter fitting using DREAM

%%
function ymodel=dream_modelfun(theta,xdata,sim,Delta,Max_iter)
    % reset parameters and initial condition
    sim=sim.set_default_parameters();
    sim=sim.set_default_init_conditions();

    % change fitted parameters to current values
    sim.parameters('a_a') = 2.72e5; % metabolic prot. transcription rate (/h)
    sim.parameters('a_r') = sim.parameters('a_a').*theta(1); % ribosome transcription rate (/h) - rescaled!
    sim.parameters('nu_max') = theta(2); % max metabolic rate (/h)
    sim.parameters('K_e') = theta(3); % elongation rate Hill constant (nM)
    sim.parameters('K_nut') = theta(3); % tRNA charging rate Hill constant (nM)
    sim.parameters('kcm') = theta(4); % chloramphenical binding rate constant (/h/nM)
    
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
end