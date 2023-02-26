%% heterologous_approx.m
% collection of functions for approximating the effects of heterologous
% gene expression.Required for Figure 3 scripts

% WARNING: these functions work only when the external input is unchanging,
% as we assume t=0 to calculate this input

%%

classdef heterologous_approx
   
    methods (Access = public)
        % setady state mRNA levels
        function m=ss_m(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                F, ... % gene of interest (GOI) ss regulation function value
                c, ... % GOI concentration
                a, ... % GOI transcription rate
                b, ... % GOI mRNA decay rate
                ss, ... % steady state of the system
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)

            % get growth rate
            l=obj.ss_l(plasmid_conc,ss,ss0,e0,sim);

            %  get mRNA level
            m=F.*c.*a./(l+b);
        end
        
        % steady state protein mass fraction
        function phi=ss_phi(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss, ... % steady state of the system
                F, ... % gene of interest (GOI) ss regulation function value
                c, ... % GOI concentration
                a, ... % GOI transcription rate
                kplus, ... % GOI mRNA-ribosome binding rate
                kminus, ... % GOI mRNA-ribosome unbinding rate
                n, ...  % GOI length in aa
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in an array
            par=sim.parameters;

            % GET the denominator for phi_i expression
            denominator=obj.denominator(plasmid_conc,ss,ss0,e0,sim);

            % get the phi value
            k = sim.form.k(e0,kplus,kminus,n,par('kcm').*ss(8));
            phi = F.*c.*a./k./denominator;
        end

        % steady state protein mass fraction OF ALL HETEROLOGOUS GENES
        function phi_het=ss_phi_het(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss, ... % steady state of the system
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in an array
            par=sim.parameters;
            ext_inp=sim.ext.input(ss,0);

            % GET the denominator for phi_i expression
            denominator=obj.denominator(plasmid_conc,ss,ss0,e0,sim);

            % get numerator (sum of all heterologous genes' contributions)
            numerator=0; % initialise numerator
            for i=1:sim.num_het
                k_i=sim.form.k(e0,par(['k+_',sim.het.names{i}]),...
                    par(['k-_',sim.het.names{i}]),par(['n_',sim.het.names{i}]),par('kcm').*ss(8));
                F_i=sim.het.regulation(sim.het.names{i},ss, ext_inp);
                numerator=numerator+F_i.*plasmid_conc.*par(['a_',sim.het.names{i}])./k_i;
            end

            % calculate the phi_het value
            phi_het=numerator./denominator;
        end

        % steady state growth rate
        function l=ss_l(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss, ... % steady state of the system
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)
            
            % store model parameters in an array
            par=sim.parameters;

            % get the l value
            F_r=sim.form.F_r(par,ss0(5)./ss0(6));
            phi_r=obj.ss_phi(plasmid_conc,ss,F_r,par('c_r'),par('a_r'),...
                par('k+_r'),par('k-_r'),par('n_r'),ss0,e0,sim);
            l=e0./par('n_r').*phi_r;
        end
        
        % denominator used in the expressions for phi and l
        function denom=denominator(obj,... % concentration of plasmid with heterolss_l_fullogous genes
                plasmid_conc,... % plasmid cocnentration
                ss,... % steady state of the system
                ss0,... % steady state without heterologous gene expression
                e0,... % steady state translation elongation rate without heterologous gene expression
                sim) % simulator (required to get parameters and regulatory functions for all genes)

            % store model parameters in an array
            par=sim.parameters;

            denom = 0; % initialise denominator
            kcmh=par('kcm').*ss(8);

            % add contribution of metabolic genes
            k_a=sim.form.k(e0,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            denom=denom+par('c_a').*par('a_a')./k_a;

            % add contribution of ribosomal gene
            k_r=sim.form.k(e0,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
            F_r=sim.form.F_r(par,ss0(5)./ss0(6));
            denom=denom+F_r.*par('c_r').*par('a_r')./k_r;

            ext_inp=sim.ext.input(ss,0);

            % add contributions of heterologous genes
            for i=1:sim.num_het
                k_i=sim.form.k(e0,par(['k+_',sim.het.names{i}]),...
                    par(['k-_',sim.het.names{i}]),par(['n_',sim.het.names{i}]),kcmh);
                F_i=sim.het.regulation(sim.het.names{i},ss,ext_inp); % right now, regulation assumed constant!
                denom=denom+F_i.*plasmid_conc.*par(['a_',sim.het.names{i}])./k_i;
            end

            % rescale to account for housekeeping genes
            denom=denom./(1-sim.parameters('phi_q'));
        end
           
        % total protein production rate constant for the population of cells at t=0
        function mu_het=ss_mu_het(obj, ...
                plasmid_conc, ... % concentration of plasmid with heterologous genes
                ss, ... % steady state of the system
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                l0, ... % steady state translation elongation rate without heterologous gene expression
                sim,... % simulator (required to get parameters and regulatory functions for all genes)
                delta) % death rate
            par=sim.parameters;
            phi_het=obj.ss_phi_het(plasmid_conc,ss,ss0,e0,sim);

            mu_het=par('M').*... % cell mass
                phi_het.*... % het prot mass fraction
                (l0.*(1-phi_het./(1-par('phi_q')))-delta); % growth rate-death rate
        end
        
        % heterologous protein fraction maximising mu_het
        function phi_het_max=phi_het_max(obj, ...
                l0, ... % steady state translation elongation rate without heterologous gene expression
                sim,... % simulator (required to get parameters and regulatory functions for all genes)
                delta) % death rate

            phi_het_max=0.5.*(1-delta./l0).*(1-sim.parameters('phi_q'));
        end

        % heterologous protein expression burden xi maximising mu_het
        function xi_max=xi_max(obj, ...
                ss, ... % steady state of the system
                ss0, ... % steady state without heterologous gene expression
                e0, ... % steady state translation elongation rate without heterologous gene expression
                l0, ... % steady state translation elongation rate without heterologous gene expression
                sim,... % simulator (required to get parameters and regulatory functions for all genes)
                delta) % death rate
            par=sim.parameters;
            
            % FIND LUMPED NATIVE GENE EXPRESSION PARAMETER
            sum_native=0; % initialise sum of native genes' lumped parameters
            kcmh=par('kcm').*ss(8); % get ribosome inactivation rate

            % add contribution of metabolic genes
            k_a=sim.form.k(e0,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            sum_native=sum_native+par('c_a').*par('a_a')./k_a;

            % add contribution of ribosomal gene
            k_r=sim.form.k(e0,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
            F_r=sim.form.F_r(par,ss0(5)./ss0(6));
            sum_native=sum_native+F_r.*par('c_r').*par('a_r')./k_r;

            % FIND OPTIMAL PROTEIN MASS FRACTION
            phi_het_max=obj.phi_het_max(l0,sim,delta);

            % FIND OPTIMAL BURDEN XI
            %xi_max=(phi_het_max/(1-phi_het_max-par('phi_q')))*sum_native; % calculate optimal burden
            xi_max=(1-delta/l0)/(1+delta/l0)*sum_native;
        end
    end

end