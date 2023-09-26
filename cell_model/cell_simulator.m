%% cell_simulator.m
% Matlab class enabling simulations of the host cell. 

% The expression of synthetic circuits can be simulated by loading the
% 'heterologous genes' and 'external input' modules (see het_modules and
% ext_inputs folders). Remember to PUSH (obj.push_het()) the associated 
% modules' parameters into the main framework every time you alter them.

%%

classdef cell_simulator
    
    properties (SetAccess = public)
        % VARIABLES USED IN SIMULATIONS
        x0; % initial condition
        t; % time of simulation
        x; % state of the system
        
        % DESCRIBE THE SYSTEM
        init_conditions; % initial conditions according to which x0 is defined when simulation starts
        parameters; % parameters describing the host cell

        % HETEROLOGOUS GENES
        het; % object describing heterologous genes
        num_het=0; % number of heterologous genes
        num_misc=0; % number of miscellaneous species modelled

        % EXTERNAL INPUTS
        ext; % object describing external inputs
        num_ext=0; % number of external inputs
        
        % DEFAULT SETTINGS FOR THE SIMULATOR
        opt = odeset('RelTol', 1e-12, 'AbsTol', 1e-16); % integration tolerances
        tf = 100; % time frame over which we conduct the simulation
        form=cell_formulae; % formulae for rate and activation functions

        % HYBRID SIMULATOR
        S; % stoichiometry matix for hybrid simulations of heterologous gene expression
        num_stoch_reactions; % number of possible stochastic reactions
        em_timestep=0.01; % Eulare-Maruyama simulation timestep (h)
        em_save_every_n_steps=1000;
    end
    
    methods (Access = public)
        %% CONSTRUCTOR
        function obj = cell_simulator(tf)
            % if non-default simulation time suggested, use it
            if nargin == 1
                obj.tf = tf;
            end

            % set default parameters and initial conditions for the host cell
            obj=obj.set_default_parameters();
            obj=obj.set_default_init_conditions();
            
            % set up heterologous genes and external inputs
            obj=obj.load_heterologous_and_external('no_het','no_ext'); % none by default

            % push parameters and initial conditions of het. system into main framework
            obj=obj.push_het();
        end
        

        % SET default parameters (defined in cell_params.m)
        function obj=set_default_parameters(obj)
            obj.parameters=cell_params();
        end
        

        % SET default initial conditions (defined in cell_init_conds.m)
        function obj=set_default_init_conditions(obj)
            obj.init_conditions=cell_init_conds(obj.parameters);
        end
        
        %% HETEROLOGOUS GENE AND EXTERNAL INPUT MODULE MANAGEMENT
        % LOAD heterologous gene and external input modules
        function obj = load_heterologous_and_external(obj,Sys,ext_sys)
            % LOAD HETEROLOGOUS GENES MODULE
            % access the relevant file
            addpath(genpath([pwd, filesep,'het_modules']));
            het_func=str2func(Sys);

            % class file describing a synthetic gene system
            obj.het=het_func();

            % get number of heterologous genes
            obj.num_het=size(obj.het.names,2);

            % get number of miscellaneous species
            obj.num_misc=size(obj.het.misc_names,2);

            % LOAD EXTERNAL INPUTS MODULE
            % access the relevant file
            addpath(genpath([pwd, filesep,'ext_modules']));
            ext_func=str2func(ext_sys);

            % class file describing a synthetic gene system
            obj.ext=ext_func();

            % get number of heterologous genes
            obj.num_ext=size(obj.ext.name,2);

            % push parameters
            obj=obj.push_het();

            % CHECK COMPAITIBILITY
            % is het gene module compatible with external input module?
            if(size(obj.ext.compatible_hets,2)~=0)
                het_in_compatible=false;
                for i=1:size(obj.ext.compatible_hets,2)
                    if strcmp(obj.het.module_name,obj.ext.compatible_hets{i})
                        het_in_compatible=true;
                        break;
                    end
                end
            else
                het_in_compatible=true;
            end

            % does external input module allow the het gene module to work?
            if(size(obj.het.prerequisite_exts,2)~=0)
                ext_in_prerequisite=false;
                for i=1:size(obj.het.prerequisite_exts,2)
                    if strcmp(obj.ext.module_name,obj.het.prerequisite_exts{i})
                        ext_in_prerequisite=true;
                        break;
                    end
                end
            else
                ext_in_prerequisite=true;
            end

            % check compaitibility
            if ~(het_in_compatible && ext_in_prerequisite)
                disp('Incompatible modules! Expect errors')
            end
        end
        

        % PUSH parameters and initial conditions of het. system into main framework
        function obj=push_het(obj)
            obj=obj.push_het_parameters();
            obj=obj.push_het_init_conditions();
        end


        % PUSH heterologous gene parameter values into main framework
        function obj = push_het_parameters(obj)
            % add heterologous genes if there are any
            if(obj.num_het>0)
                for key=keys(obj.het.parameters)
                    obj.parameters(key{1})=obj.het.parameters(key{1});
                end
            end
        end
        

        % PUSH initial conditions for heterologous genes into main framework
        function obj = push_het_init_conditions(obj)
            % add heterologous genes if there are any
            if(obj.num_het>0)
                for key=keys(obj.het.init_conditions)
                    obj.init_conditions(key{1})=obj.het.init_conditions(key{1});
                end
            end
        end
        
        %% SIMULATION
        % CALL the simulator, save the outcome
        function obj = simulate_model(obj)
            obj = obj.set_x0; % set initial condition
            [obj.t, obj.x] = ode15s(@obj.ss_model, [0, obj.tf], [obj.x0], obj.opt);
        end
        

        % DEFINE initial condition according to obj.init_conditions
        function obj = set_x0(obj)
            % NATIVE GENES
            obj.x0 = [
                      % mRNAs;
                      obj.init_conditions('m_a'); % metabolic gene transcripts
                      obj.init_conditions('m_r'); % ribosomal gene transcripts

                      % proteins
                      obj.init_conditions('p_a'); % metabolic proteins
                      obj.init_conditions('R'); % non-inactivated ribosomes

                      % tRNAs
                      obj.init_conditions('tc'); % charged
                      obj.init_conditions('tu'); % uncharged

                      % free ribosomes inactivated by chloramphenicol
                      obj.init_conditions('Bcm');

                      % culture medium'S nutrient quality and chloramphenicol concentration
                      obj.init_conditions('s'); % nutrient quality
                      obj.init_conditions('h'); % chloramphenicol levels
                      ];

            % ...ADD HETEROLOGOUS GENES IF THERE ARE ANY
            if(obj.num_het>0)
                x0_het=zeros(2*obj.num_het,1); % initialise
                for i=1:obj.num_het
                    % mRNA
                    x0_het(i)=obj.init_conditions(['m_',obj.het.names{i}]);
                    % protein
                    x0_het(i+obj.num_het)=obj.init_conditions(['p_',obj.het.names{i}]);
                end
                obj.x0=[obj.x0;x0_het]; % concantenate
            end

            % ...ADD MISCELLANEOUS SPECIES IF THERE ARE ANY
            if(obj.num_misc>0)
                x0_misc=zeros(obj.num_misc,1); % initialise
                for i=1:obj.num_misc
                    x0_misc(i)=obj.init_conditions(obj.het.misc_names{i});
                end
                obj.x0=[obj.x0;x0_misc]; % concantenate
            end
        end
        

        % ODEs
        function dxdt = ss_model(obj, t, x)
            % denote obj. parameters as par for convenience
            par = obj.parameters;
            
            % give the state vector entries meaningful names
            m_a = x(1); % metabolic gene mRNA
            m_r = x(2); % ribosomal gene mRNA
            p_a = x(3); % metabolic proteins
            R = x(4); % non-inactivated ribosomes
            tc = x(5); % charged tRNAs
            tu = x(6); % uncharged tRNAs
            Bcm = x(7); % inactivated ribosomes
            s = x(8); % nutrient quality (constant)
            h = x(9); % chloramphenicol concentration (constant)
            x_het=x(10:(9+2*obj.num_het+obj.num_misc)); % heterologous genes and miscellaneous synthetic species

            % CALCULATE PHYSIOLOGICAL VARIABLES
            % translation elongation rate
            e=obj.form.e(par,tc);

            % ribosome inactivation rate due to chloramphenicol
            kcmh=par('kcm').*h;

            % ribosome dissociation constants
            k_a=obj.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=obj.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

            % heterologous genes rib. dissoc. constants
            k_het=ones(obj.num_het,1); % initialise with default value 1
            if(obj.num_het>0)
                for i=1:obj.num_het
                    k_het(i)=obj.form.k(e,...
                    obj.parameters(['k+_',obj.het.names{i}]),...
                    obj.parameters(['k-_',obj.het.names{i}]),...
                    obj.parameters(['n_',obj.het.names{i}]),...
                    kcmh);
                end
            end

            T=tc./tu... % ratio of charged to uncharged tRNAs 
                .*(1-par('is_fixed_T'))+par('fixed_T').*par('is_fixed_T'); % OR a fixed value (to enable comparison with flux-parity regulation)
            D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:obj.num_het)./k_het))./...
                (1-par('phi_q')); % denominator in ribosome competition calculations
            B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

            nu=obj.form.nu(par,tu,s); % tRNA charging rate

            l=obj.form.l(par,e,B); % growth rate

            psi=obj.form.psi(par,T); % tRNA synthesis rate - MUST BE SCALED BY GROWTH RATE

            % GET RNAP ACTIVITY
            rnap_act=l;

            % GET EXTERNAL INPUT
            ext_inp=obj.ext.input(x,t);

            % GET RATE OF EFFECTIVE NUTR. QUAL. CHANGE (for upshifts)
            if(par('is_upshift')==1)
                dsdt = (par('s_postshift') - s) * ...
                    (e./par('n_a')).*(m_a./k_a./D).*R./p_a;
            else
                dsdt = 0;
            end

            % DEFINE DX/DT FOR...
            % ...THE HOST CELL
            dxdt = [
                    % mRNAs
                    rnap_act.*par('c_a').*par('a_a')-(par('b_a')+l).*m_a-kcmh.*(m_a./k_a./D).*R;
                    rnap_act.*obj.form.F_r(par,T).*par('c_r').*par('a_r')-(par('b_r')+l).*m_r-kcmh.*(m_r./k_r./D).*R;
                    % ,metabolic protein a
                    (e./par('n_a')).*(m_a./k_a./D).*R-l.*p_a;
                    % ribosomes
                    (e./par('n_r')).*(m_r./k_r./D).*R-l.*R-kcmh.*B;
                    % tRNAs
                    nu.*p_a-l.*tc-e.*B;
                    psi*l-l.*tu-nu.*p_a+e.*B;
                    % ribosomes inactivated by chloramphenicol
                    kcmh.*B-l.*Bcm;
                    % nutrient quality
                    dsdt;
                    % chloramphenicol concentration assumed constant
                    0;
                    ];

            % ...HETEROLOGOUS GENES
            if(obj.num_het>0)
                dxdt_het=zeros(2*obj.num_het,1); % initialise
                % calculate
                for i=1:obj.num_het
                    % mRNA
                    dxdt_het(i)=rnap_act.*obj.het.regulation(obj.het.names{i},x,ext_inp)...
                        .*par(['c_',obj.het.names{i}]).*par(['a_',obj.het.names{i}])...
                        -(par(['b_',obj.het.names{i}])+l).*x_het(i)...
                        +obj.het.extra_m_term(obj.het.names{i},x,ext_inp)...
                        -kcmh.*(x_het(i)./k_het(i)./D).*R;

                    % protein
                    dxdt_het(i+obj.num_het)=(e./par(['n_',obj.het.names{i}])).*(x_het(i)./k_het(i)./D).*R...
                        -l.*x_het(i+obj.num_het)+...
                        obj.het.extra_p_term(obj.het.names{i},x,ext_inp);
                end
                dxdt=[dxdt;dxdt_het]; % concantenate
            end

            % ...MISCELLANEOUS SPECIES
            if(obj.num_misc>0)
                dxdt=[dxdt;obj.het.misc_ode(t,x,ext_inp,l)];
            end
        end

        %% STOCHASTIC SIMULATION - CHEMICAL LANGEVIN EQUATION
        
        % SIMULATE the system using the Euler-Maruyama method
        function obj = simulate_model_cle(obj)
            obj = obj.set_x0; % set initial condition
            obj =obj.generate_stoichiometry_matrix; % generate the stoichiometry matrix
            
            obj.t=0:(obj.em_timestep*obj.em_save_every_n_steps):obj.tf; % define the time axis
            obj.x=zeros(size(obj.x0,1),size(obj.t,2)); % initialise the trajectory
            obj.x(:,1)=obj.x0; % set the initial condition
            sqrt_timestep=sqrt(obj.em_timestep); % pre-calculkate the sqyuare root of time step for speed

            % simulate!
            t=obj.t(1);
            x=obj.x(:,1);
            record_cntr=2;
            for i=1:(size(0:obj.em_timestep:obj.tf,2)-1)
                % find the deterministic change
                det=obj.ss_model(t,x)*obj.em_timestep; 

                % find the stochastic change
                v=obj.generate_propensity_vector(t,x);
                disp(i)
                dw=normrnd(0,sqrt_timestep,size(v));
                stoch=obj.S*(sqrt(v).*dw);
                x=x+det+stoch;
                t=t+obj.em_timestep;

                % saving the simulation
                if(rem(i,obj.em_save_every_n_steps)==0)
                    obj.x(:,record_cntr)=x;
                    obj.t(record_cntr)=t;
                    disp(obj.t(record_cntr))
                    record_cntr=record_cntr+1;
                end
            end

            % transpose the records for consistecy with deterministic simulations
            obj.t=obj.t.';
            obj.x=obj.x.';
        end

        % CLE ODEs 
        function dxdt = ss_model_cle(obj, t, x)
            % DETERMINISTIC COMPONENT
            dxdt=obj.ss_model(t,x); % deterministic component

            % STOCHASTIC COMPONENT
            v=obj.generate_propensity_vector(t,x); % get the propensity vector for all stochastic reactions
            w=normrnd(0,1,size(v)); % get Gaussian noise
            dxdt=dxdt+(obj.S*(sqrt(v).*w)); % add the stochastic component
        end

        % generating the stoichiometry matrix for heterologous genes (also find number of stochastic reactions)
        function obj = generate_stoichiometry_matrix(obj)
            par=obj.parameters;
            
            % find total number of reactions that can occur
            obj.num_stoch_reactions = (2+obj.num_het)*3+... % synthesis/degradation/dilution of mRNA - 2 native genes + heterologous
                obj.num_het*2+... % synthesis/dilution of protein - 2 native genes + heterologous
                4; % tRNA aminoacylation, aa-tRNA dilution, uncraged tRNA syntesis, uncharged tRNA dilution
            if(strcmp(obj.het.module_name,'pi_controller'))
                obj.num_stoch_reactions=obj.num_stoch_reactions+3; % mutual annihiulation in AIF motif; degradation and dilution of the resultant complex
            end

            % initialise the stoichiometry matrix with zeros
            obj.S=zeros(9+2*obj.num_het+obj.num_misc, ... % number of species in the cell - INCLUDING HOST CELL VARIABLES
                obj.num_stoch_reactions);

            reaction_cntr=1; % start the reaction counter

            % DEFINE STOICHIOMETRIES FOR NATIVE GENES

            % mRNA ODEs
            for i=1:2 % metabolic and ribosomal genes
                obj.S(i,reaction_cntr)=1; % mRNA synthesis
                reaction_cntr=reaction_cntr+1;
                obj.S(i,reaction_cntr)=-1; % mRNA degradation
                reaction_cntr=reaction_cntr+1;
                obj.S(i,reaction_cntr)=-1; % mRNA dilution
                reaction_cntr=reaction_cntr+1;
            end

            % protein ODEs
            for i=1:2 % metabolic and ribosomal genes
                obj.S(2+i,reaction_cntr)=1; % protein synthesis
                obj.S(5,reaction_cntr)=-par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (-tc)
                obj.S(6,reaction_cntr)=par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (+tu)
                reaction_cntr=reaction_cntr+1;
                obj.S(2+i,reaction_cntr)=-1; % protein dilution
                reaction_cntr=reaction_cntr+1;
            end

            % tRNA ODEs
            obj.S(5,reaction_cntr)=1; % tRNA aminoacylation
            reaction_cntr=reaction_cntr+1;
            obj.S(5,reaction_cntr)=-1; % aa-tRNA dilution
            reaction_cntr=reaction_cntr+1;
            obj.S(6,reaction_cntr)=1; % (uncharged) tRNA synthesis
            reaction_cntr=reaction_cntr+1;
            obj.S(6,reaction_cntr)=-1; % (uncharged) tRNA dilution
            reaction_cntr=reaction_cntr+1;
            
            % DEFINE STOICHIOMETRIES FOR HETEROGLOUS GENES
            % mRNA - reactions common for all genes
            for i=1:obj.num_het
                obj.S(9+i,reaction_cntr)=1; % mRNA synthesis
                reaction_cntr=reaction_cntr+1;
                obj.S(9+i,reaction_cntr)=-1; % mRNA degradation
                reaction_cntr=reaction_cntr+1;
                obj.S(9+i,reaction_cntr)=-1; % mRNA dilution
                reaction_cntr=reaction_cntr+1;
            end

            % mRNA annihilation in AIF motif
            if(strcmp(obj.het.module_name,'pi_controller'))
                obj.S(9+2,reaction_cntr)=-1; obj.S(9+3,reaction_cntr)=-1; obj.S(9+2*obj.num_het+1,reaction_cntr)=1; % mutual annihilation of m_act and m_anti, formation of bound complex
                reaction_cntr=reaction_cntr+1;
            end

            % proteins - reactions for all genes
            for i=1:obj.num_het
                obj.S(9+obj.num_het+i,reaction_cntr)=1; % protein synthesis
                obj.S(5,reaction_cntr)=-par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (-tc)
                obj.S(6,reaction_cntr)=par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (+tu)
                reaction_cntr=reaction_cntr+1;
                obj.S(9+obj.num_het+i,reaction_cntr)=-1; % protein dilution
                reaction_cntr=reaction_cntr+1;
            end

            % m_anti-m_act bound complex dilution and degradation in AIF motif
            if(strcmp(obj.het.module_name,'pi_controller'))
                obj.S(9+2*obj.num_het+1,reaction_cntr)=-1; % degradation
                reaction_cntr=reaction_cntr+1;
                obj.S(9+2*obj.num_het+1,reaction_cntr)=-1; % dilution
                reaction_cntr=reaction_cntr+1;
            end
        end

        % generating the propensity vector at a given time
        function v=generate_propensity_vector(obj,t,x)
            % denote obj. parameters as par for convenience
            par = obj.parameters;
            
            % give the state vector entries meaningful names
            m_a = x(1); % metabolic gene mRNA
            m_r = x(2); % ribosomal gene mRNA
            p_a = x(3); % metabolic proteins
            R = x(4); % non-inactivated ribosomes
            tc = x(5); % charged tRNAs
            tu = x(6); % uncharged tRNAs
            Bcm = x(7); % inactivated ribosomes
            s = x(8); % nutrient quality (constant)
            h = x(9); % chloramphenicol concentration (constant)
            x_het=x(10:(9+2*obj.num_het+obj.num_misc)); % heterologous genes and miscellaneous synthetic species

            % CALCULATE PHYSIOLOGICAL VARIABLES
            % translation elongation rate
            e=obj.form.e(par,tc);

            % ribosome inactivation rate due to chloramphenicol
            kcmh=par('kcm').*h;

            % ribosome dissociation constants
            k_a=obj.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
            k_r=obj.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);

            % heterologous genes rib. dissoc. constants
            k_het=ones( size(x_het(1:obj.num_het)) ); % initialise with default value 1
            if(obj.num_het>0)
                for i=1:obj.num_het
                    k_het(i)=obj.form.k(e,...
                    obj.parameters(['k+_',obj.het.names{i}]),...
                    obj.parameters(['k-_',obj.het.names{i}]),...
                    obj.parameters(['n_',obj.het.names{i}]),...
                    kcmh);
                end
            end

            T=tc./tu... % ratio of charged to uncharged tRNAs 
                .*(1-par('is_fixed_T'))+par('fixed_T').*par('is_fixed_T'); % OR a fixed value (to enable comparison with flux-parity regulation)
            D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:obj.num_het)./k_het))./...
                (1-par('phi_q')); % denominator in ribosome competition calculations
            B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)

            nu=obj.form.nu(par,tu,s); % tRNA charging rate

            l=obj.form.l(par,e,B); % growth rate

            psi=obj.form.psi(par,T); % tRNA synthesis rate - MUST BE SCALED BY GROWTH RATE

            % GET RNAP ACTIVITY
            rnap_act=l;

            % GET EXTERNAL INPUT
            ext_inp=obj.ext.input(x,t);

            % initialise the propensity vector
            v=zeros(obj.num_stoch_reactions,1);
            
            % DEFINE PROPENSITIES FOR NATIVE GENES
            v=zeros(obj.num_stoch_reactions,1); % initialise the propensity vector
            reaction_cntr=1; % start the reaction counter

            % metabolic mRNA
            v(reaction_cntr) = rnap_act.*1.*par('c_a').*par('a_a'); % mRNA synthesis
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = par('b_a').*m_a; % mRNA degradation
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*m_a; % mRNA dilution
            reaction_cntr=reaction_cntr+1;
            % ribosomal mRNA
            v(reaction_cntr) = rnap_act.*obj.form.F_r(par,T).*par('c_r').*par('a_a'); % mRNA synthesis
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = par('b_a').*m_r; % mRNA degradation
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*m_r; % mRNA dilution
            reaction_cntr=reaction_cntr+1;

            % metabolic proteins
            v(reaction_cntr) = e./par('n_a').*(m_a./k_a./D).*R; % protein synthesis
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*p_a; % protein dilution
            reaction_cntr=reaction_cntr+1;
            % ribosomes
            v(reaction_cntr) = e./par('n_r').*(m_r./k_r./D).*R; % protein synthesis
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*R; % protein dilution
            reaction_cntr=reaction_cntr+1;

            % charged tRNAs
            v(reaction_cntr) = nu.*p_a; % tRNA aminoacylation
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*tc; % aa-tRNA dilution
            reaction_cntr=reaction_cntr+1;
            % uncharged tRNAs
            v(reaction_cntr) = psi.*l; % (uncharged) tRNA synthesis
            reaction_cntr=reaction_cntr+1;
            v(reaction_cntr) = l.*tu; % (uncharged) tRNA dilution
            reaction_cntr=reaction_cntr+1;

            % DEFINE PROPENSITIES FOR HETEROLOGOUS GENES

            % mRNA - reactions common for all genes
            for i=1:obj.num_het
                v(reaction_cntr) = rnap_act.*obj.het.regulation(obj.het.names{i},x,ext_inp)...
                        .*par(['c_',obj.het.names{i}]).*par(['a_',obj.het.names{i}]); % mRNA synthesis
                reaction_cntr=reaction_cntr+1;
                v(reaction_cntr) = par(['b_',obj.het.names{i}]).*x_het(i); % mRNA degradation
                reaction_cntr=reaction_cntr+1;
                v(reaction_cntr) = l.*x_het(i); % mRNA dilution
                reaction_cntr=reaction_cntr+1;
            end

            % mRNA annihilation in AIF motif
            if(strcmp(obj.het.module_name,'pi_controller'))
                v(reaction_cntr) = par('kb_anti').*x_het(2).*x_het(3); % mutual annihilation of m_act and m_anti, formation of bound complex
                reaction_cntr=reaction_cntr+1;
            end

            % proteins - reactions for all genes
            for i=1:obj.num_het
                v(reaction_cntr) = (e./par(['n_',obj.het.names{i}])).*(x_het(i)./k_het(i)./D).*R; % protein synthesis
                reaction_cntr=reaction_cntr+1;
                v(reaction_cntr) = l.*x_het(i+obj.num_het); % protein dilution
                reaction_cntr=reaction_cntr+1;
            end

            % m_anti-m_act bound complex dilution in AIF motif
            if(strcmp(obj.het.module_name,'pi_controller'))
                v(reaction_cntr) = par('b_bound').*x_het(2*obj.num_het+1); % degradation
                reaction_cntr=reaction_cntr+1;
                v(reaction_cntr) = l.*x_het(2*obj.num_het+1); % dilution
                reaction_cntr=reaction_cntr+1;
            end
        end
        
        % check the condition for discrete events happening
        function [value,isterminal,direction]=check_de_condition(obj,t,x)
            value=x(end)-obj.de_condition; % the last element of x is the integral of discrete event propensities over time
            isterminal=1; % need to halt ODE integration upon encountering a discrete event
            direction=+1; % the value function increases with time until it reaches zero
        end

%         %% HYBRID TAU-LEAP SIMULATION - HOST CELL VARIABLES DETERMINISTIC, HETEROLOGOUS CIRUIT STOCHASTIC
%         % hybrid simulator with tau-leaping
%         function obj = simulate_model_hybrid_tauleap(obj)
%             % set initial condition
%             obj = obj.set_x0; 
%             obj.x0=[obj.x0(1:9); ...        % host cell variables continuous and deterministic
%                    round(obj.x0(10:end))];   % heterologous circuit variables discrete and stochastic => must be integer
% 
%             % add the discrete event condition check to ODE options
%             detopt=obj.opt; % back up deterministic simulation options
%             obj.opt=odeset(obj.opt,'Events', @obj.check_de_condition);
%             
%             % genereate the stocihiometry matrix for heterologous genes
%             obj=obj.generate_Stoichiometry_matrix();
% 
%             % start at t=0
%             t=0;
% 
%             % initialise time and state records
%             obj.t=[0]; obj.x=[obj.x0];
% 
%             % simulate
%             while (t_start<obj.tf)
%                 disp(t_start)
%                                 
%                 % if we didn't just run into the end of the simulation overall
%                 if(t<obj.tf)
%                     % find number of times stochastic reactions happened over tau hours
%                     tau=obj.tau_dt;
%                     v=obj.generate_het_propensity_vector(obj,t,x);
%                     num_het_reactions=poissrnd(tau./v);
%                     dx=num_het_reactions*obj.S; % get the change in x
%                     [t_over_tau, x_over_tau] = ode15s(@obj.ss_model, [t, t+tau], [obj.x(end,:)], obj.opt);
% 
% 
%                     % if the updated x would turn out to be negative, decrease the time step
%                     while(any(obj.x(end,:)<dx))
%                         tau=tau/2;
%                         v=obj.generate_het_propensity_vector(obj,t,x);
%                         num_het_reactions=poissrnd(tau./v);
%                         dx=num_het_reactions*obj.S; % get the change in x
%                     end
% 
%                     % evolve the deterministic part of the system (native genes)
%                     [t_over_tau, x_over_tau] = ode15s(@obj.ss_model, [t, t+tau], [obj.x(end,:)], obj.opt);
% 
%                     % find which discrete event happened
%                     v=obj.generate_het_propensity_vector(t_until_de(end),x_until_de(end,:)); % get propensities of all reactions
%                     random_draw=rand*sum(v); % draw a random variable between 0 and 1, scale by total propensity
%                     which_reaction=find(random_draw>[0;cumsum(v)],1); % find which reaction occured according to the cdf and random draw
%                     
%                     % update x according to which reaction happened
%                     x_until_de(end,:)=x_until_de(end-1,:)+[obj.S(:,which_reaction);0].';
%                 end
% 
%                 % append the trajectory to our previous results
%                 obj.t=[obj.t; t_until_de];
%                 obj.x=[obj.x; x_until_de(1:end-1,:)];
%                 
%                 % next time, start where we finished
%                 t_start=t_until_de(end);                
%             end
% 
%             % undo the addition of a discrete event condtion check
%             obj.opt=detopt;
%         end
% 
%         % generating a stoichiometry matrix for native genes
%         function obj = generate_Stoichiometry_matrix(obj)
%             par=obj.parameters;
% 
%             % initialise the stoichiometry matrix with zeros
%             obj.S=zeros(9+2*obj.num_het+obj.num_misc, ... % number of species in the cell - INCLUDING HOST CELL VARIABLES
%                 2*3+2*2+4); % mRNA synthesis/degradation/dilution; protein synthesis/dilution; tRNA charging, synthesis, dilution of charged and uncharged tRNAs
% 
%             reaction_cntr=1; % start the reaction counter
% 
%             % mRNA ODEs
%             for i=1:2 % metabolic and ribosomal genes
%                 obj.S(i,reaction_cntr)=1; % mRNA synthesis
%                 reaction_cntr=reaction_cntr+1;
%                 obj.S(i,reaction_cntr)=-1; % mRNA degradation
%                 reaction_cntr=reaction_cntr+1;
%                 obj.S(i,reaction_cntr)=-1; % mRNA dilution
%                 reaction_cntr=reaction_cntr+1;
%             end
% 
%             % protein ODEs
%             for i=1:2 % metabolic and ribosomal genes
%                 obj.S(2+i,reaction_cntr)=1; % protein synthesis
%                 obj.S(5,reaction_cntr)=-par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (-tc)
%                 obj.S(6,reaction_cntr)=par(['n_',obj.het.names{i}]); % includes tRNA unchraging during translation (+tu)
%                 reaction_cntr=reaction_cntr+1;
%                 obj.S(2+i,reaction_cntr)=-1; % protein dilution
%                 reaction_cntr=reaction_cntr+1;
%             end
% 
%             % tRNA ODEs
%             obj.S(5,reaction_cntr)=1; % tRNA aminoacylation
%             reaction_cntr=reaction_cntr+1;
%             obj.S(5,reaction_cntr)=-1; % aa-tRNA dilution
%             reaction_cntr=reaction_cntr+1;
%             obj.S(6,reaction_cntr)=1; % (uncharged) tRNA synthesis
%             reaction_cntr=reaction_cntr+1;
%             obj.S(6,reaction_cntr)=-1; % (uncharged) tRNA dilution
%             reaction_cntr=reaction_cntr+1;
%         end
%         
%         % generating a propensity vector for native genes
%         function v=generate_nat_propensity_vector(obj,t,x)
%             % denote obj. parameters as par for convenience
%             par = obj.parameters;
%             
%             % give the state vector entries meaningful names
%             m_a = x(1); % metabolic gene mRNA
%             m_r = x(2); % ribosomal gene mRNA
%             p_a = x(3); % metabolic proteins
%             R = x(4); % non-inactivated ribosomes
%             tc = x(5); % charged tRNAs
%             tu = x(6); % uncharged tRNAs
%             Bcm = x(7); % inactivated ribosomes
%             s = x(8); % nutrient quality (constant)
%             h = x(9); % chloramphenicol concentration (constant)
%             x_het=x(10:(9+2*obj.num_het+obj.num_misc)); % heterologous genes and miscellaneous synthetic species
% 
%             % CALCULATE PHYSIOLOGICAL VARIABLES
%             % translation elongation rate
%             e=obj.form.e(par,tc);
% 
%             % ribosome inactivation rate due to chloramphenicol
%             kcmh=par('kcm').*h;
% 
%             % ribosome dissociation constants
%             k_a=obj.form.k(e,par('k+_a'),par('k-_a'),par('n_a'),kcmh);
%             k_r=obj.form.k(e,par('k+_r'),par('k-_r'),par('n_r'),kcmh);
% 
%             % heterologous genes rib. dissoc. constants
%             k_het=ones( size(x_het(1:obj.num_het)) ); % initialise with default value 1
%             if(obj.num_het>0)
%                 for i=1:obj.num_het
%                     k_het(i)=obj.form.k(e,...
%                     obj.parameters(['k+_',obj.het.names{i}]),...
%                     obj.parameters(['k-_',obj.het.names{i}]),...
%                     obj.parameters(['n_',obj.het.names{i}]),...
%                     kcmh);
%                 end
%             end
% 
%             T=tc./tu... % ratio of charged to uncharged tRNAs 
%                 .*(1-par('is_fixed_T'))+par('fixed_T').*par('is_fixed_T'); % OR a fixed value (to enable comparison with flux-parity regulation)
%             D=1+(m_a./k_a+m_r./k_r+sum(x_het(1:obj.num_het)./k_het))./...
%                 (1-par('phi_q')); % denominator in ribosome competition calculations
%             B=R.*(1-1./D); % actively translating ribosomes (inc. those translating housekeeping genes)
% 
%             nu=obj.form.nu(par,tu,s); % tRNA charging rate
% 
%             l=obj.form.l(par,e,B); % growth rate
% 
%             psi=obj.form.psi(par,T); % tRNA synthesis rate - MUST BE SCALED BY GROWTH RATE
% 
%             % GET RNAP ACTIVITY
%             rnap_act=l;
% 
%             % GET EXTERNAL INPUT
%             ext_inp=obj.ext.input(x,t);
%             
%             % DEFINE PRPENSITY VECTOR V FOR NATIVE GENES
%             v=zeros(obj.num_stoch_reactions,1); % initialise the propensity vector
%             reaction_cntr=1; % start the reaction counter
% 
%             % metabolic mRNA
%             v(reaction_cntr) = rnap_act.*1.*par('c_a').*par('a_a'); % mRNA synthesis
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = par('b_a').*m_a; % mRNA degradation
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*m_a; % mRNA dilution
%             reaction_cntr=reaction_cntr+1;
%             % ribosomal mRNA
%             v(reaction_cntr) = rnap_act.*obj.form.F_r(par,T).*par('c_r').*par('a_a'); % mRNA synthesis
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = par('b_a').*m_r; % mRNA degradation
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*m_r; % mRNA dilution
%             reaction_cntr=reaction_cntr+1;
% 
%             % metabolic proteins
%             v(reaction_cntr) = e./par('n_a').*(m_a./k_a./D).*R; % protein synthesis
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*p_a; % protein dilution
%             reaction_cntr=reaction_cntr+1;
%             % ribosomes
%             v(reaction_cntr) = e./par('n_r').*(m_r./k_r./D).*R; % protein synthesis
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*R; % protein dilution
%             reaction_cntr=reaction_cntr+1;
% 
%             % charged tRNAs
%             v(reaction_cntr) = nu.*p_a; % tRNA aminoacylation
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*tc; % aa-tRNA dilution
%             reaction_cntr=reaction_cntr+1;
%             % uncharged tRNAs
%             v(reaction_cntr) = psi.*l; % (uncharged) tRNA synthesis
%             reaction_cntr=reaction_cntr+1;
%             v(reaction_cntr) = l.*tu; % (uncharged) tRNA dilution
%             reaction_cntr=reaction_cntr+1;
%         end

        %% VARIED
        % PLOT simulation results
        function plot_simulation(obj,species,plot_type)
            % PLOT FOR THE HOST CELL
            if strcmp(species,'native')
                % PROTEIN COMPOISITION BY MASS
                if strcmp(plot_type,'protein masses')
                    % Protein massess (aa)
                    str_var = {'q','a','R','B_{cm}','het'};
                    amino_q=obj.parameters('M')*ones(size(obj.t))*obj.parameters('phi_q');
                    amino_a=obj.x(:,3)*obj.parameters('n_a');
                    amino_r=obj.x(:,4)*obj.parameters('n_r');
                    amino_Bcm=obj.x(:,7)*obj.parameters('n_r');
                    amino_het=zeros(size(obj.t)); % initialise masses of heterologous proteins
                    if(obj.num_het~=0)
                        x_het=obj.x(:,10: (9+2*obj.num_het) );
                        for i=1:obj.num_het
                            amino_het=sum((amino_het+x_het(:,(obj.num_het+i))...
                                *obj.parameters(['n_',obj.het.names{i}])),2);
                        end
                    end
                    Faa = figure('Position',[0 0 540 480]);
                    set(Faa, 'defaultLineLineWidth', 2)
                    set(Faa, 'defaultAxesFontSize', 16)
                    hold all
                    patch([obj.t; flip(obj.t)],[amino_q+amino_a+amino_r+amino_Bcm+amino_het; flip(amino_a+amino_r+amino_Bcm+amino_het)], [0.75, 0.75, 0.75])
                    patch([obj.t; flip(obj.t)],[amino_a+amino_r+amino_Bcm+amino_het; flip(amino_r+amino_Bcm)], [0.9290, 0.6940, 0.1250])
                    patch([obj.t; flip(obj.t)],[amino_r+amino_Bcm+amino_het; flip(amino_Bcm+amino_het)], [0.4940, 0.1840, 0.5560])
                    patch([obj.t; flip(obj.t)],[amino_Bcm+amino_het; flip(amino_het)], [0.4660 0.6740 0.1880])
                    patch([obj.t; flip(obj.t)],[amino_het; zeros(size(obj.t))], [0 0.4470 0.7410])
        
                    xlabel('Time (h)');
                    ylabel('Protein mass (aa)');
                    legend(str_var{[1,2,3,4,5]});
                    title('Prot. mass (aa) evolution');
                    % ylim([0 obj.parameters('M')])
                    hold off
                
                % CONCENTRATIONS OF MRNAs, PROTEINS, TRNAs
                elseif strcmp(plot_type,'concentrations')
                    Fm = figure('Position',[0 0 1000 800]);
                    set(Fm, 'defaultLineLineWidth', 2)
                    set(Fm, 'defaultAxesFontSize', 14)
                    
                    % mRNA levels
                    subplot(2,2,1)
                    hold on
                    plot(obj.t,obj.x(:,1),'Color',[0.9290, 0.6940, 0.1250]);
                    plot(obj.t,obj.x(:,2),'Color',[0.4940, 0.1840, 0.5560]);
                    xlabel('Time (h)');
                    ylabel('mRNA concentration (nM)');
                    legend('a','r');
                    title('mRNA levels');
                    hold off
                    
                    % protein levels
                    subplot(2,2,2)
                    hold on
                    plot(obj.t,obj.x(:,3),'Color',[0.9290, 0.6940, 0.1250]);
                    plot(obj.t,obj.x(:,4),'Color',[0.4940, 0.1840, 0.5560]);
                    xlabel('Time (h)');
                    ylabel('Protein concentration (nM)');
                    legend('a','r');
                    title('protein levels');
                    hold off
                    
                    % tRNA levels
                    subplot(2,2,3)
                    plot(obj.t,obj.x(:,5),obj.t,obj.x(:,6))
                    xlabel('Time (h)');
                    ylabel('tRNA concentration (nM)');
                    legend('charged','uncharged');
                    title('tRNA charging');
                end

            % PLOT FOR HETEROLOGOUS PROTEINS
            elseif strcmp(species,'heterologous')
                % colour scheme for heterologous proteins
                colours=[[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [1 0 1]; [0 1 0]];% colours for plots
                % TOTAL MASSES OF PROTEINS IN CELL
                if strcmp(plot_type,'protein masses')
                    x_het=obj.x(:,10: (9+2*obj.num_het) );

                    Faa_het = figure('Position',[0 0 540 480]);
                    hold all
                    amino=zeros(size(obj.t,1),obj.num_het);

                    % record masses evolution
                    for i=1:size(obj.het.names,2)
                        amino(:,i)=obj.parameters(['n_',obj.het.names{i}])*x_het(:,obj.num_het+i);
                    end
                    % make plots
                    for i=1:obj.num_het
                        if(i~=obj.num_het)
                            patch([obj.t; flip(obj.t)],[sum(amino(:,i:end),2);...
                                flip(sum(amino(:,(i+1):end),2))], colours(rem(i-1,size(colours,1))+1,:))
                        else
                            patch([obj.t; flip(obj.t)],[amino(:,end); zeros(size(obj.t))], colours(rem(i-1,size(colours,1))+1,:))
                        end
                    end
                    % format plot
                    xlabel('Time (h)');
                    ylabel('Protein mass (aa)');
                    legend(obj.het.names);
                    title('Prot. mass (aa) evolution');
                    hold off
                
                % MRNA AND PROTEIN CONCENTRATIONS
                elseif strcmp(plot_type,'concentrations')
                    Fc_het = figure('Position',[0 0 900 400]);
                    x_het=obj.x(:,10: (9+2*obj.num_het) );
                    
                    % protein levels
                    subplot(1,2,1)
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,x_het(:,(i+obj.num_het)),'Color',colours(rem(i-1,size(colours,1))+1,:));
                    end
                    xlabel('Time (h)');
                    ylabel('protein concentration (nM)');
                    legend(obj.het.names);
                    title('Protein levels');
                    hold off
    
                    % mRNA levels
                    subplot(1,2,2)
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,x_het(:,i),'Color',colours(rem(i-1,size(colours,1))+1,:));
                    end
                    xlabel('Time (h)');
                    ylabel('mRNA concentration (nM)');
                    legend(obj.het.names);
                    title('mRNA levels');
                    hold off
                
                % VALUES OF GENE TRANSCRIPTION REGULATION FUNCTIONS
                elseif strcmp(plot_type,'regulation')
                    F_reg = figure('Position',[0 0 500 400]);
                    % get external inputs
                    ext_inps=zeros(size(obj.t,1),1); % initialise
                    if ~strcmp(obj.ext.module_name,'no_ext') % in case of no ext. input, no need to calculate it
                        for i=1:size(obj.t,1)
                            ext_inps(i)=obj.ext.input(obj.x(i),obj.t(i));
                        end
                    end

                    % calculate regulation function values
                    F_het = zeros([size(obj.t,1) obj.num_het]);
                    for i=1:size(obj.t,1)
                        for j=1:obj.num_het
                            F_het(i,j) = obj.het.regulation(obj.het.names{j},obj.x(i,:),ext_inps(i));
                        end
                    end
                    
                    % plot
                    hold on
                    for i=1:obj.num_het
                        plot(obj.t,F_het(:,i),'Color',colours(rem(i-1,size(colours,1))+1,:))
                    end
                    xlabel('Time (h)');
                    ylabel('F_x');
                    legend(obj.het.names);
                    title('Transcription regulation');
                    ylim([0 1])
                    hold off
                    
                end
             elseif strcmp(species,'miscellaneous')
                Fmisc = figure('Position',[0 0 500 400]);

                hold on
                for i=1:obj.num_misc
                    plot(obj.t,obj.x(:,9+obj.num_het*2+i))
                end
                xlabel('Time (h)');
                ylabel('Concentration, nM');
                legend(obj.het.misc_names);
                title('Miscellaneous species');
                hold off
            elseif strcmp(species,'external')
                Fext = figure('Position',[0 0 500 400]);

                if(obj.num_ext>0)
                    ext_inps=zeros(size(obj.t,1),1); % initialise
                    for i=1:size(obj.t,1)
                        ext_inps(i)=obj.ext.input(obj.x(i),obj.t(i));
                    end
                end
                
                plot(obj.t,ext_inps,'k');
                xlabel('Time (h)');
                ylabel('External input');
                legend(obj.ext.name);
                title('External input over time');
                hold off
            end
        end
    

        % RETURN total heterologous mRNA transcription (gene conc. times transcription rate) for a given state of the system
        function tht = total_heterologous_transcription(obj, x, t)
            tht=zeros(1,size(x,1));
            if(obj.num_het>0)
                ext_inp=obj.ext.input(x,t);
                for i=1:obj.num_het
                    F_i=obj.het.regulation(obj.het.names{i},x,ext_inp);
                    tht=tht+F_i.*obj.parameters(['a_',obj.het.names{i}]).*obj.parameters(['c_',obj.het.names{i}]);
                end
            end
        end
    end
end