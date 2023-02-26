%% two_constit.m
% Describes heterologous genes expressed in the cell
% Here, TWO CONSTITUTIVE GENES

%%

classdef two_constit
    properties (SetAccess = public)
        module_name='two_constit'; % name of the heterologous gene module
        names; % names of all heterologous genes
        misc_names; % names of modelled miscellanous species (e.g. compound of interest that het. proteins synthesise)
        parameters; % parameters of heterologous genes
        init_conditions; % initial conditions for all modelled species

        % 'Compatible' external input signal modules (e.g.
        % optogenetic circuits only work with light signals)
        % If this is left empty, the module is assumed compatible with any set of synthetic genes.
        prerequisite_exts;
    end

    methods (Access = public)
        % CONSTRUCTOR
        function obj = two_constit(obj)

            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            obj.prerequisite_exts={}; % not affected by any external signals => irrelevant

            % SPECIFY GENE NAMES HERE
            obj.names={'xtra1','xtra2'};
            
            % SPECIFY MISCELLANEOUS SPECIES TO MODEL
            obj.misc_names={}; % no miscellaneous species

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------


             % initialise the map storing all genes' parameters and fill with default gene expression parameters
            obj.parameters=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');
             for i=1:size(obj.names,2)
                obj.parameters(['c_',obj.names{i}]) = 100; % copy no. (nM)
                obj.parameters(['a_',obj.names{i}]) = 100; % max transcription rate (/h)
                obj.parameters(['b_',obj.names{i}]) = 6; % mRNA decay rate (/h)
                obj.parameters(['k+_',obj.names{i}]) = 60; % ribosome binding rate (/h/nM)
                obj.parameters(['k-_',obj.names{i}]) = 60; % ribosome unbinding rate (/h)
                obj.parameters(['n_',obj.names{i}]) = 300; % protein length (aa)
             end

            % initialise the map storing initial conditions and fill with default values
            obj.init_conditions=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');
             for i=1:size(obj.names,2)
                obj.init_conditions(['m_',obj.names{i}]) = 0;
                obj.init_conditions(['p_',obj.names{i}]) = 0;
             end
             for i=1:size(obj.misc_names,2)
                obj.init_conditions(obj.misc_names{i}) = 0;
             end


            % -------------------------------------------------------------
            % SPECIFY GENE INFO--------------------------------------------

            % SPECIFY NON-DEFAULT PARAMETERS HERE
            % all parameters assumed to be default
            obj.parameters('a_xtra2')=200;

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % regulation function for heterologous gene transcription
        % for the gene called gene_name -
        % depends on the system's state x and external input ext_inp
        function F = regulation(obj,gene_name,x,ext_inp)
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            misc=x(10+size(obj.names,2)*2:end); % get miscellaneous species info

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'xtra1')
                F=1; % constitutive expression
            elseif strcmp(gene_name,'xtra2')
                F=1; % constitutive expression

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
            else
                F=1; % constitutive expression by default
            end
        end

        % any extra terms in mRNA ODEs apart from synthesis & growth dilution
        function extra_term=extra_m_term(obj,gene_name,x,ext_inp)
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            misc=x(10+size(obj.names,2)*2:end); % get miscellaneous species info

            % -------------------------------------------------------------
            % SPECIFY TERMS------------------------------------------------

            if strcmp(gene_name,'xtra_1')
                extra_term=0; % no extra terms in ODEs
            elseif strcmp(gene_name,'xtra_2')
                extra_term=0; % no extra terms in ODEs

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------

            else
                extra_term=0; % none by default
            end
        end

        % any extra terms in protein ODEs apart from synthesis & growth dilution
        function extra_term=extra_p_term(obj,gene_name,x,ext_inp)
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            misc=x(10+size(obj.names,2)*2:end); % get miscellaneous species info

            % -------------------------------------------------------------
            % SPECIFY REGULATORY FUNCTIONS---------------------------------

            if strcmp(gene_name,'xtra_1')
                extra_term=0; % no extra terms in ODEs
            elseif strcmp(gene_name,'xtra_2')
                extra_term=0; % no extra terms in ODEs

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
            
            else
                extra_term=0; % none by default
            end
        end

        % ODEs for miscellaneous species in the system
        function dxdt=misc_ode(obj,t,x,ext_inp,l)
            het_par=obj.parameters;
            x_het=x( 10 : 9+size(obj.names,2)*2 ); % get heterologous gene info
            misc=x(10+size(obj.names,2)*2:end); % get miscellaneous species info

            % -------------------------------------------------------------
            % SPECIFY  MISCELLANEOUS SPECIES' ODES-------------------------
            
            dxdt=[]; % no miscellaneous species

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
    end
end