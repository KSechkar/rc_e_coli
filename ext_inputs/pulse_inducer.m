%% pulse_inducer.m
% Describes an external input (e.g. chemical or light stimulus applied to cells)
% Here, A PULSE IN CHEMICAL INDUCER CONCENTRATION

%%

classdef pulse_inducer
    properties (SetAccess = public)
        module_name='pulse_inducer'; % name of the external input module
        name; % name of input
        input_func_parameters; % parameters of the input functions
        
        % 'Compatible' heterologous gene expression modules (e.g. if administered 
        % signal is claculated using fluorescent readouts, only modules with F.P. expression can be used with it)
        % If this is left empty, the module is assumed compatible with any set of synthetic genes.
        compatible_hets; 
    end

    % regulatory functions
    methods (Access = public)
        function obj =pulse_inducer(obj)
            obj.module_name='pulse_inducer';
            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT INFO----------------------------------
            
            % SPECIFY INPUT NAMES HERE
            obj.name={'inducer'}; % chemical inducer

            % SPECIFY COMPATIBLE HETEROLOGOUS MODULES
            obj.compatible_hets={};


            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------


            % initialise map storing parameterss
            obj.input_func_parameters=containers.Map('KeyType', 'char', ...
                    'ValueType', 'double');

            % -------------------------------------------------------------
            % SPECIFY INPUT INFO--------------------------------------------

            % SPECIFY PARAMETERS HERE
            obj.input_func_parameters('inducer_base_level')=100; % baseline conc of inducer (nM)
            obj.input_func_parameters('pulse_value_prop')=1.25; % specify value during the pulse relative to baseline
            obj.input_func_parameters('pulse_start_time')=10; % specify start time of the pulse
            obj.input_func_parameters('pulse_duration')=1; % pulse duration (h)

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % applied input as a function of the system's state x and time t
        function inp = input(obj,x,t)

            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT FUNCTION------------------------------
            ifpar=obj.input_func_parameters;

            if (t>ifpar('pulse_start_time') && t<ifpar('pulse_start_time')+ifpar('pulse_duration'))
                inp=ifpar('inducer_base_level').*ifpar('pulse_value_prop'); % within the pulse time window, apply sepcified concentration
            else
                inp=ifpar('inducer_base_level'); % otherwise, apply baseline concentration
            end

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end

    end
end