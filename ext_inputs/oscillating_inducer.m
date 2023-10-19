%% oscillating_inducer.m
% Describes an external input (e.g. chemical or light stimulus applied to cells)
% Here, A SINE WAVE

%%

classdef oscillating_inducer
    properties (SetAccess = public)
        module_name='oscillating_inducer'; % name of the external input module
        name; % name of input
        input_func_parameters; % parameters of the input functions
        
        % 'Compatible' heterologous gene expression modules (e.g. if administered 
        % signal is claculated using fluorescent readouts, only modules with F.P. expression can be used with it)
        % If this is left empty, the module is assumed compatible with any set of synthetic genes.
        compatible_hets; 
    end

    % regulatory functions
    methods (Access = public)
        function obj =oscillating_inducer(obj)
            obj.module_name='oscillating_inducer';
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
            obj.input_func_parameters('wave_period')=2; % period of the sine wave (h)
            obj.input_func_parameters('wave_amplitude')=0.5; % amplitude of the sine wave
            obj.input_func_parameters('oscillation_start_time')=72; % start time of oscillations
             obj.input_func_parameters('phase_shift')=pi/2; % phase shift of the sine wave

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end
        
        % applied input as a function of the system's state x and time t
        function inp = input(obj,x,t)

            % -------------------------------------------------------------
            % SPECIFY EXTERNAL INPUT FUNCTION------------------------------
            ifpar=obj.input_func_parameters;

            if (t>ifpar('oscillation_start_time'))
                inp=ifpar('wave_amplitude').*(1+sin( ...
                    (t-ifpar('oscillation_start_time')).*(2.*pi./ifpar('wave_period'))-ifpar('phase_shift') ...
                    )); % within the pulse time window, apply sepcified concentration
            else
                inp=0; % otherwise, apply baseline concentration
            end

            % END OF USER SPEC---------------------------------------------
            % -------------------------------------------------------------
        end

    end
end