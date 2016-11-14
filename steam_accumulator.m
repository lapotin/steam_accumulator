classdef steam_accumulator < handle
    %STEAMACCUMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % verbosity - Determines the kind of debug statements sent to
        % stdout.
        verbosity = 0
        % i - Total loop count for the accumulator.
        i
        % max_iter - Pre-assignment size for vectors.
        max_iter
        % time_step - The time step of the assigned to each loop.
        time_step
        % tank_length - Length of the accumulator piping [m].
        tank_length
        % tank_volume - Volume of the accumulator piping [m^3].
        tank_volume
        % water_level - Water level in the tank [%].
        water_level
        % high_pressure - Design starting pressure for discharge [bar].
        high_pressure
        % low_pressure - Design low pressure for discharge [bar].
        low_pressure
        % high_quality - Design high steam quality after charging.
        high_quality
        % low_quality - Design low steam quality after discharging.
        low_quality
        % minimum_efficiency - Design minimum cycle efficiency.
        minimum_efficiency
        % minimum_water_mass - Design minimum water mass required at the
        % beginning of a discharge [kg].
        minimum_water_mass
        % volume_defect_detected - True if a volume_defect was detected in
        % any loop.
        volume_defect_detected = false
        % initial_total_enthalpy
        initial_total_enthalpy
        initial_steam_mass
        initial_water_mass
        initial_water_level
        prepared_for_evaluation = false
        p = []
        x = []
        v_1 = []
        v_2 = []
        rho_1 = []
        rho_2 = []
        rho_mixture = []
        v_mixture = []
        t_1 = []
        t_2 = []
        m_1 = []
        m_2 = []
        m_total = []
        vol_1 = []
        vol_2 = []
        vol_total = []
        vol_defect = []
        h_1 = []
        h_2 = []
        q_loss = []
        q_loss_1 = []
        q_loss_2 = []
        q_21 = []
        m_dot_1b = []
        m_dot_2b = []
        mh_dot_1b = []
        mh_dot_2b = []
        r = []
        m_c = []
        m_e = []
        m_dot_pt_1 = []
        m_dot_pt_2 = []
        dv1dh = []
        dv2dh = []
        dv1dp = []
        dv2dp = []
        term1 = []
        term2 = []
        term3 = []
        term4 = []
        term5 = []
        term6 = []
        dpdt = []
        dh_1dt = []
        dh_2dt = []
        time = []
        m_dot_turb = []
        eta = []
        m_in = []
        m_out = []
        loop_time = []
        x_1 = []
        x_2 = []
    end
    properties(Constant)
        TAU = 85.0  % relaxation time [s]
        EPSILON = 0.02  % used to disrupt equilibrium in temperature between
                        % liquid and steam phases when setting up initial
                        % conditions
        MPA_PER_BAR = 0.1   % MPa per bar [MPa/bar]
        BAR_PER_MPA = 10    % bar per MPa [bar/MPa]
        W_PER_KW = 1000 % watts per kilowatt [W/kW]
        KW_PER_W = 0.001    % kilowatts per watt [kW/W]
        MASS_TOLERANCE = 0.01   % tolerance used to determine if a mass
                                % imbalance is present [%]
        VOLUME_TOLERANCE = 0.05 % tolerance used to determine if a volume
                                % imbalance is present
        PIPE_RADIUS = 0.4064    % radius of the natural gas pipeline being used
                                % to construct the accumulator
        PIPE_THICKNESS = 0.15875    % thickness of the natural gas pipeline
        INSULATION_THICKNESS = 0.2032   % insulation thickness
        K_INSULATION = 0.079    % thermal conductivity of the insulation
        K_PIPE = 41.0   % thermal conductivity of the pipe
        H_AIR = 15.0    % heat transfer coefficient of the air
        PASCALS_PER_BAR = 1e5   % Pascals per bar [Pa/bar]
        BAR_PER_PASCAL = 1e-5   % bar per Pascal [bar/Pa]
        ATM_PRESSURE_IN_BAR = 1.01325   % atmospheric pressure [bar]
    end
    methods
        function obj = steam_accumulator(p0, x0, l_tank, dt, max_iter, verbosity)
            if nargin > 0
                obj.verbosity = verbosity;
                obj.i = 1;
                obj.max_iter = max_iter;
                obj.time_step = dt;
                obj.tank_length = l_tank;
                obj.initialize_arrays();
                obj.setup_initial_conditions(p0, x0);
            end
        end
        function value = get.tank_volume(obj)
            value = pi * (obj.PIPE_RADIUS^2) * obj.tank_length;
        end
        function value = get.initial_total_enthalpy(obj)
            value = obj.h_1(1) * obj.m_1(1) + obj.h_2(1) * obj.m_2(1);
        end
        function value = get.initial_water_mass(obj)
           value = obj.m_1(1); 
        end
        function value = get.initial_steam_mass(obj)
           value = obj.m_2(1); 
        end
        function value = get.initial_water_level(obj)
           value = obj.vol_1(1) / obj.vol_total(1); 
        end
        function [mass_charged, charge_time, final_quality, volume_defect] = charge(obj,...
                final_pressure, final_quality, charge_pressure,...
                charge_flowrate, time_step)
            mass_charged = 0.0;
            obj.time_step = time_step;
            if(obj.verbosity > 0)
                fprintf('Commencing charge to %.2f bar.\n', final_pressure);
            end
            start_time = tic;
            i_begin = obj.i;
            while obj.p(obj.i) < final_pressure
                delta_t = obj.time(obj.i) - obj.time(i_begin);
                flow_multiplier = delta_t / 600.0;
                if flow_multiplier > 1.0
                    flow_multiplier = 1.0;
                end
                charge_steam_flow = 0.0;
                charge_liquid_flow = 0.0;
                h_1_in = IAPWS_IF97('hL_p', charge_pressure * obj.MPA_PER_BAR);
                h_2_in = IAPWS_IF97('hV_p', charge_pressure * obj.MPA_PER_BAR);
                if obj.p(obj.i) < final_pressure
                    charge_steam_flow = charge_flowrate * flow_multiplier;
                end                              
                if(obj.verbosity > 1 && mod(obj.i, 100) == 0)
                    fprintf('%10s = %.4f.\n', 'Tank Level', obj.water_level(obj.i));
                    fprintf('Steam flow: %.2f kg/s.\n', charge_steam_flow);
                    fprintf('Liquid flow: %.2f kg/s.\n', charge_liquid_flow);
                end
                mass_charged = mass_charged + (charge_liquid_flow + charge_steam_flow) * obj.time_step;
                obj.increment_time_step(charge_liquid_flow, 0.0, charge_steam_flow, 0.0, h_1_in, h_2_in);
            end
            if(obj.verbosity > 0)
                fprintf('Charge function complete in %.2f seconds.\n', toc(start_time));
                fprintf('Time to charge %.2f hours.\n', (obj.i - i_begin) * obj.time_step / 3600);
            end
            charge_time = (obj.i - i_begin) * obj.time_step;
            final_quality = obj.x(obj.i);
            volume_defect = obj.volume_defect_detected;
        end
        function [mass_discharged, final_quality, volume_defect] = discharge(obj, power,...
                time, minimum_pressure, time_step)
            mass_discharged = 0.0;
            if(obj.verbosity > 0)
                fprintf('Discharging at %.1f MW for %.1f hour(s).\n', power, time / 3600);
            end
            start_time = tic;
            obj.time_step = time_step;
            i_begin = obj.i;
            while (obj.i - i_begin) * obj.time_step < time
                [obj.eta(obj.i), obj.m_dot_turb(obj.i)] = obj.run_turbine(obj.p(obj.i), 0.9, power, 5);
                mass_discharged = mass_discharged + obj.m_dot_turb(obj.i) * obj.time_step;
                obj.increment_time_step_test(0.0, 0.0, 0.0, obj.m_dot_turb(obj.i), 0.0, 0.0);
            end
            if(obj.verbosity > 0)
                fprintf('Discharge function complete in %.2f seconds.\n', toc(start_time));
                fprintf('Discharge complete (%.2f MW, %.2f hours).\n', power, (obj.i - i_begin) * obj.time_step / 3600);
            end
            final_quality = obj.x(obj.i);
            volume_defect = obj.volume_defect_detected;
        end
        function [mass_discharged, final_quality, volume_defect] = discharge_test(obj, power,...
                time, time_step)
            mass_discharged = 0.0;
            if(obj.verbosity > 0)
                fprintf('Discharging at %.1f MW for %.1f hour(s).\n', power, time / 3600);
            end
            start_time = tic;
            obj.time_step = time_step;
            i_begin = obj.i;
            while (obj.i - i_begin) * obj.time_step < time
                [obj.eta(obj.i), obj.m_dot_turb(obj.i), m_dot_1, h_1, m_dot_2, h_2] = obj.run_turbine_test(obj.p(obj.i), obj.h_2(obj.i), 0.9, power, 5);
                mass_discharged = mass_discharged + obj.m_dot_turb(obj.i) * obj.time_step - m_dot_1 * obj.time_step;
                obj.increment_time_step(m_dot_1, 0.0, 0.0, obj.m_dot_turb(obj.i), h_1, 0.0);
            end
            if(obj.verbosity > 0)
                fprintf('Discharge function complete in %.2f seconds.\n', toc(start_time));
                fprintf('Discharge complete (%.2f MW, %.2f hours).\n', power, (obj.i - i_begin) * obj.time_step / 3600);
            end
            final_quality = obj.x(obj.i);
            volume_defect = obj.volume_defect_detected;
        end
        function obj = wait(obj, time, time_step)
            obj.time_step = time_step;
            begin_i = obj.i;
            while (obj.i - begin_i) * obj.time_step < time
                obj.increment_time_step(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            end
        end
        function obj = get_plots(obj)
           time_plot = obj.time(1:obj.i);
           
           p_plot = obj.p(1:obj.i);
           figure(1);
           plot(time_plot, p_plot);
           hold on;
           xlabel('Time [s]');
           ylabel('Pressure [bar]');
           xlim([1, max(time_plot)]);
           
           x_plot = obj.x(1:obj.i);
           figure(2);
           plot(time_plot, x_plot * 100);
           hold on;
           xlabel('Time [s]');
           ylabel('Quality [%]');
           xlim([1, max(time_plot)]);
           
           t_1_plot = obj.t_1(1:obj.i);
           t_2_plot = obj.t_2(1:obj.i);
           figure(3);
           plot(time_plot, t_1_plot - 273.15, 'r')
           hold on;
           plot(time_plot, t_2_plot - 273.15)           
           xlabel('Time [s]');
           ylabel('Phase Temperature [C]');
           legend('Liquid', 'Steam');
           xlim([1, max(time_plot)]);
           
           vol_1_plot = obj.vol_1(1:obj.i);
           vol_2_plot = obj.vol_2(1:obj.i);
           figure(4);
           plot(time_plot, vol_1_plot, 'r');
           hold on;
           plot(time_plot, vol_2_plot);
           plot(time_plot, vol_1_plot + vol_2_plot, 'g');
           xlabel('Time [s]');
           ylabel('Phase Volume [m^3]');
           legend('Liquid', 'Steam', 'Total');
           xlim([1, max(time_plot)]);
           
           eta_plot = obj.eta(1:obj.i);
           eta_non_zeros = eta_plot(eta_plot~=0);
           average_eta = mean(eta_non_zeros);
           average_eta_plot = average_eta * ones(1, obj.i);
           figure(6);
           plot(time_plot, eta_plot * 100);
           hold on;
           plot(time_plot, average_eta_plot * 100, 'r');
           xlabel('Time [s]');
           ylabel('Efficiency [%]');
           xlim([1, max(time_plot)]);
           
           h_1_plot = obj.h_1(1:obj.i);
           m_1_plot = obj.m_1(1:obj.i);
           h_2_plot = obj.h_2(1:obj.i);
           m_2_plot = obj.m_2(1:obj.i);
           figure(7);
           plot(time_plot, times(h_1_plot, m_1_plot), 'r');
           hold on;
           plot(time_plot, times(h_2_plot, m_2_plot));
           plot(time_plot, times(h_1_plot, m_1_plot) + times(h_2_plot, m_2_plot), 'g');
           xlabel('Time [s]');
           ylabel('Phase Enthalpy [kJ]');
           legend('Liquid', 'Steam', 'Total');
           xlim([1, max(time_plot)]);
           
           figure(8);
           plot(time_plot, m_1_plot, 'r');
           hold on;
           plot(time_plot, m_2_plot);
           plot(time_plot, m_1_plot + m_2_plot, 'g');
           xlabel('Time [s]');
           ylabel('Phase Mass [kg]');
           legend('Liquid', 'Steam', 'Total');
           xlim([1, max(time_plot)]);
           
           m_dot_pt_1_plot = obj.m_dot_pt_1(1:obj.i);
           m_dot_pt_2_plot = obj.m_dot_pt_2(1:obj.i);
           figure(9);
           plot(time_plot, m_dot_pt_1_plot, 'r');
           hold on;
           plot(time_plot, m_dot_pt_2_plot);
           plot(time_plot, m_dot_pt_1_plot + m_dot_pt_2_plot, 'g');
           xlabel('Time [s]');
           ylabel('Mass Change from Evaporation/Condensation [kg/s]');
           legend('Liquid', 'Steam', 'Total');
           xlim([1, max(time_plot)]);
           
           water_level_plot = obj.water_level(1:obj.i);
           figure(10);
           plot(time_plot, water_level_plot);
           hold on;
           xlabel('Time [s]');
           ylabel('Tank Water Level [%]');
           xlim([1, max(time_plot)]);
           
           v_1_plot = obj.v_1(1:obj.i);
           v_2_plot = obj.v_2(1:obj.i);
           figure(11);
           plot(time_plot, v_1_plot, 'r');
           hold on;
           plot(time_plot, v_2_plot);
           plot(time_plot, v_1_plot + v_2_plot, 'g');
           xlabel('Time [s]');
           ylabel('Specific Volume [m^3/kg]');
           legend('Liquid', 'Steam', 'Total');
           xlim([1, max(time_plot)]);
           
           x_1_plot = obj.x_1(1:obj.i);
           x_2_plot = obj.x_2(1:obj.i);
           figure(12);
           plot(time_plot, x_1_plot, 'r');
           hold on;
           plot(time_plot, x_2_plot);
           xlabel('Time [s]');
           ylabel('Quality');
           legend('Liquid', 'Steam');
           xlim([1, max(time_plot)]);
           
           q_21_plot = obj.q_21(1:obj.i);
           figure(13);
           plot(time_plot, q_21_plot);
           xlabel('Time [s]');
           ylabel('Heat transfer [kJ/s]');
           xlim([1, max(time_plot)]);
           
        end
        function obj = reset(obj, pressure, quality)
            obj.i = 1;
            obj.initialize_arrays();
            obj.setup_initial_conditions(pressure, quality);
        end
        function obj = soft_reset(obj)
            if obj.i <= 1
               return; 
            end
            obj.p = obj.p(obj.i:end); 
            obj.x = obj.x(obj.i:end);
            obj.v_1 = obj.v_1(obj.i:end); 
            obj.v_2 = obj.v_2(obj.i:end);
            obj.rho_1 = obj.rho_1(obj.i:end);
            obj.rho_2 = obj.rho_2(obj.i:end);
            obj.rho_mixture = obj.rho_mixture(obj.i:end);
            obj.v_mixture =obj.v_mixture(obj.i:end); 
            obj.t_1 = obj.t_1(obj.i:end);
            obj.t_2 = obj.t_2(obj.i:end);
            obj.m_1 = obj.m_1(obj.i:end);
            obj.m_2 = obj.m_2(obj.i:end); 
            obj.m_total = obj.m_total(obj.i:end); 
            obj.vol_1 = obj.vol_1(obj.i:end); 
            obj.vol_2 = obj.vol_2(obj.i:end);
            obj.vol_total = obj.vol_total(obj.i:end);
            obj.vol_defect = obj.vol_defect(obj.i:end); 
            obj.h_1 = obj.h_1(obj.i:end);
            obj.h_2 = obj.h_2(obj.i:end); 
            obj.q_loss = obj.q_loss(obj.i:end); 
            obj.q_loss_1 = obj.q_loss_1(obj.i:end);
            obj.q_loss_2 = obj.q_loss_2(obj.i:end);
            obj.q_21 = obj.q_21(obj.i:end);
            obj.m_dot_1b = obj.m_dot_1b(obj.i:end); 
            obj.m_dot_2b = obj.m_dot_2b(obj.i:end); 
            obj.mh_dot_1b = obj.mh_dot_1b(obj.i:end);
            obj.mh_dot_2b = obj.mh_dot_2b(obj.i:end); 
            obj.r = obj.r(obj.i:end); 
            obj.m_c = obj.m_c(obj.i:end); 
            obj.m_e = obj.m_e(obj.i:end);
            obj.m_dot_pt_1 = obj.m_dot_pt_1(obj.i:end);
            obj.m_dot_pt_2 = obj.m_dot_pt_2(obj.i:end); 
            obj.dv1dh = obj.dv1dh(obj.i:end); 
            obj.dv2dh = obj.dv2dh(obj.i:end);  
            obj.dv1dp = obj.dv1dp(obj.i:end); 
            obj.dv2dp = obj.dv2dp(obj.i:end);  
            obj.term1 = obj.term1(obj.i:end);  
            obj.term2 = obj.term2(obj.i:end); 
            obj.term3 = obj.term3(obj.i:end);  
            obj.term4 = obj.term4(obj.i:end); 
            obj.term5 = obj.term5(obj.i:end);  
            obj.term6 = obj.term6(obj.i:end);  
            obj.dpdt = obj.dpdt(obj.i:end); 
            obj.dh_1dt = obj.dh_1dt(obj.i:end); 
            obj.dh_2dt = obj.dh_2dt(obj.i:end);
            obj.time = obj.time(obj.i:end);            
            obj.m_dot_turb = obj.m_dot_turb(obj.i:end); 
            obj.eta = obj.eta(obj.i:end);
            obj.m_in = obj.m_in(obj.i:end);
            obj.m_out = obj.m_out(obj.i:end);
            obj.loop_time = obj.loop_time(obj.i:end);
            obj.water_level = obj.water_level(obj.i:end);
            obj.i = 1;
            obj.time(obj.i) = 0.0;
            obj.time_step = 1.0;
        end
        function obj = prep_for_eval(obj, power, time, charge_pressure, charge_flowrate)
            for j = 1:3
                [mass_discharged, final_qualitym, volume_defect] = obj.discharge(power, time, obj.low_pressure, 1.0);
                [mass_charged, charge_time, final_quality, volume_defect] = obj.charge(obj.high_pressure, obj.high_quality, charge_pressure, charge_flowrate, 1.0);                
            end
            obj.soft_reset();
            obj.prepared_for_evaluation = true;
        end
        function [avg_mass_charged, avg_charge_time, avg_mass_discharged,...
                avg_min_quality, avg_max_quality, valid_model] =...
                evaluate_accumulator(obj, power, time, charge_pressure, charge_flowrate)
%             if obj.prepared_for_evaluation == false
%                 obj.prep_for_eval(power, time, charge_pressure, charge_flowrate);
%             end
            n = 1;
            mass_charged = zeros(1, n);
            mass_discharged = zeros(1, n);
            charge_time = zeros(1, n);
            volume_defect_d = zeros(1, n);
            volume_defect_c = zeros(1, n);
            max_quality = zeros(1, n);
            min_quality = zeros(1, n);
            parfor j = 1:n
                par_object = obj;
                [mass_discharged(j), min_quality(j), volume_defect_d(j)] = par_object.discharge_test(power, time, obj.time_step);
                [mass_charged(j), charge_time(j), max_quality(j), volume_defect_c(j)] = par_object.charge(obj.high_pressure, obj.high_quality, charge_pressure, charge_flowrate, obj.time_step);                 
            end
            avg_mass_charged = mean(mass_charged);
            avg_charge_time = mean(charge_time);
            avg_mass_discharged = mean(mass_discharged);
            avg_min_quality = mean(min_quality);
            avg_max_quality = mean(max_quality);
            valid_model = all(volume_defect_d == false) && all(volume_defect_c == false);
        end
    end
    methods (Access = private)
        function obj = initialize_arrays(obj)
            obj.p = zeros(1, obj.max_iter); 
            obj.x = zeros(1, obj.max_iter);
            obj.v_1 = zeros(1, obj.max_iter); 
            obj.v_2 = zeros(1, obj.max_iter);
            obj.rho_1 = zeros(1, obj.max_iter);
            obj.rho_2 = zeros(1, obj.max_iter);
            obj.rho_mixture = zeros(1, obj.max_iter);
            obj.v_mixture = zeros(1, obj.max_iter); 
            obj.t_1 = zeros(1, obj.max_iter); 
            obj.t_2 = zeros(1, obj.max_iter); 
            obj.m_1 = zeros(1, obj.max_iter); 
            obj.m_2 = zeros(1, obj.max_iter); 
            obj.m_total = zeros(1, obj.max_iter); 
            obj.vol_1 = zeros(1, obj.max_iter); 
            obj.vol_2 = zeros(1, obj.max_iter); 
            obj.vol_total = zeros(1, obj.max_iter);
            obj.vol_defect = zeros(1, obj.max_iter); 
            obj.h_1 = zeros(1, obj.max_iter); 
            obj.h_2 = zeros(1, obj.max_iter); 
            obj.q_loss = zeros(1, obj.max_iter); 
            obj.q_loss_1 = zeros(1, obj.max_iter);
            obj.q_loss_2 = zeros(1, obj.max_iter);
            obj.q_21 = zeros(1, obj.max_iter); 
            obj.m_dot_1b = zeros(1, obj.max_iter); 
            obj.m_dot_2b = zeros(1, obj.max_iter); 
            obj.mh_dot_1b = zeros(1, obj.max_iter);
            obj.mh_dot_2b = zeros(1, obj.max_iter); 
            obj.r = zeros(1, obj.max_iter); 
            obj.m_c = zeros(1, obj.max_iter); 
            obj.m_e = zeros(1, obj.max_iter); 
            obj.m_dot_pt_1 = zeros(1, obj.max_iter);
            obj.m_dot_pt_2 = zeros(1, obj.max_iter); 
            obj.dv1dh = zeros(1, obj.max_iter); 
            obj.dv2dh = zeros(1, obj.max_iter); 
            obj.dv1dp = zeros(1, obj.max_iter); 
            obj.dv2dp = zeros(1, obj.max_iter); 
            obj.term1 = zeros(1, obj.max_iter); 
            obj.term2 = zeros(1, obj.max_iter);
            obj.term3 = zeros(1, obj.max_iter); 
            obj.term4 = zeros(1, obj.max_iter); 
            obj.term5 = zeros(1, obj.max_iter); 
            obj.term6 = zeros(1, obj.max_iter); 
            obj.dpdt = zeros(1, obj.max_iter); 
            obj.dh_1dt = zeros(1, obj.max_iter);
            obj.dh_2dt = zeros(1, obj.max_iter);
            obj.time = zeros(1, obj.max_iter); 
            obj.m_dot_turb = zeros(1, obj.max_iter); 
            obj.eta = zeros(1, obj.max_iter); 
            obj.m_in = zeros(1, obj.max_iter); 
            obj.m_out = zeros(1, obj.max_iter);
            obj.loop_time = zeros(1, obj.max_iter);
            obj.water_level = zeros(1, obj.max_iter);
            obj.x_1 = zeros(1, obj.max_iter);
            obj.x_2 = zeros(1, obj.max_iter);
        end
        function obj = setup_initial_conditions(obj, p0, x0)
            if(obj.verbosity > 0)
                fprintf('Setting up initial conditions for p = %0.2f bar and x = %0.2f and length = %0.2f m...\n', p0, x0, obj.tank_length);
            end
            obj.time(obj.i) = 0.0;
            obj.p(obj.i) = p0;
            obj.x(obj.i) = x0;
            if(obj.verbosity > 1)
                fprintf('%10s = %10f\n', 'p', obj.p(obj.i));
                fprintf('%10s = %10f\n', 'x', obj.x(obj.i));
            end
            obj.t_1(obj.i) = IAPWS_IF97('Tsat_p', obj.p(obj.i) * obj.MPA_PER_BAR);
            obj.t_2(obj.i) = obj.t_1(obj.i);
            obj.h_1(obj.i) = IAPWS_IF97('h_pT', obj.p(obj.i) * obj.MPA_PER_BAR, obj.t_1(obj.i) - obj.EPSILON);
            obj.h_2(obj.i) = IAPWS_IF97('h_pT', obj.p(obj.i) * obj.MPA_PER_BAR, obj.t_2(obj.i) + obj.EPSILON);
            obj.x_1(obj.i) = 0.0;
            obj.x_2(obj.i) = 1.0;
%             obj.v_1(obj.i) = IAPWS_IF97('vL_p', obj.p(obj.i) * obj.MPA_PER_BAR);
%             obj.v_2(obj.i) = IAPWS_IF97('vV_p', obj.p(obj.i) * obj.MPA_PER_BAR);
            obj.v_1(obj.i) = IAPWS_IF97('v_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_1(obj.i));
            obj.v_2(obj.i) = IAPWS_IF97('v_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_2(obj.i));            
            obj.rho_1(obj.i) = 1 / obj.v_1(obj.i);
            obj.rho_2(obj.i) = 1 / obj.v_2(obj.i);
            obj.rho_mixture(obj.i) = 1 / (obj.x(obj.i) / obj.rho_2(obj.i) + (1 - obj.x(obj.i)) / obj.rho_1(obj.i));
            obj.v_mixture(obj.i) = 1 / obj.rho_mixture(obj.i);
            obj.m_total(obj.i) = obj.rho_mixture(obj.i) * obj.tank_volume;
            obj.m_1(obj.i) = obj.m_total(obj.i) * (1 - obj.x(obj.i));
            obj.m_2(obj.i) = obj.m_total(obj.i) * obj.x(obj.i);           
            obj.vol_1(obj.i) = obj.m_1(obj.i) / obj.rho_1(obj.i);
            obj.vol_2(obj.i) = obj.m_2(obj.i) / obj.rho_2(obj.i);
            obj.vol_total(obj.i) = obj.vol_1(obj.i) + obj.vol_2(obj.i);
            obj.vol_defect(obj.i) = abs(obj.vol_total(obj.i) - obj.tank_volume) / obj.tank_volume;
            if(obj.vol_defect(obj.i) > obj.VOLUME_TOLERANCE)
                obj.volume_defect_detected = true;
            end
            if(obj.verbosity > 0 && obj.vol_defect(obj.i) > obj.VOLUME_TOLERANCE)
                disp('Volume defect detected. Results unreliable!');
            end
            obj.water_level(obj.i) = obj.vol_1(obj.i) / obj.vol_total(obj.i);
            
        end
        function obj = increment_time_step(obj, m_dot_1_in, m_dot_1_out, m_dot_2_in, m_dot_2_out, h_1_in, h_2_in)
            start_time = tic;
            obj.m_in(obj.i) = (m_dot_1_in + m_dot_2_in) * obj.time_step;
            obj.m_out(obj.i) = (m_dot_1_out + m_dot_2_out) * obj.time_step;
            h_f = IAPWS_IF97('hL_p', obj.p(obj.i) * obj.MPA_PER_BAR);
            h_g = IAPWS_IF97('hV_p', obj.p(obj.i) * obj.MPA_PER_BAR);
            obj.r(obj.i) = h_g - h_f;
            obj.q_loss(obj.i) = obj.heat_loss();
            obj.q_loss_1(obj.i) = (obj.vol_1(obj.i) / obj.vol_total(obj.i)) * obj.q_loss(obj.i);
            obj.q_loss_2(obj.i) = (obj.vol_2(obj.i) / obj.vol_total(obj.i)) * obj.q_loss(obj.i);
            obj.q_21(obj.i) = obj.heat_transfer_21();
            obj.m_dot_1b(obj.i) = m_dot_1_in - m_dot_1_out;
            obj.m_dot_2b(obj.i) = m_dot_2_in - m_dot_2_out;
            obj.mh_dot_1b(obj.i) = m_dot_1_in * h_1_in - m_dot_1_out * obj.h_1(obj.i);
            obj.mh_dot_2b(obj.i) = m_dot_2_in * h_2_in - m_dot_2_out * obj.h_2(obj.i);
            obj.m_c(obj.i) = 0.0;
            obj.m_e(obj.i) = 0.0;
            if(obj.h_1(obj.i) > h_f)
                obj.m_e(obj.i) = obj.rho_1(obj.i) * obj.vol_1(obj.i) * (obj.h_1(obj.i) - h_f) / (obj.TAU * obj.r(obj.i));
            else
                obj.m_c(obj.i) = obj.rho_1(obj.i) * obj.vol_1(obj.i) * (h_f - obj.h_1(obj.i)) / (obj.TAU * obj.r(obj.i));
            end
%             if(obj.h_2(obj.i) < h_g)
%                 obj.m_c(obj.i) = obj.m_c(obj.i) + obj.rho_2(obj.i) * obj.vol_2(obj.i) * (h_g - obj.h_2(obj.i)) / (obj.TAU * obj.r(obj.i));
%             else
%                 obj.m_e(obj.i) = obj.m_e(obj.i) + obj.rho_2(obj.i) * obj.vol_2(obj.i) * (obj.h_2(obj.i) - h_g) / (obj.TAU * obj.r(obj.i));
%             end
            obj.m_dot_pt_1(obj.i) = obj.m_c(obj.i) - obj.m_e(obj.i);
            obj.m_dot_pt_2(obj.i) = obj.m_e(obj.i) - obj.m_c(obj.i);
            obj.m_1(obj.i + 1) = obj.m_1(obj.i) + (obj.m_dot_1b(obj.i) + obj.m_dot_pt_1(obj.i)) * obj.time_step;
            obj.m_2(obj.i + 1) = obj.m_2(obj.i) + (obj.m_dot_2b(obj.i) + obj.m_dot_pt_2(obj.i)) * obj.time_step;
            obj.m_total(obj.i + 1) = obj.m_1(obj.i + 1) + obj.m_2(obj.i + 1);
            obj.dv1dh(obj.i) = IAPWS_IF97('dvdh_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_1(obj.i));
            obj.dv2dh(obj.i) = IAPWS_IF97('dvdh_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_2(obj.i));
            obj.dv1dp(obj.i) = IAPWS_IF97('dvdp_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_1(obj.i));
            obj.dv2dp(obj.i) = IAPWS_IF97('dvdp_ph', obj.p(obj.i) * obj.MPA_PER_BAR, obj.h_2(obj.i));
            obj.term1(obj.i) = (obj.h_1(obj.i) * obj.dv1dh(obj.i) - obj.v_1(obj.i)) * (obj.m_1(obj.i + 1) - obj.m_1(obj.i)) / obj.time_step;
            obj.term2(obj.i) = (obj.h_2(obj.i) * obj.dv2dh(obj.i) - obj.v_2(obj.i)) * (obj.m_2(obj.i + 1) - obj.m_2(obj.i)) / obj.time_step;
            obj.term3(obj.i) = obj.dv1dh(obj.i) * (obj.mh_dot_1b(obj.i) + obj.m_dot_pt_1(obj.i) * h_g + obj.q_21(obj.i) - obj.q_loss_1(obj.i));
            obj.term4(obj.i) = obj.dv2dh(obj.i) * (obj.mh_dot_2b(obj.i) + obj.m_dot_pt_2(obj.i) * h_g - obj.q_21(obj.i) - obj.q_loss_2(obj.i));
            obj.term5(obj.i) = (obj.dv1dp(obj.i) + obj.v_1(obj.i) * obj.dv1dh(obj.i) * 1000) * obj.m_1(obj.i);
            obj.term6(obj.i) = (obj.dv2dp(obj.i) + obj.v_2(obj.i) * obj.dv2dh(obj.i) * 1000) * obj.m_2(obj.i);
            obj.dpdt(obj.i) = ((obj.term1(obj.i) + obj.term2(obj.i) - obj.term3(obj.i) - obj.term4(obj.i)) / (obj.term5(obj.i) + obj.term6(obj.i))) * obj.BAR_PER_MPA;
            obj.p(obj.i + 1) = obj.p(obj.i) + obj.dpdt(obj.i) * obj.time_step;
            if(obj.verbosity > 1 && mod(obj.i, 100) == 0)
                fprintf('%10s = %10.1f\n', 'Time', obj.time(obj.i) + obj.time_step);
                fprintf('%10s = %10f\n', 'p', obj.p(obj.i + 1));
            end
            obj.dh_1dt(obj.i) = (obj.mh_dot_1b(obj.i) + obj.m_dot_pt_1(obj.i) * h_g + obj.q_21(obj.i) - obj.q_loss_1(obj.i) + obj.m_1(obj.i) * obj.v_1(obj.i) * obj.dpdt(obj.i) * 100 - obj.h_1(obj.i) * (obj.m_1(obj.i + 1) - obj.m_1(obj.i)) / obj.time_step) / obj.m_1(obj.i);
            obj.h_1(obj.i + 1) = obj.h_1(obj.i) + obj.dh_1dt(obj.i) * obj.time_step;
            obj.x_1(obj.i + 1) = IAPWS_IF97('x_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_1(obj.i + 1));
            obj.dh_2dt(obj.i) = (obj.mh_dot_2b(obj.i) + obj.m_dot_pt_2(obj.i) * h_g - obj.q_21(obj.i) - obj.q_loss_2(obj.i) + obj.m_2(obj.i) * obj.v_2(obj.i) * obj.dpdt(obj.i) * 100 - obj.h_2(obj.i) * (obj.m_2(obj.i + 1) - obj.m_2(obj.i)) / obj.time_step) / obj.m_2(obj.i);
            obj.h_2(obj.i + 1) = obj.h_2(obj.i) + obj.dh_2dt(obj.i) * obj.time_step;     
            obj.x_2(obj.i + 1) = IAPWS_IF97('x_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_2(obj.i + 1));
            obj.t_1(obj.i + 1) = IAPWS_IF97('T_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_1(obj.i + 1));
            obj.t_2(obj.i + 1) = IAPWS_IF97('T_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_2(obj.i + 1));
            obj.v_1(obj.i + 1) = IAPWS_IF97('v_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_1(obj.i + 1));
            obj.v_2(obj.i + 1) = IAPWS_IF97('v_ph', obj.p(obj.i + 1) * obj.MPA_PER_BAR, obj.h_2(obj.i + 1));
            obj.rho_1(obj.i + 1) = 1 / obj.v_1(obj.i + 1);
            obj.rho_2(obj.i + 1) = 1 / obj.v_2(obj.i + 1);
            obj.vol_1(obj.i + 1) = obj.m_1(obj.i + 1) * obj.v_1(obj.i + 1);
            obj.vol_2(obj.i + 1) = obj.m_2(obj.i + 1) * obj.v_2(obj.i + 1);
            obj.vol_total(obj.i + 1) = obj.vol_1(obj.i + 1) + obj.vol_2(obj.i + 1);
            obj.vol_defect(obj.i + 1) = abs(obj.vol_total(obj.i + 1) - obj.tank_volume) / obj.tank_volume;
            if(obj.vol_defect(obj.i) > obj.VOLUME_TOLERANCE)
                obj.volume_defect_detected = true;
            end
            if(obj.verbosity > 0 && mod(obj.i, 100) == 0 && obj.vol_defect(obj.i + 1) > obj.VOLUME_TOLERANCE)
                disp('Volume defect detected. Results unreliable!');
            end
            obj.water_level(obj.i + 1) = obj.get_water_level(obj.i + 1);
            obj.x(obj.i + 1) = obj.get_quality(obj.i + 1);
            if(obj.verbosity > 1 && mod(obj.i, 100) == 0)
                fprintf('%10s = %10f\n', 'x', obj.x(obj.i + 1));
            end
            obj.time(obj.i + 1) = obj.time(obj.i) + obj.time_step;
            obj.loop_time(obj.i) = toc(start_time);
            obj.i = obj.i + 1;            
        end
        function value = heat_transfer_21(obj)
            value = 5e4 * (obj.t_2(obj.i) - obj.t_1(obj.i)) * obj.vol_1(obj.i) * obj.KW_PER_W;
        end
        function value = heat_loss(obj)
           t_acc = IAPWS_IF97('Tsat_p', obj.p(obj.i) * obj.MPA_PER_BAR);
           t_inf = 313.15;
           r_insul_o = obj.PIPE_RADIUS + obj.PIPE_THICKNESS + obj.INSULATION_THICKNESS;
           r_pipe_i = obj.PIPE_RADIUS - obj.PIPE_THICKNESS;
           sa_tank = 2 * pi * obj.PIPE_RADIUS * obj.tank_length;
           q = (t_acc - t_inf) / (r_insul_o * log(obj.PIPE_RADIUS / r_pipe_i) / obj.K_INSULATION + r_insul_o * log(r_insul_o / obj.PIPE_RADIUS) / obj.K_PIPE + 1 / obj.H_AIR);
           value = q * sa_tank * obj.KW_PER_W;
           value = 0.0;
        end   
        function value = get_quality(obj, loop)
            quality_1 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_1(loop));
            quality_2 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_2(loop));
            m_1 = obj.m_1(loop);
            m_2 = obj.m_2(loop);
            if quality_1 > 0.0
                m_1 = obj.m_1(loop) * (1 - quality_1);
                m_2 = obj.m_2(loop) + (obj.m_1(loop) * quality_1);
            end
            if quality_2 < 1.0
                m_1 = m_1 + (obj.m_2(loop) * (1 - quality_2));
                m_2 = m_2 * quality_2;
            end
            value = m_2 / (m_1 + m_2);
%             value = obj.m_2(loop) / (obj.m_1(loop) + obj.m_2(loop));
        end
        function value = get_water_level(obj, loop)
            quality_1 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_1(loop));
            quality_2 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_2(loop));
            m_1 = obj.m_1(loop);
            m_2 = obj.m_2(loop);
            if quality_1 > 0.0
                m_1 = obj.m_1(loop) * (1 - quality_1);
                m_2 = obj.m_2(loop) + (obj.m_1(loop) * quality_1);
            end
            if quality_2 < 1.0
                m_1 = m_1 + (obj.m_2(loop) * (1 - quality_2));
                m_2 = m_2 * quality_2;
            end
            vol_1 = m_1 * obj.v_1(loop);
            vol_2 = m_2 * obj.v_2(loop);
            value = vol_1 / (vol_1 + vol_2);
            value = obj.vol_1(loop) / obj.vol_total(loop);
        end
        function value = adjust_masses_for_new_qualities(obj, loop)
           quality_1 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_1(loop));
           quality_2 = IAPWS_IF97('x_ph', obj.p(loop) * obj.MPA_PER_BAR, obj.h_2(loop));
           m_1 = obj.m_1(loop);
           m_2 = obj.m_2(loop);
           if quality_1 > 0.0
               m_1 = obj.m_1(loop) * (1 - quality_1);
               m_2 = obj.m_2(loop) + (obj.m_1(loop) * quality_1);
           end
           if quality_2 < 1.0
              m_1 = m_1 + (obj.m_2(loop) * (1 - quality_2));
              m_2 = m_2 * quality_2;
           end
           value = m_2 / (m_1 + m_2);
        end
    end
    methods (Static)
        function [eff, m_dot_stm] = run_turbine(p_turb_in, p_turb_exh, turb_power, cond_depression)
            turb_power = turb_power * 1e6;
            MPA_PER_BAR = 0.1;
            h_turb_in = IAPWS_IF97('hV_p', p_turb_in * MPA_PER_BAR);
            s_turb_in = XSteam('sV_p', p_turb_in);
            t_exh = IAPWS_IF97('Tsat_p', p_turb_exh * MPA_PER_BAR);
            h_exh = XSteam('h_ps', p_turb_exh, s_turb_in);
            w_turb = h_turb_in - h_exh;            
            m_dot_stm = turb_power * 1e-3 / w_turb;
            h_pump_out = IAPWS_IF97('h_pT', p_turb_in * MPA_PER_BAR, t_exh - cond_depression);
            h_cond_out = IAPWS_IF97('h_pT', p_turb_exh * MPA_PER_BAR, t_exh - cond_depression);
            q_cond = h_exh - h_cond_out;
            w_pump = h_pump_out - h_cond_out;
            q_in = w_turb - w_pump + q_cond;
            eff = w_turb / q_in;
        end
        function [eff, m_dot_total, m_dot_1, h_1, m_dot_2, h_2] = run_turbine_test(p_turb_in, h_turb_in, p_turb_exh, turb_power, cond_depression)
            turb_power = turb_power * 1e6;
            MPA_PER_BAR = 0.1;
            x_turb_in = IAPWS_IF97('x_ph', p_turb_in * MPA_PER_BAR, h_turb_in);
            hV_turb_in = IAPWS_IF97('hV_p', p_turb_in * MPA_PER_BAR);
            h_2 = hV_turb_in;
            hL_turb_in = IAPWS_IF97('hL_p', p_turb_in * MPA_PER_BAR);
            h_1 = hL_turb_in;
            sV_turb_in = XSteam('sV_p', p_turb_in);
            t_exh = IAPWS_IF97('Tsat_p', p_turb_exh * MPA_PER_BAR);
            h_exh = XSteam('h_ps', p_turb_exh, sV_turb_in);
            w_turb = hV_turb_in - h_exh;            
            m_dot_stm = turb_power * 1e-3 / w_turb;
            m_dot_total = m_dot_stm / x_turb_in;
            m_dot_1 = m_dot_total * (1 - x_turb_in);
            m_dot_2 = m_dot_stm;
            h_pump_out = IAPWS_IF97('h_pT', p_turb_in * MPA_PER_BAR, t_exh - cond_depression);
            h_cond_out = IAPWS_IF97('h_pT', p_turb_exh * MPA_PER_BAR, t_exh - cond_depression);
            q_cond = h_exh - h_cond_out;
            w_pump = h_pump_out - h_cond_out;
            q_in = w_turb - w_pump + q_cond;
            eff = w_turb / q_in;
        end
        function value = get_minimum_pressure(initial_pressure, max_power, min_efficiency, verbosity)
            min_pressure = 0.0;
            efficiency = 0.0;
            steam_flow = 0.0;
            while efficiency < min_efficiency
                if min_pressure < initial_pressure
                    min_pressure = min_pressure + 1.0;
                end
                if (min_pressure >= initial_pressure)
                    assert(true, 'Minimum efficiency (%.2f) is higher than efficiency at initial pressure (%.2f bar).\nTry setting the initial pressure higher or the minimum efficiency lower.', min_efficiency, initial_pressure);
                end
                [efficiency, steam_flow] = steam_accumulator.run_turbine(min_pressure, 0.9, max_power, 5);
            end
            if verbosity > 0
                fprintf('Minimum pressure allowable for desired minimum efficiency (%.2f) is %.2f bar\n', min_efficiency, min_pressure);
            end
            value = min_pressure;
        end
        function value = get_minimum_steam_mass(initial_pressure, minimum_pressure, power, duration, verbosity)
            i = 0;
            efficiency = 0.0;
            steam_mass = 0.0;
            steam_flow = 0.0;
            current_pressure = initial_pressure;
            while i <= duration
                current_pressure = initial_pressure + (minimum_pressure - initial_pressure) * ( i / duration);
                [efficiency, steam_flow] = steam_accumulator.run_turbine(current_pressure, 0.9, power, 5);
                steam_mass = steam_mass + steam_flow;
                i = i + 1;
            end
            if verbosity > 0
                fprintf('Minimum steam mass is %.6e kg\n', steam_mass); 
            end
            value = steam_mass;
        end
        function value = get_minimum_liquid_mass(steam_mass, initial_pressure, final_pressure, verbosity)
            MPA_PER_BAR = 0.1;
            A = 11.934;
            B = 3985;
            C = 234.1;
            average_pressure = (initial_pressure + final_pressure) / 2;
            c_p_avg = XSteam('CpL_p', average_pressure);
            h_f_ref = IAPWS_IF97('hL_p', average_pressure * MPA_PER_BAR);
            h_g_ref = IAPWS_IF97('hV_p', average_pressure * MPA_PER_BAR);
            t_ref = IAPWS_IF97('Tsat_p', average_pressure * MPA_PER_BAR) - 273.15;
            r_ref = h_g_ref - h_f_ref;
            a = ((B / (A - log(average_pressure))) - C + 273.15) / 647;
            b = (t_ref + 273.15) / 647;
            c = ((1 - a) / (1 - b))^(0.38);
            numerator = steam_mass * r_ref * c;
            d = 1 / (A - log(initial_pressure));
            e = 1 / (A - log(final_pressure));
            denominator = c_p_avg * B * (d - e);
            min_liquid_mass = numerator / denominator;
            if verbosity > 0
                fprintf('Minimum liquid mass is %.6e kg\n', min_liquid_mass); 
            end
            value = min_liquid_mass;
        end
        function value = get_minimum_tank_length(liquid_mass, pressure, quality, verbosity)
            PIPE_RADIUS = 0.4064;
            MPA_PER_BAR = 0.1;
            v_1 = IAPWS_IF97('vL_p', pressure * MPA_PER_BAR);
            v_2 = IAPWS_IF97('vV_p', pressure * MPA_PER_BAR);
            rho_1 = 1 / v_1;
            rho_2 = 1 / v_2;
            rho_mixture = 1 / (quality / rho_2 + (1 - quality) / rho_1);
            steam_mass = liquid_mass * quality / (1 - quality);
            m_total = liquid_mass + steam_mass;
            length = m_total / (rho_mixture * pi * PIPE_RADIUS^2);
            if verbosity > 0
                fprintf('Minimum tank length is %.2f m\n', length);
            end
            value = length;
        end   
        function accumulator = size_accumulator(initial_pressure, initial_quality, minimum_efficiency, max_power, max_duration, time_step, max_iter, verbosity)
            minimum_pressure = steam_accumulator.get_minimum_pressure(initial_pressure, max_power, minimum_efficiency, verbosity);
            required_steam_mass = steam_accumulator.get_minimum_steam_mass(initial_pressure, minimum_pressure, max_power, max_duration, verbosity);
            required_liquid_mass = steam_accumulator.get_minimum_liquid_mass(required_steam_mass, initial_pressure, minimum_pressure, verbosity);
            tank_length = steam_accumulator.get_minimum_tank_length(required_liquid_mass, initial_pressure, initial_quality, verbosity);
            accumulator = steam_accumulator(initial_pressure, initial_quality, tank_length, time_step, max_iter, verbosity);
            accumulator.high_pressure = initial_pressure;
            accumulator.low_pressure = minimum_pressure;
            accumulator.high_quality = initial_quality;
            accumulator.low_quality = 0.01;
            accumulator.minimum_efficiency = minimum_efficiency;
            accumulator.minimum_water_mass = required_liquid_mass;
        end
    end
end

