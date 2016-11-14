function out_values = net_revenue_optimization(power, time)
    var0 = [power, time];
    % linear inequality constraints matrix
    A = [];
    % linear inequality contraints vector
    b = [];
    % linear equality contraints matrix
    A_eq = [];
    % linear equality contraints vector
    b_eq = [];
    % upper and lower bounds (power, time)
    lb = [50.0 - eps, 0.5 - eps];
    ub = [1000.0 + eps, 6.0 + eps];
    % non-linear constraints function handle
    nonlcon = [];
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter', 'FinDiffRelStep', 0.01, 'TolX', 1e-6,...
        'TolFun', 0.001, 'TolCon', 1e-6, 'UseParallel', true);
    N = 1;
    [var, net_revenue] = fmincon(@Opti, var0, A, b, A_eq, b_eq, lb, ub, nonlcon, options);
    function net_revenue = Opti(var)
        SECONDS_PER_HOUR = 3600.0;
        HOURS_PER_DAY = 24.0;
        DAYS_PER_YEAR = 365.0;
        charge_pressure = 83.0;
        
        initial_pressure = 70.0;
        initial_quality = 0.06;
        minimum_efficiency = 0.28;
        required_power = var(1,1);
        required_time = var(1,2) * SECONDS_PER_HOUR;
        charge_flow = 1200.0;
        
        fprintf('Determining net revenue for %.2f MW and %0.2f hours.\n', required_power, var(1,2));
        
        acc = steam_accumulator.size_accumulator(initial_pressure, initial_quality, minimum_efficiency, required_power, required_time,...
            10.0, 1000000, 0);
        [avg_m_c, avg_charge_time, avg_m_d, avg_min_quality, avg_max_quality, valid_model] = acc.evaluate_accumulator(required_power,...
            required_time, charge_pressure, charge_flow);
        
        life = 40.0;
        interest = 0.07;
        period = 6.0;
        peak_amplitude = 25.0;
        avg_elec_price = 34.0;
        storage_cycles_per_year = (1 / period) * HOURS_PER_DAY * DAYS_PER_YEAR;
        n = 0.7;
        CSP_MW = 50.0;
        CSP_MWh = 375.0;
        CSP_power_cost = 31.4;
        total_power_cost = CSP_power_cost * (required_power / CSP_MW)^n;
        cost_per_MW = total_power_cost / required_power;
        scaling_CSP_energy_cost = 11.17;
        
        pipe_cost = 218.0;
        insulation_cost = 200.0;
        insulation_thickness = 0.203;
        
        cost_per_sqft = 50.0;
        L_b = 200.0;
        H_b = 10.0;
        W_b = acc.tank_length / (L_b * H_b);
        SA_b = 4 * (L_b * H_b) * 2 * W_b * H_b;
        building_cost = L_b * W_b * 10.764 * cost_per_sqft / 1e6;
        
        unscaled_energy_cost = (pipe_cost * acc.tank_length + insulation_cost * insulation_thickness * SA_b) / 1e6;
        
        scaled_energy_cost = scaling_CSP_energy_cost * ((required_power * required_time) / CSP_MWh)^n;
        
        total_energy_cost = scaled_energy_cost + unscaled_energy_cost + building_cost;
        
        cost_per_MWh = total_energy_cost / (required_power * required_time);
        
        total_capital_cost = total_power_cost + total_energy_cost;
        
        power_om = total_power_cost * 0.05;
        energy_om = total_energy_cost * 0.05;
        total_om = power_om + energy_om;
        
        c1 = (3/4) * period - (avg_charge_time / 2);
        c2 = (3/4) * period + (avg_charge_time / 2);
        d1 = (period / 4) - (required_time / 2);
        d2 = (period / 4) + (required_time / 2);
        
        y = @(t)peak_amplitude * sin((2 * pi() * t) / period) + avg_elec_price;
        integral_charge = integral(y, c1, c2);
        integral_discharge = integral(y, d1, d2);
        average_discharge_price = integral_discharge / (d2 - d1);
        average_charge_price = integral_charge / (c2 - c1);
        delta_price = average_discharge_price - average_charge_price;
        
        charge_revenue = average_charge_price * avg_charge_time * storage_cycles_per_year * required_power / 1e6;
        discharge_revenue = average_discharge_price * required_time * storage_cycles_per_year * required_power / 1e6;
        capital_cost = total_capital_cost * (interest + (interest / ((1 + interest)^life - 1)));
        
        net_revenue = (discharge_revenue - charge_revenue - capital_cost - total_om) * -1; 
        
        fprintf('Net revenue for %.2f MW and %0.2f hours is $ %.3e MM.\n', required_power, var(1,2), (discharge_revenue - charge_revenue - capital_cost - total_om));
        
        out_values(N, 1) = var(1,1);
        out_values(N, 2) = var(1,2);
        out_values(N, 3) = acc.tank_length;
        out_values(N, 4) = avg_charge_time / SECONDS_PER_HOUR;
        out_values(N, 5) = avg_max_quality;
        out_values(N, 6) = net_revenue * -1;
        N = N + 1;
        save out_values;
    end
end