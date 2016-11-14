clear all 
close all
clc

load = load_forecast(1);


% POWER = 500.0;
% HOURS_TO_DISCHARGE = 6.0;
% HOURS_TO_WAIT = 6.0;
% 
% CHARGE_PRESSURE = 70.0;
% CHARGE_FLOWRATE = 600.0;
% 
% START_PRESSURE = 70.0;
% INITIAL_PRESSURE = 35.0;
% MINIMUM_PRESSURE = 40.0;
% 
% INITIAL_QUALITY = 0.06;
% MINIMUM_QUALITY = 0.03;
% QUALITY = 0.06;
% 
% MINIMUM_EFFICIENCY = 0.25;
% 
% MAX_ITER = 1000000;
% 
% VERBOSITY = 3;
% 
% TIME_STEP = 10.0;
% 
% SECONDS_PER_HOUR = 3600.0;
% 
% start_time = tic;
% 
% acc = steam_accumulator.size_accumulator(START_PRESSURE, QUALITY,...
%     MINIMUM_EFFICIENCY, POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR,...
%     TIME_STEP, MAX_ITER, VERBOSITY);

% [avg_m_c, avg_charge_time, avg_m_d, avg_min_quality, avg_max_quality, valid_model] = acc.evaluate_accumulator(POWER,...
%             HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, CHARGE_PRESSURE, CHARGE_FLOWRATE);

% for i = 1:10
%     acc.discharge_test(POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, TIME_STEP);
% 
%     acc.charge(START_PRESSURE, INITIAL_QUALITY, CHARGE_PRESSURE, CHARGE_FLOWRATE,...
%     TIME_STEP);
% 
%     acc.wait(HOURS_TO_WAIT * SECONDS_PER_HOUR, TIME_STEP);
% end
% 
% acc.get_plots();

% fprintf('Total time = %.2f seconds.\n', toc(start_time));