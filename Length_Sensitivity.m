clc;
clear;

POWER = 500e6;
INITIAL_PRESSURE = 70.0;
INITIAL_QUALITY = 0.035;
MINIMUM_PRESSURE = 40.0;
SECONDS_PER_HOUR = 3600.0;
HOURS_TO_DISCHARGE = 6.0;
HOURS_TO_WAIT = 6.0;
MINIMUM_EFFICIENCY = 0.25;
TIME_STEP = 1.0;
MAX_ITER = 200000;
VERBOSITY = 1;

tank_length_multiplier = linspace(1.0, 2.0, 3);

accumulators = steam_accumulator.empty(length(tank_length_multiplier), 0);
tank_lengths = zeros(length(tank_length_multiplier));

initial_length = steam_accumulator.size_accumulator(INITIAL_PRESSURE, INITIAL_QUALITY, MINIMUM_EFFICIENCY, POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, TIME_STEP, MAX_ITER, VERBOSITY).tank_length;
parfor i = 1:length(tank_length_multiplier)
    new_length = initial_length * tank_length_multiplier(i);
    tank_lengths(i) = new_length;
    acc = steam_accumulator(INITIAL_PRESSURE, INITIAL_QUALITY, new_length, TIME_STEP, MAX_ITER, VERBOSITY);
    acc.discharge(POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, MINIMUM_PRESSURE, TIME_STEP);
    acc.charge(INITIAL_PRESSURE, INITIAL_QUALITY, TIME_STEP);
    acc.discharge(POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, MINIMUM_PRESSURE, TIME_STEP);
    acc.charge(INITIAL_PRESSURE, INITIAL_QUALITY, TIME_STEP);
    acc.discharge(POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, MINIMUM_PRESSURE, TIME_STEP);
    acc.charge(INITIAL_PRESSURE, INITIAL_QUALITY, TIME_STEP);
    acc.discharge(POWER, HOURS_TO_DISCHARGE * SECONDS_PER_HOUR, MINIMUM_PRESSURE, TIME_STEP);
    acc.charge(INITIAL_PRESSURE, INITIAL_QUALITY, TIME_STEP);
    accumulators(i) = acc;
end
i_all = [accumulators.i];
max_i = max(i_all);
quality_results = zeros(length(tank_length_multiplier), max_i);

acc_with_max_i = findobj(accumulators, 'i', max_i);
time_plot = acc_with_max_i.time(1:max_i);

figure(1)
for k = 1:length(tank_length_multiplier)
    quality_plot = accumulators(k).x(1:max_i);
    quality_results(k, :) = quality_plot;
    plot(time_plot, quality_plot);
    if(k == 1)
        hold on;
    end
end
xlabel('Time [s]');
ylabel('Quality [%]');
xlim([1, max(time_plot)]);

z = quality_results;
surf(time_plot, tank_lengths, z);