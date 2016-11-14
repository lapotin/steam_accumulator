classdef load_forecast
    %LOAD_FORECAST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        seed
        monthly_average_loads = [];
    end
    properties (Constant)
        monthly_peak_loads = [57691.0, 56550.0, 45965.0, 49607.0, 63158.0, 63325.0, 66802.0, 70588.0, 59636.0, 54437.0, 49799.0, 51092.0];
    end
    methods
        function obj = load_forecast(seed_value)
            if nargin > 0
                seed = seed_value;
                s = rng(seed_value);
                obj.generate_monthly_average_loads();
            end
        end
        function value = generate_monthly_average_loads(obj)
            random = rand(1, 12);
            disp(random);
            value = obj.monthly_peak_loads * random;
            disp(value);
        end
    end
    
    methods (Static)
    end
end

