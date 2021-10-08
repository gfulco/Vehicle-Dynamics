function [mu,dmu] = burckhardt(lambda,road_condition)
%% mu = BURCKHARDT(lambda,road_condition)
%   Burckhardt Friction model.
%   Return friction coefficient mu and its slope dmu for lambda \in [0,1]
%   and road_condition, according to:
%   road_condition = 1; --> Dry Asphalt
%   road_condition = 2; --> Wet Asphalt
%   road_condition = 3; --> Cobblestone
%   road_condition = 4; --> Snow

switch road_condition
    case 1 % Dry Asphalt
        theta = [1.28 23.99 0.52];
    case 2 % Wet Asphalt
        theta = [0.86 33.82 0.35];
    case 3 % Cobblestone
        theta = [1.37 6.46 0.67];
    case 4 % Snow
        theta = [0.19 94.13 0.06];
    otherwise % Dry Asphalt
        theta = [1.28 23.99 0.52];
end
   mu  = theta(1)*(1 - exp(-lambda*theta(2))) - lambda*theta(3);
   dmu = theta(1)*theta(2)*exp(-lambda*theta(2)) - theta(3);
end