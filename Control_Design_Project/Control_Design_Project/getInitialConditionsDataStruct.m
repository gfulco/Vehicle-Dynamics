function initial_conditions_data = getInitialConditionsDataStruct( varargin )
%% initial_conditions_data = getInitialConditionsDataStruct( init_LongSpeed , init_LatSpeed )
%   ON INPUT:
%     init_LongSpeed [km/h]: initial longitudinal speed 
%     init_LatSpeed  [km/h]: initial longitudinal speed 
%     init_Omega     [rad/s]: initial yaw-rate
%     init_Phi       [rad]: initial roll angle
%     init_delta     [rad]: initial steering angle
%     init_phi_dot   [rad/2]: initial roll-rate

%   ___      _ _   _      _    ___             _ _ _   _             
%  |_ _|_ _ (_) |_(_)__ _| |  / __|___ _ _  __| (_) |_(_)___ _ _  ___
%   | || ' \| |  _| / _` | | | (__/ _ \ ' \/ _` | |  _| / _ \ ' \(_-<
%  |___|_||_|_|\__|_\__,_|_|  \___\___/_||_\__,_|_|\__|_\___/_||_/__/
%    

if nargin == 0
  initial_conditions_data.longSpeed = 50/3.6;
  initial_conditions_data.latSpeed  = 0;
  initial_conditions_data.Omega = 0;    
  initial_conditions_data.Phi = 0;      
  initial_conditions_data.delta = 0;  
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 1
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = 0;
  initial_conditions_data.Omega = 0;    
  initial_conditions_data.Phi = 0;      
  initial_conditions_data.delta = 0;
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 2
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = varargin{2}/3.6;
  initial_conditions_data.Omega = 0;    
  initial_conditions_data.Phi = 0;      
  initial_conditions_data.delta = 0;
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 3
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = varargin{2}/3.6;
  initial_conditions_data.Omega = varargin{3};    
  initial_conditions_data.Phi = 0;      
  initial_conditions_data.delta = 0;
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 4
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = varargin{2}/3.6;
  initial_conditions_data.Omega = varargin{3};    
  initial_conditions_data.Phi = varargin{4};      
  initial_conditions_data.delta = 0;
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 5
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = varargin{2}/3.6;
  initial_conditions_data.Omega = varargin{3};    
  initial_conditions_data.Phi = varargin{4};      
  initial_conditions_data.delta = varargin{5};
  initial_conditions_data.Phi_dot = 0;   
elseif nargin == 6
  initial_conditions_data.longSpeed = varargin{1}/3.6;
  initial_conditions_data.latSpeed  = varargin{2}/3.6;
  initial_conditions_data.Omega = varargin{3};    
  initial_conditions_data.Phi = varargin{4};      
  initial_conditions_data.delta = varargin{5};
  initial_conditions_data.Phi_dot = varargin{6};   
else
  error('Wrong number of input parameters.')
end

end
                                                                  

