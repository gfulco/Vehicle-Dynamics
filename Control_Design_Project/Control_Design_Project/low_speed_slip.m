function [kappa_eff , alpha_eff] = low_speed_slip(Vx,Vsx,Vsy,kappa,alpha,Fz,gamma,Vlow,pacejkaParam)


% Coefficients

pCy1 = pacejkaParam.pCy1;
pDy1 = pacejkaParam.pDy1;
pDy2 = pacejkaParam.pDy2;
pDy3 = pacejkaParam.pDy3;

pKy1 = pacejkaParam.pKy1;
pKy2 = pacejkaParam.pKy2;
pKy3 = pacejkaParam.pKy3;

pCx1 = pacejkaParam.pCx1;
pDx1 = pacejkaParam.pDx1;
pDx2 = pacejkaParam.pDx2;
pDx3 = pacejkaParam.pDx3;

pKx1 = pacejkaParam.pKx1;
pKx2 = pacejkaParam.pKx2;
pKx3 = pacejkaParam.pKx3;

Fz0 = pacejkaParam.Fz0;


% Pacejka equations

Fz01 = Fz0;
dfz = Fz / Fz01 - 1;

Cx = pCx1;
mu__x = (dfz * pDx2 + pDx1) * (-pDx3 * gamma ^ 2 + 1);
Dx = mu__x * Fz;
Kxk = Fz * (dfz * pKx2 + pKx1) * exp(-(pKx3 * dfz));
Bx = Kxk / Cx / Dx;
gamma__s = gamma;
Cy = pCy1;
mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1);
Dy = mu__y * Fz;
Kya = Fz01 * pKy1 * sin(0.2e1 * atan((Fz / Fz01 / pKy2))) * (1 - pKy3 * abs(gamma__s));
By = Kya / Cy / Dy;






%   _                  _ _           _ _           _   ___ _ _      
%  | |   ___ _ _  __ _(_) |_ _  _ __| (_)_ _  __ _| | / __| (_)_ __ 
%  | |__/ _ \ ' \/ _` | |  _| || / _` | | ' \/ _` | | \__ \ | | '_ \
%  |____\___/_||_\__, |_|\__|\_,_\__,_|_|_||_\__,_|_| |___/_|_| .__/
%                |___/                                        |_|   

kVlow0 = 770;
CFk = Bx*Cx*Dx;

if abs(Vx) <= Vlow
  kVlow = 1/2*kVlow0*(1+cos(pi*abs(Vx)/Vlow));
else
  kVlow = 0;
end


kappa_eff = kappa + kVlow*(Vsx)/CFk;


%   _         _                _   ___ _ _      
%  | |   __ _| |_ ___ _ _ __ _| | / __| (_)_ __ 
%  | |__/ _` |  _/ -_) '_/ _` | | \__ \ | | '_ \
%  |____\__,_|\__\___|_| \__,_|_| |___/_|_| .__/
%                                         |_|   

CFalpha = By*Cy*Dy; 

if abs(Vx) <= Vlow
  kVlow = 1/2*kVlow0*(1+cos(pi*abs(Vx)/Vlow));
else
  kVlow = 0;
end


alpha_eff = alpha + kVlow*(Vsy)/CFalpha;







