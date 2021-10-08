function [Fx,Fy,Mz] = pacejka_model(kappa,Fz,alpha,gamma,pacejkaParam)
%% [Fx_comb,Fy_comb,Mz] = PACEJKA_MODEL(kappa,Fz,alpha,gamma)
%% Pacejka Model

% Required Inputs (Units and Sign Convention)
% kappa: slip ratio (%, i.e. 0.10 = 10%) where a + slip ratio generates a driving force
% Fz: normal force (N)
% Dfz  = (Fz - Fz0)/Fz0
% Fz0 = Nominal Load as defined in Coefficient Tables in Appendix I
% alpha: slip angle (rad) where a + slip angle generates a positive lateral force
% gamma: camber angle (rad) where a + camber angle indicates the top of the tire tipped right (same sign convention for left and right tires)
% lambda: multiplier, which should be left at 1.0 for most applications.
% 
% Resulting Outputs (Units and Sign Convention)
% Fy : Lateral Force (N)
% Fx : Longitudinal Force (N)
% Mz : Aligning Torque (N-m)

% Coefficients

rBy1 = pacejkaParam.rBy1;
rBy2 = pacejkaParam.rBy2;
rBy3 = pacejkaParam.rBy3;
rBy4 = pacejkaParam.rBy4;
rCy1 = pacejkaParam.rCy1;
rEy1 = pacejkaParam.rEy1;
rEy2 = pacejkaParam.rEy2;
rHy1 = pacejkaParam.rHy1;
rHy2 = pacejkaParam.rHy2;
rVy1 = pacejkaParam.rVy1;
rVy2 = pacejkaParam.rVy2;
rVy3 = pacejkaParam.rVy3;
rVy4 = pacejkaParam.rVy4;
rVy5 = pacejkaParam.rVy5;
rVy6 = pacejkaParam.rVy6;
rBx1 = pacejkaParam.rBx1;
rBx2 = pacejkaParam.rBx2;
rCx1 = pacejkaParam.rCx1;
rEx1 = pacejkaParam.rEx1;
rEx2 = pacejkaParam.rEx2;
rHx1 = pacejkaParam.rHx1;
pCy1 = pacejkaParam.pCy1;
pDy1 = pacejkaParam.pDy1;
pDy2 = pacejkaParam.pDy2;
pDy3 = pacejkaParam.pDy3;
pEy1 = pacejkaParam.pEy1;
pEy2 = pacejkaParam.pEy2;
pEy3 = pacejkaParam.pEy3;
pEy4 = pacejkaParam.pEy4;
pHy1 = pacejkaParam.pHy1;
pHy2 = pacejkaParam.pHy2;
pHy3 = pacejkaParam.pHy3;
pKy1 = pacejkaParam.pKy1;
pKy2 = pacejkaParam.pKy2;
pKy3 = pacejkaParam.pKy3;
pVy1 = pacejkaParam.pVy1;
pVy2 = pacejkaParam.pVy2;
pVy3 = pacejkaParam.pVy3;
pVy4 = pacejkaParam.pVy4;
qBz1 = pacejkaParam.qBz1;
qBz2 = pacejkaParam.qBz2;
qBz3 = pacejkaParam.qBz3;
qBz4 = pacejkaParam.qBz4;
qBz5 = pacejkaParam.qBz5;
qBz9 = pacejkaParam.qBz9;
qBz10 =pacejkaParam.qBz1;
qCz1 = pacejkaParam.qCz1;
qDz1 = pacejkaParam.qDz1;
qDz2 = pacejkaParam.qDz2;
qDz3 = pacejkaParam.qDz3;
qDz4 = pacejkaParam.qDz4;
qDz6 = pacejkaParam.qDz6;
qDz7 = pacejkaParam.qDz7;
qDz8 = pacejkaParam.qDz8;
qDz9 = pacejkaParam.qDz9;
qEz1 = pacejkaParam.qEz1;
qEz2 = pacejkaParam.qEz2;
qEz3 = pacejkaParam.qEz3;
qEz4 = pacejkaParam.qEz4;
qEz5 = pacejkaParam.qEz5;
qHz1 = pacejkaParam.qHz1;
qHz2 = pacejkaParam.qHz2;
qHz3 = pacejkaParam.qHz3;
qHz4 = pacejkaParam.qHz4;
R0   = pacejkaParam.R0;
pCx1 = pacejkaParam.pCx1;
pDx1 = pacejkaParam.pDx1;
pDx2 = pacejkaParam.pDx2;
pDx3 = pacejkaParam.pDx3;
pEx1 = pacejkaParam.pEx1;
pEx2 = pacejkaParam.pEx2;
pEx3 = pacejkaParam.pEx3;
pEx4 = pacejkaParam.pEx4;
pHx1 = pacejkaParam.pHx1;
pHx2 = pacejkaParam.pHx2;
pKx1 = pacejkaParam.pKx1;
pKx2 = pacejkaParam.pKx2;
pKx3 = pacejkaParam.pKx3;
pVx1 = pacejkaParam.pVx1;
pVx2 = pacejkaParam.pVx2;
Fz0  = pacejkaParam.Fz0 ;


%% Pacejka equations

Fz01 = Fz0;
dfz = Fz / Fz01 - 1;
SHx = dfz * pHx2 + pHx1;
SVx = Fz * (dfz * pVx2 + pVx1);

%********
SVx = 0; 
SHx = 0;
%********

kappa__x = kappa + SHx; 



Cx = pCx1*1.3;
mu__x = (dfz * pDx2 + pDx1) * (-pDx3 * gamma ^ 2 + 1);
Dx = mu__x * Fz;
Kxk = Fz * (dfz * pKx2 + pKx1) * exp(-(pKx3 * dfz));
Ex = (dfz ^ 2 * pEx3 + dfz * pEx2 + pEx1) * (1 - pEx4);
Bx = Kxk / Cx / Dx;
gamma__s = gamma;
SHy = dfz * pHy2 + pHy3 * gamma__s + pHy1;
SVy = Fz * (pVy1 + pVy2 * dfz + (dfz * pVy4 + pVy3) * gamma__s);
%********
SVy = 0; 
SHy = 0;
%********
alpha__y = alpha + SHy;
Cy = pCy1;
mu__y = (dfz * pDy2 + pDy1) * (-pDy3 * gamma__s ^ 2 + 1);
Dy = mu__y * Fz;
Ey = (dfz * pEy2 + pEy1) * (-pEy4 * gamma__s - pEy3 + 1);
Kya = Fz01 * pKy1 * sin(0.2e1 * atan((Fz / Fz01 / pKy2))) * (1 - pKy3 * abs(gamma__s));
By = Kya / Cy / Dy;
gamma__z = gamma;
SHf = SHy + SVy / Kya;
SHt = qHz1 + qHz2 * dfz + (dfz * qHz4 + qHz3) * gamma__z;
alpha__t = alpha + SHt;
alpha__r = alpha + SHf;
Bt = (dfz ^ 2 * qBz3 + dfz * qBz2 + qBz1) * (1 + qBz4 * gamma__z + qBz5 * abs(gamma__z));
Ct = qCz1;
Dt = Fz * (dfz * qDz2 + qDz1) * (qDz4 * gamma__z ^ 2 + qDz3 * gamma__z + 1) * R0 / Fz0;
Et = (dfz ^ 2 * qEz3 + dfz * qEz2 + qEz1) * (0.1e1 + (qEz5 * gamma__z + qEz4) * atan((Bt * Ct * alpha__t)));
Br = By * Cy * qBz10 + qBz9;
Dr = Fz * (qDz6 + qDz7 * dfz + (dfz * qDz9 + qDz8) * gamma__z) * R0;
SHxa = rHx1;
Bxa = rBx1 * (kappa ^ 2 * rBx2 ^ 2 + 1) ^ (-0.1e1 / 0.2e1);
Cxa = rCx1;
Dxa = 0.1e1 / cos(Cxa * atan((Bxa * SHxa)));
Gxa = Dxa * cos(Cxa * atan((Bxa * (alpha + SHxa))));
SHyk = rHy1;
DVyk = mu__y * Fz * (dfz * rVy2 + gamma * rVy3 + rVy1) * (alpha ^ 2 * rVy4 ^ 2 + 1) ^ (-0.1e1 / 0.2e1);
SVyk = DVyk * sin(rVy5 * atan((rVy6 * kappa)));
Byk = rBy1 * (1 + rBy2 ^ 2 * (alpha - rBy3) ^ 2) ^ (-0.1e1 / 0.2e1);
Cyk = rCy1;
Dyk = 0.1e1 / cos(Cyk * atan((Byk * SHyk)));
Gyk = Dyk * cos(Cyk * atan((Byk * (kappa + SHyk))));

%% Output

% Fx
Fx0 = (-Dx * sin(Cx * atan(-Bx * kappa__x + Ex * (Bx * kappa__x - atan(Bx * kappa__x)))) + SVx);
% fprintf(1,' kappa_x = %12.6f, Fx = %12.6f, SHx = %12.6f, SVx = %12.6f\n', kappa__x, Fx0,SHx,SVx);
Fx_comb =  Fx0* Gxa;

%  offsetFx = (Fz * (-0.281372732980853925e4 + Fz * 0.810191533417750187e0) * sin(atan(exp(Fz * (-0.616138223684905249e-3) + 0.548143673871933856e0) * (Fz * (-0.101430342064303437e-1) + 0.480851687881085454e5) / (-0.281372732980853925e4 + Fz * 0.810191533417750187e0) * (Fz * 0.101894089591398932e-5 - 0.552431757689196754e-4) * 0.693990264255749834e0 + (Fz ^ 2 * (-0.333165990314313833e-6) + Fz * 0.350069496442759621e-3 + 0.865566227723819082e-1) * (exp(Fz * (-0.616138223684905249e-3) + 0.548143673871933856e0) * (Fz * (-0.101430342064303437e-1) + 0.480851687881085454e5) / (-0.281372732980853925e4 + Fz * 0.810191533417750187e0) * (Fz * 0.101894089591398932e-5 - 0.552431757689196754e-4) * (-0.693990264255749834e0) - atan(exp(Fz * (-0.616138223684905249e-3) + 0.548143673871933856e0) * (Fz * (-0.101430342064303437e-1) + 0.480851687881085454e5) / (-0.281372732980853925e4 + Fz * 0.810191533417750187e0) * (Fz * 0.101894089591398932e-5 - 0.552431757689196754e-4) * (-0.693990264255749834e0)))) * 0.144094240439009846e1) * 0.112404512366744451e-2 + Fz * (Fz * (-0.261373983233011906e-1) - 0.539868736564886476e2) * 0.112404512366744451e-2) / cos(atan(0.293163545351097649e-1) * 0.926013692178949799e0) * cos(atan(-0.293163545351097649e-1) * 0.926013692178949799e0);

Fx = Fx_comb;% - offsetFx;

% Fy
Fy0 = (-Dy * sin(Cy * atan(-By * alpha__y + Ey * (By * alpha__y - atan(By * alpha__y)))) + SVy);
Fy_comb = Fy0  * Gyk + SVyk;
% offsetFy = Fz * (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * sin(atan(sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * (-0.152511103219044525e8) + (0.729166688663670237e-1 + Fz * (-0.365078573864376527e-4)) * (sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * 0.152511103219044525e8 - atan(sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * 0.152511103219044525e8))) * 0.154027309068454099e1) * 0.112404512366744451e-2 + Fz * (0.781490003250828664e2 + Fz * (-0.401680086940341061e-1)) * 0.112404512366744451e-2 + Fz * (0.509629360051086987e6 + Fz * (-0.580349220782825483e2)) * sin(atan(0.0e0) * 0.906119916069003523e-1) * 0.126347744004056056e-5;

Fy = Fy_comb;% - offsetFy;
% Mz0
Mz0 = Dr * ((Br ^ 2 * alpha__r ^ 2 + 1) ^ (-0.1e1 / 0.2e1)) * cos(alpha);
Mz0 = 0;
% fprintf('Mz_0 = %12.6f\n',Mz0)
Mz_comb = -((-Dy * sin(Cy * atan(-By * alpha__y + Ey * (By * alpha__y - atan(By * alpha__y)))) + SVy) * Gyk + SVyk) * Dt * cos(Ct * atan(-Bt * alpha__t + Et * (Bt * alpha__t - atan(Bt * alpha__t)))) * cos(alpha) + Mz0;

% offsetMz0 = (Fz * (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * sin(atan(sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * (-0.152511103219044525e8) + (0.729166688663670237e-1 + Fz * (-0.365078573864376527e-4)) * (sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * 0.152511103219044525e8 - atan(sin(0.2e1 * atan(Fz * 0.110511927043874963e-2)) * (abs(0) * 0.134159743861844283e1 - 1) / Fz / (-0.240954039813793634e4 + Fz * 0.274390567365253968e0) * (Fz * (-0.132636164465592824e-5) - 0.371083630262798898e-2) * 0.152511103219044525e8))) * 0.154027309068454099e1) * 0.112404512366744451e-2 + Fz * (0.781490003250828664e2 + Fz * (-0.401680086940341061e-1)) * 0.112404512366744451e-2) * Fz * (0.224309955325285813e3 + Fz * (-0.549111816154102678e-1)) * cos(atan(-(Fz ^ 2 * abs(0) * 0.219500142448391398e-5 + 0.589923028704022467e1 + abs(0) * Fz * (-0.241301554085827930e-2) + (abs(0) * (-0.166097760625274660e1)) + Fz ^ 2 * (-0.779590214501757490e-5) + Fz * 0.857021449786109918e-2) * (-0.299697690075865792e-2 + Fz * (-0.862586787573066365e-5)) + (-0.117925961287417191e1 + Fz ^ 2 * atan((Fz ^ 2 * abs(0) * 0.173726997801674199e1 + 0.466904283376112022e7 + abs(0) * Fz * (-0.190982083604184959e4) + (abs(0) * (-0.131460804412892833e7)) + Fz ^ 2 * (-0.617019496981862137e1) + Fz * 0.678303721638747720e4) * (-0.266624251789853561e1 + Fz * (-0.767395160043653126e-2)) * 0.237313360298639370e-8) * (-0.705810235338949284e-6) + atan((Fz ^ 2 * abs(0) * 0.173726997801674199e1 + 0.466904283376112022e7 + abs(0) * Fz * (-0.190982083604184959e4) + (abs(0) * (-0.131460804412892833e7)) + Fz ^ 2 * (-0.617019496981862137e1) + Fz * 0.678303721638747720e4) * (-0.266624251789853561e1 + Fz * (-0.767395160043653126e-2)) * 0.237313360298639370e-8) * Fz * 0.941721905257777124e-3 + atan((Fz ^ 2 * abs(0) * 0.173726997801674199e1 + 0.466904283376112022e7 + abs(0) * Fz * (-0.190982083604184959e4) + (abs(0) * (-0.131460804412892833e7)) + Fz ^ 2 * (-0.617019496981862137e1) + Fz * 0.678303721638747720e4) * (-0.266624251789853561e1 + Fz * (-0.767395160043653126e-2)) * 0.237313360298639370e-8) * (-0.194289083926614176e0) + Fz ^ 2 * (-0.428399520995642626e-5) + Fz * 0.571588782542656153e-2) * ((Fz ^ 2 * abs(0) * 0.219500142448391398e-5 + 0.589923028704022467e1 + abs(0) * Fz * (-0.241301554085827930e-2) + (abs(0) * (-0.166097760625274660e1)) + Fz ^ 2 * (-0.779590214501757490e-5) + Fz * 0.857021449786109918e-2) * (-0.299697690075865792e-2 + Fz * (-0.862586787573066365e-5)) - atan((Fz ^ 2 * abs(0) * 0.219500142448391398e-5 + 0.589923028704022467e1 + abs(0) * Fz * (-0.241301554085827930e-2) + (abs(0) * (-0.166097760625274660e1)) + Fz ^ 2 * (-0.779590214501757490e-5) + Fz * 0.857021449786109918e-2) * (-0.299697690075865792e-2 + Fz * (-0.862586787573066365e-5))))) * 0.167097884314218659e1) * cos(0.0e0) * (-0.288830942793272139e-6) + Fz * (-0.125523877321407813e2 + Fz * 0.171847774709159469e-1) * cos(0.0e0) * 0.256956715270377836e-3;

Mz = Mz_comb;% - offsetMz0;

% fprintf('kappa = %6.3f, Fx_comb=%6.3f, Fy_comb=%6.3f, Mz=%6.3f, Fz=%6.3f\n',kappa__x, Fx_comb, Fy_comb, Mz, Fz);







