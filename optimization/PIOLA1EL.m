function stress = PIOLA1EL(lamda)
global DFGRD0 STRETCH_NEW STATEV DTIME ...
       MU  ALPHA  KELAS  ...
       MUVIS ALPHAVIS  KVIS ETADEV ETAVOL

DFGRD1=zeros(3,3);
DFGRD1(1,1) = STRETCH_NEW;
DFGRD1(2,2) = lamda;
DFGRD1(3,3) = lamda;
[~, ~, ~, FIRSTPIOLA, ~] = VISCOEL_OGDEN(DFGRD0, DFGRD1, STATEV, DTIME, ...
                                                             MU, ALPHA, KELAS, ...
                                                             MUVIS, ALPHAVIS, KVIS, ETADEV, ETAVOL);
stress = FIRSTPIOLA(2,2);

 