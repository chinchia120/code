clc;
clear all;

%Q19.12
syms EA NA EB NB EC NC EMk1 NMk1 EMk2 NMk2 LAB LBC ThetaA ThetaB ThetaC SLAB SLBC SThetaA SThetaB SThetaC;

F_LAB = sqrt((EB-EA)^2+(NB-NA^2));
F_LBC = sqrt((EC-EB)^2+(NC-NB^2));
F_ThetaA = atan((EB-EA)/(NB-NA))-atan((EMk1-EA)/(NMk1-NA));
F_ThetaB = atan((EC-EB)/(NC-NB))-atan((EA-EB)/(NA-NB));
F_ThetaC = atan((EMk2-EC)/(NMk2-NC))-atan((EB-EC)/(NB-NC));

L = [LAB ; LBC ; ThetaA ; ThetaB ; ThetaC];
A = [diff(F_LAB,EB,1) diff(F_LAB,NB,1) ; 
    diff(F_LBC,EB,1) diff(F_LBC,NB,1) ; 
    diff(F_ThetaA,EB,1) diff(F_ThetaA,NB,1) ; 
    diff(F_ThetaB,EB,1) diff(F_ThetaB,NB,1) ; 
    diff(F_ThetaC,EB,1) diff(F_ThetaC,NB,1)];
P = diag([1/SLAB^2 1/SLBC^2 1/SThetaA^2 1/SThetaB^2 1/SThetaC^2]);

X_0 = simple(subs(inv(A.'*P*A)*A.'*P*L,[EA NA EC NC EMk1 NMk1 EMk2 NMk2 LAB LBC ThetaA ThetaB ThetaC SLAB SLBC SThetaA SThetaB SThetaC] , [5572.32 , 6208.30 , 9552.58 , 6349.45 , 6212.18 , 4956.83 , 11547.42 , 6518.82 , 2717.95 , 2589.28 , dms2rad([254 53 08]) , dms2rad([262 46 20]) , dms2rad([134 34 14]) , 0.024 , 0.024 , 5.3/206265 , 5.1/206265 , 5.2/206265]))





