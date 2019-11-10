AuxStress(1,1) = K1*FACStress1/SQR*CT2*(1-ST2*S3T2);
AuxStress(2,2) = K1*FACStress1/SQR*CT2*(1+ST2*S3T2);
AuxStress(1,2) = K1*FACStress1/SQR*ST2*CT2*C3T2;
AuxStress(2,1) = AuxStress(1,2);

u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);

u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);

AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;

AuxEps(1,1) = AuxGradDisp(1,1);
AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
AuxEps(1,2) = AuxEps(2,1);
AuxEps(2,2) = AuxGradDisp(2,2);