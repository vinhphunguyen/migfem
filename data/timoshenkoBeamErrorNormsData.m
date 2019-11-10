
h=[12/4;12/8;12/16;12/32];

%% p=q=2
energy_norm2=[0.280367922203176;0.065949927784214;0.015513554945101;0.003676401977117];
disp_norm2  =[0.002529672714234;2.941750524748104e-04;3.222393506930819e-05;3.519741260737019e-06];

%% p=q=4

energy_norm3=[0.050260635432111;0.005171726586377;0.001815128452680;5.510207982135889e-04];
disp_norm3  =[4.174278223106901e-04;3.321230106750417e-05;9.724656549152937e-06;
1.774991938805558e-06];

figure
hold on
plot(log10(h),log10(energy_norm2),'s-','LineWidth',1.1)
%plot(log10(h),log10(energy_norm3),'*-','LineWidth',1.1)
energy_slope2 = polyfit(log10(h),log10(energy_norm2),1);
%energy_slope3 = polyfit(log10(h),log10(energy_norm3),1);
xlabel('log_{10}(h)')
ylabel('log_{10}(e_{energy})')
rate2  = energy_slope2(1)
%rate3  = energy_slope3(1)

figure
hold on
plot(log10(h),log10(disp_norm2),'s-','LineWidth',1.1)
%plot(log10(h),log10(disp_norm3),'*-','LineWidth',1.1)
disp_slope2 = polyfit(log10(h),log10(disp_norm2),1);
disp_slope3 = polyfit(log10(h),log10(disp_norm3),1);
rate2  = disp_slope2(1)
%rate3  = disp_slope3(1)
xlabel('log_{10}(h)')
ylabel('log_{10}(e_{displacement})')