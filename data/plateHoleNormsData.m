%loglog(h3,energy_norm3,'ob-','LineWidth',1.05)

h=[1.9039;1.0097;0.5194;0.2633;0.1326];
%h=[0.0625;0.0312;0.0156;0.0078];
hd=[120;360;1224;4488];
energy_norm=[3.4757e-04;1.2063e-04;3.2917e-05;8.2415e-06;2.0553e-06];
disp_norm=[6.1544e-07;7.9872e-08;8.1171e-09;8.6544e-10;1.0177e-10];

h=h(2:end);
energy_norm=energy_norm(2:end);
disp_norm=disp_norm(2:end);

figure
plot(log10(h),log10(energy_norm),'s-','LineWidth',1.1)
slope1 = polyfit(log10(h),log10(energy_norm),1);
xlabel('log_{10}(h)')
ylabel('log_{10}(e_{energy})')
rate1  = slope1(1)

figure
plot(log10(h),log10(disp_norm),'s-','LineWidth',1.1)
slope2 = polyfit(log10(h),log10(disp_norm),1);
rate2  = slope2(1)
xlabel('log_{10}(h)')
ylabel('log_{10}(e_{displacement})')
