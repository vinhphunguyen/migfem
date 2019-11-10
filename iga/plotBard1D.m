exact     = csvread('bar1d-exact.csv');
p1Elems16 = csvread('bar1d-p1-16elems-C00.csv');
p2Elems16 = csvread('bar1d-p2-16elems-C00.csv');
p3Elems16 = csvread('bar1d-p3-16elems-C00.csv');
%p4Elems16 = csvread('bar1d-p4-32elems-C0.csv');


figure
hold on
plot(exact(:,1),exact(:,2),'b-','LineWidth',1.9);
plot(p1Elems16(:,1),p1Elems16(:,2),'r-','LineWidth',1.8);
plot(p2Elems16(:,1),p2Elems16(:,2),'c-','LineWidth',1.8);
plot(p3Elems16(:,1),p3Elems16(:,2),'k-','LineWidth',1.8);
%plot(p4Elems16(:,1),p4Elems16(:,2),'g-','LineWidth',1.8);
%plot(p5Elems16(:,1),p5Elems16(:,2),'ma-','LineWidth',1.8);
%set(gca,'Color',[0.8 0.8 0.8]);
legend('exact solution','p=1','p=2','p=3')
xlabel('x')
ylabel('u')