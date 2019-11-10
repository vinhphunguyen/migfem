addpath ../nurbs-geopdes/inst/
addpath ../C_files/

% for exporting to EPS file

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);

%% original curve (quadratic)

p = 2;
          
controlPts = [0 0.5 2   3 4;
              0 1.0 0.5 1.4 0];          

weights    = ones(5,1);

uKnot      = [0 0 0 1 2 3 3 3];          
%uKnot      = 1/max(uKnot)*uKnot;

% make the curve 
curve      = nrbmak(controlPts, uKnot);

[sCurve,iKnot] = plot2DNURBSCurve (uKnot, controlPts', p, weights, 80);
[BsplineVals, NURBSderivs,xi] = getAllShapeGrads (uKnot,p,weights);

figure
hold on
plot(sCurve(1,:),sCurve(2,:),'b-','LineWidth',1.8);
plot(controlPts(1,:),controlPts(2,:),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);          
plot(iKnot(1,:),iKnot(2,:),'rs',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',11);            
axis on

figure
plot(xi, BsplineVals,'LineWidth',1.4)

%% refined curve 1

curve1      = nrbkntins(curve,[1/3]);
controlPts1 = curve1.coefs(1:2,:);
uKnot1      = curve1.knots;
weights1    = curve1.coefs(4,:);

[sCurve1,iKnot1] = plot2DNURBSCurve (uKnot1, controlPts1', p, weights1, 80);
[BsplineVals1, NURBSderivs,xi] = getAllShapeGrads (uKnot1,p,weights1);


figure
hold on

plot(sCurve1(1,:),sCurve1(2,:),'b-','LineWidth',1.8);
plot(controlPts1(1,:),controlPts1(2,:),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);           

% plot(iKnot1(1,:),iKnot1(2,:),'rs',...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',11);                
axis on

figure
plot(xi, BsplineVals1,'LineWidth',1.4)

%% refined curve 2

curve2      = nrbkntins(curve1,[2/3]);
controlPts2 = curve2.coefs(1:2,:);
uKnot2      = curve2.knots;
weights2    = curve2.coefs(4,:);

[sCurve2,iKnot2] = plot2DNURBSCurve (uKnot2, controlPts2', p, weights2, 80);
[BsplineVals2, NURBSderivs,xi] = getAllShapeGrads (uKnot2,p,weights2);


figure
hold on

plot(sCurve2(1,:),sCurve2(2,:),'b-','LineWidth',1.8);
plot(controlPts2(1,:),controlPts2(2,:),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);                       
axis on

figure
plot(xi, BsplineVals2,'LineWidth',1.4)