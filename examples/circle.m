addpath ../C_files/
addpath ../fem_util/

clear all

knots      = [0 0 0 1 1 2 2 3 3 4 4 4];
controlPts = [1 0; 1 1; 0 1; -1 1;-1 0;-1 -1;0 -1; 1 -1; 1 0];
p          = 2;

weights    = ones(1,9)';
weights(2) = 1/sqrt(2);
weights(4) = 1/sqrt(2);
weights(6) = 1/sqrt(2);
weights(8) = 1/sqrt(2);

%knots      = 1/max(knots)*knots;
noPts      = 100;
xi         = linspace(0,max(knots),noPts);
sCurve     = zeros(2,noPts);



for i=1:noPts
  sCurve(1,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,1), weights);
  sCurve(2,i) = NURBSinterpolation(xi(i), p, knots, controlPts(:,2), weights);
end

hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);
plot(sCurve(1,:),sCurve(2,:),'b-','LineWidth',1.8);
axis([-1 1 -1 1])
axis equal
%axis tight
axis off

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'circle.eps',opts)



