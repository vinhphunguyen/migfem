%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This M file demonstrates h-refinement in isogeometric
% analysis using knot insertion.
%
% V.P. Nguyen
% Johns Hopkins University, 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../C_files/
addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input: control points, knots, order and weights
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

knots      = [0 0 0 0.5 1 1 1];
controlPts = [0 1; 0.5 3; 1.5 1; 2.5 3];
p          = 2;
weights    = [1 1 1 1]'; 
n          = length(knots) - p - 1;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% global h-refinement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% new knots

uniqueKnots = unique(knots);
newKnotsX   = [0.125 0.25 0.375 0.625 0.75 0.875];
nonewkX     = size(newKnotsX,2);
weightedPts = [controlPts(:,1).*weights ...
               controlPts(:,2).*weights weights];
    
[newKnots,newControlPts] = ...
            RefineKnotVectCurve(n-1,p,knots,weightedPts,newKnotsX,nonewkX-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotting ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uniqueKnotsNew = unique(newKnots);

noPts      = 80;
xiArr      = linspace(0,max(knots),noPts);
sCurve     = zeros(2,noPts);
nCurve     = zeros(2,noPts);
iKnot      = zeros(2,length(uniqueKnots));
iKnotNew   = zeros(2,length(uniqueKnotsNew));


for i=1:noPts
  xi          = xiArr(i);   
  sCurve(1,i) = NURBSinterpolation(xi, p, knots, controlPts(:,1), weights);
  sCurve(2,i) = NURBSinterpolation(xi, p, knots, controlPts(:,2), weights);
  
  nCurve(1,i) = NURBSinterpolation(xi, p, newKnots, ...
                             newControlPts(:,1), newControlPts(:,3));
  nCurve(2,i) = NURBSinterpolation(xi, p, newKnots, ...
                             newControlPts(:,2), newControlPts(:,3));
end

for i=1:length(uniqueKnots)
  iKnot(1,i) = NURBSinterpolation(uniqueKnots(i), p, knots, ...
                       controlPts(:,1), weights);
  iKnot(2,i) = NURBSinterpolation(uniqueKnots(i), p, knots, ...
                       controlPts(:,2), weights);
end

for i=1:length(uniqueKnotsNew)
  iKnotNew(1,i) = NURBSinterpolation(uniqueKnotsNew(i), p, ...
                       newKnots, newControlPts(:,1), newControlPts(:,3));
  iKnotNew(2,i) = NURBSinterpolation(uniqueKnotsNew(i), p, ...
                       newKnots, newControlPts(:,2), newControlPts(:,3));
end

hold on
plot(controlPts(:,1),controlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',8);
plot(iKnot(1,:),iKnot(2,:),'rs',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',11);
plot(sCurve(1,:),sCurve(2,:),'b-','LineWidth',1.8);
%axis([0 1 0 1])
axis tight

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
exportfig(gcf,'refine1.eps',opts)

figure

hold on
plot(newControlPts(:,1),newControlPts(:,2),'r-o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',8);
% plot(iKnotNew(1,:),iKnotNew(2,:),'rs',...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',11);
plot(nCurve(1,:),nCurve(2,:),'b-','LineWidth',1.8);
%plot(nCurve(1,:),nCurve(2,:),'b-','LineWidth',1.8);
%axis([0 1 0 1])
axis tight
exportfig(gcf,'refine2.eps',opts)




