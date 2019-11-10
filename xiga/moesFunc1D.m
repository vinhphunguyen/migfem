%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for one dimensional elasticity problems.
%
% One dimensional bar
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/code/xfem-efg-matlab/fem_util');
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-util/

clc
clear all

global p

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PRE-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D bar [0,1]

refinement = 3;

% knotVec    = [0 0 0 1 1 1];
% controlPts = [0 0;
%               0.5 0;
% 	          1 0];
% p       = 2;

knotVec    = [0 0 0 0 1 1 1 1];
controlPts = [0 0;
              1/3 0;
              2/3 0;
	          1 0];
p       = 3;

noGPs   = 3;

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));

%%
generateIGA1DMesh


noCtrPts = size(controlPts,1);% no of control points
noElems  = size(elConn,1);    % no of elements

%% material interface

x0 = 0.45;

ls = controlPts(:,1) - x0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature( 26, 'GAUSS', 1 ); % 2 point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

gps = [];
phi = [];

% Loop over elements (knot spans)

for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);
   lset  = ls(conn);
   alset = abs(lset);
   
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));% coord in parameter space  
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
      
      N       = [];
      dNdxi   = [];
        
      % compute derivative of basis functions w.r.t parameter coord
      
      for in=1:noFns
       [Ni,dNi]  = NURBSbasis (conn(in),p,Xi,knotVec,controlPts(:,3)');
       N         = [N Ni];
       dNdxi     = [dNdxi dNi];
      end
      
      func = dot(N,alset)-abs(dot(N,lset));
      
      % compute the external force, kind of body force
      X       = N * controlPts(conn,1:2);
      
      gps = [gps;X(1)];      
      phi = [phi;func];
    end
end

uu=unique(knotVec);


n3  = plot(gps,phi,'ko-');
hold on
n0=plot(controlPts(:,1),controlPts(:,2)+0.04,'b*');
n1=plot(uu,0.04*ones(length(uu)),'rs-');
set(n1,'MarkerSize',16,'LineWidth',1.01);
set(n0,'MarkerSize',12,'LineWidth',1.01);
set(n3,'MarkerSize',4,'LineWidth',1.001);

% exact solution u = -x^3/6 +1/6x for BCs: u[0]=u[1]=0
% exact solution u = -x^3/6 +7/6x for BCs: u[0]=0;u[1]=1

% x      = linspace(0,1,120);
% uExact1 = -x.^3/6 + 1/6*x;
% uExact2 = -x.^3/6 + 7/6*x;
% %plot(x,uExact1,'r-')
% plot(x,uExact2,'r-')
% %legend('IGA','exact')
% set(gca,'Color',[0.8 0.8 0.8]);

% open the file with write permission

% data1 = [x' uExact1'];
% data2 = [pts dis];

%csvwrite('bar1d1.csv',data1);
%csvwrite('bar1d2.csv',data2);








