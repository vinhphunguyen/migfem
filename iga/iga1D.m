%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for one dimensional elasticity problems.
%
% One dimensional bar example as described in our paper 
% "An introduction to isogeometric analysis"
% VP Nguyen, RN Simpson, SPA Bordas, T Rabczuk
%
% Vinh Phu Nguyen,
% Johns Hopkins University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
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

% 1D bar [0,1]

refinement = 4;

barC0Data
%barC1Data;  % quadratic elements
%barC2Data;  % cubic elements

generateIGA1DMesh


noCtrPts = size(controlPts,1);% no of control points
noElems  = size(elConn,1);    % no of elements

% initialization
K = zeros(noCtrPts,noCtrPts); % global stiffness matrix 
u = zeros(noCtrPts,1);        % displacement vector
f = zeros(noCtrPts,1);        % external force vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  noGPs, 'GAUSS', 1 ); % 2 point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)
for e=1:noElems
   xiE   = elRange(e,:); % [xi_i,xi_i+1]
   conn  = elConn(e,:);
   noFns = length(conn);
 
   % loop over Gauss points 
    for gp=1:size(W,1)                        
      pt      = Q(gp,:);                          
      wt      = W(gp);                            
      Xi      = 0.5 * ( ( xiE(2) - xiE(1) ) * pt + xiE(2) + xiE(1));% coord in parameter space  
      J2      = 0.5 * ( xiE(2) - xiE(1) ); 
            
      [N dNdxi] = NURBS1DBasisDers(Xi,p,knotVec,controlPts(:,3)');

      % compute the jacobian of physical and parameter domain mapping
      % then the derivative w.r.t spatial physical coordinates
      jacob1 = dNdxi*controlPts(conn,1:2);
      J1     = norm (jacob1);
      dNdx   = (1/J1)*dNdxi;
      
      % compute elementary stiffness matrix and
      % assemble it to the global matrix
      K(conn,conn) = K(conn,conn) + dNdx' * dNdx * J1 * J2 * wt;
      
      % compute the external force, kind of body force
      X       = N * controlPts(conn,1:2);
      bx      = X(1);
      f(conn) = f(conn) + bx * N' * J1 * J2 * wt;
    end
end

% Solve the equation
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
udofs  = [1 noCtrPts];  % global indecies  of the fixed x displacements
uFixed = [0 0]'; % BCs: u[0]=u[1]=0
%uFixed = [0 1]'; % BCs: u[0]=0;u[1]=1;

f=f-K(:,udofs)*uFixed;  % modify the  force vector

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(:,udofs)=0;
K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
f(udofs)=bcwt*uFixed;

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noXi = 20;
xi   = linspace(0,1,noXi);
pts  = zeros(noXi,1);
dis  = zeros(noXi,1);

for i=1:noXi
  dis(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     U', controlPts(:,3)');
  pts(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     controlPts(:,1)',controlPts(:,3)');
end

plot(pts,dis,'k-s','LineWidth',1.01);
hold on

% exact solution u = -x^3/6 +1/6x for BCs: u[0]=u[1]=0
% exact solution u = -x^3/6 +7/6x for BCs: u[0]=0;u[1]=1

x      = linspace(0,1,120);
uExact1 = -x.^3/6 + 1/6*x;
uExact2 = -x.^3/6 + 7/6*x;
plot(x,uExact1,'r-','LineWidth',1.1)
%plot(x,uExact2,'r-')
%legend('IGA','exact')
%set(gca,'Color',[0.8 0.8 0.8]);
legend('numerical solution','exact solution')
xlabel('x')
ylabel('u')

% open the file with write permission

data1 = [x' uExact1'];
data2 = [pts dis];

csvwrite('bar1d1.csv',data1);
csvwrite('bar1d2.csv',data2);








