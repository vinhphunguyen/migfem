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

alpha = 50;

% 1D bar [0,1]

refinement = 2;

barCkData


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

noGPs = 10;

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

      if (X(1)-0.42>1e-8) && (X(1)- 0.58<1e-8)
          bx = (2*alpha^2 -4*alpha^4*(X(1)-0.5)^2)*exp(-alpha^2*(X(1)-0.5)^2);
      else
          bx = 0;
      end
      
      f(conn) = f(conn) + bx * N' * J1 * J2 * wt;
    end
end

% Solve the equation
bcwt=mean(diag(K)); % a measure of the average  size of an element in K
                    % used to keep the  conditioning of the K matrix
udofs  = [1 noCtrPts];  % global indecies  of the fixed x displacements
uFixed = [0 1]'; % BCs: u[0]=0;u[1]=1;

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

noXi = 1600;
xi   = linspace(0,1,noXi);
pts  = zeros(noXi,1);
dis  = zeros(noXi,1);

for i=1:noXi
  dis(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     U', controlPts(:,3)');
  pts(i,1) = NURBSinterpolation(xi(i), p, knotVec, ...
                     controlPts(:,1)',controlPts(:,3)');
end


% exact solution u = x + exp()

x      = linspace(0,1,1600);
uExact1 = x + exp(-alpha^2*(x-0.5).^2);

plot(pts,dis,'b--','LineWidth',1.01);
hold on
plot(x,uExact1,'r-','LineWidth',1.1)
%legend('IGA','exact')
%set(gca,'Color',[0.8 0.8 0.8]);
legend('numerical solution','exact solution')
xlabel('x')
ylabel('u')

% open the file with write permission

data1 = [x' uExact1'];
data2 = [pts dis];

%csvwrite('bar1d-exact.csv'     ,data1);
csvwrite('bar1d-p3-16elems-C00.csv',data2);








