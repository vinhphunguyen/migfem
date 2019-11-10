% run this after running igaThickCylinder.m
% to check FEM and IGA results. 
% FEM use trilinear B8 elements.

element = elementV;
noNodes = size(node,1);
noElems = size(elementV,1);
noDofs     = noNodes* 3;

% Boundary nodes
% find boundary nodes for bounjdary conditions


xConsNodes  = find(node(:,1)==0)';
yConsNodes2 = find(node(:,2)==-R)';
yConsNodes1 = find(node(:,2)==R)';
yConsNodes  = [yConsNodes2 yConsNodes1];

% essential boundary conditions

u0 = -1;

uFixed     = zeros(size(xConsNodes));
vFixed     = [zeros(size(yConsNodes2)) ...
              u0*ones(size(yConsNodes1))];


udofs = xConsNodes;             % global indecies  of the fixed x disps
vdofs = yConsNodes+numnode;    % global indecies  of the fixed y disps


% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
u = zeros(noDofs,1);        % displacement vector
f = zeros(noDofs,1);        % external force vector

%jacob = zeros(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,Q]=quadrature(  3, 'GAUSS', 3 ); 

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over elements (knot spans)

for e=1:noElems
    sctr   = element(e,:);  %  element scatter vector
    sctrB  = [sctr sctr+noNodes sctr+2*noNodes]; % scatters a B matrix
    nn     = length(sctr);
    
    B      = zeros(6,3*nn);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
    
        [N, dNdxi] = lagrange_basis('B8',pt); 
        J0         = node(sctr,:)'*dNdxi;
        invJ0      = inv(J0);
        dNdx       = dNdxi*invJ0;
        detJ0      = det(J0);
        
        % B matrix
        
        B(1,1:nn)         = dNdx(:,1)';
        B(2,nn+1:2*nn)    = dNdx(:,2)';
        B(3,2*nn+1:3*nn)  = dNdx(:,3)';
        
        B(4,1:nn)         = dNdx(:,2)';
        B(4,nn+1:2*nn)    = dNdx(:,1)';
        
        B(5,2*nn+1:3*nn)  = dNdx(:,2)';
        B(5,nn+1:2*nn)    = dNdx(:,3)';
        
        B(6,1:nn)         = dNdx(:,3)';
        B(6,2*nn+1:3*nn)  = dNdx(:,1)';
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        
        K(sctrB,sctrB) = K(sctrB,sctrB) + B' * C * B * detJ0 * wt;
    end
end

% Computing external force

disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])


bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;

K(:,udofs)=0;
K(:,vdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));


% SOLVE SYSTEM

disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ux    = U(1:noNodes);
Uy    = U(noNodes+1:2*noNodes);
Uz    = U(2*noNodes+1:noDofs);


plot_mesh(node,element,'B8','g.-');
%plot_field(node,elementV,'B8',Uz); do not support B8 elements!!!
view(3)

sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess3d(node,element,'B8','../results/thickCylinderFEM',...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Uz]);



