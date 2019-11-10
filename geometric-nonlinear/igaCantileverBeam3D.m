%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
% Geometrical nonlinearity is allowed. However strain is small.
% Total Lagrangian formulation is used.
%
% The nodal displacement are stored as [ux1 uy1 uz1 ux2 uy2 uz3...]
%
% Thin bending strip is bent into a circle. External force is a follower
% load which mut be evaluated in the current configuration.
%
% Each converged load step is written to a VTU file.
% A PVD file that collects all these VTU files is created at the end
% of the loading process.
%
% Vinh Phu Nguyen,
% Cardiff University, UK
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q r element index elRangeU elRangeV elRangeW uKnot vKnot wKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP noPtsX noPtsY noPtsZ C

E0          = 1.2e5;  % Young modulus
nu0         = 0.0;  % Poisson ratio



%% COMPUTE ELASTICITY MATRIX

C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/(1+nu0)*eye(3);


tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input file

cantileverBeam3DData

vtuFileName = 'cantileverBeam3D';

noCtrPts       = noPtsX * noPtsY *noPtsZ;
noDofs         = noCtrPts * 3;
noGPs          = p+1;

% find boundary nodes for boundary conditions

leftNodes    =  find(controlPts(:,1)==0);
fixedXNodes  =  leftNodes;
fixedYNodes  =  leftNodes;
fixedZNodes  =  leftNodes;

% forced nodes
EPS=1e-8;
forcedNodes1  =  find(abs(controlPts(:,1)-a)<EPS);
forcedNodes2  =  find(abs(controlPts(:,3)-t)<EPS);
forcedNodes3  =  find(abs(controlPts(:,3))  <EPS);
forcedNodesP  =  intersect(forcedNodes1,forcedNodes3);
forcedNodesM  =  intersect(forcedNodes1,forcedNodes2);

% forced dofs (x/z-direction)
forcedDofsX    =  3*forcedNodes1-2;
forcedDofsZ    =  3*forcedNodes1;

% recorded point for load-displacement curve plot
recordedPnt    = forcedNodesP(1);
recordedPtXDof = 3*recordedPnt-2; % x-direction dof
recordedPtZDof = 3*recordedPnt;   % z-direction dof

I      = b*t^3/12;
moment = 2*pi*E0*I/a; % final bending moment to make the strip a circle
force  = -2*pi*(t/2)*E0/a;

noSteps     = 20;
deltaF      = force/noSteps;

% build connectivity ...

generateIGA3DMesh

% build visualization mesh only once

buildVisualization3dMesh;

% initialization

u    = zeros(noDofs,1);        % displacement vector
fext = zeros(noDofs,1);        % external force vector


%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*3;

Ke0    = zeros(nElDof,nElDof); % element Ke

%% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';
wFixed     = zeros(size(fixedZNodes))';

udofs      = 3*fixedXNodes-2;
vdofs      = 3*fixedYNodes-1;
wdofs      = 3*fixedZNodes;

plotMesh3(controlPts,weights, uKnot,vKnot,wKnot,...
    p,q,r,40,'r-','try.eps');

%% load-displacement data

rxforce = [0];
rzforce = [0];
rxdisp  = [0];
rzdisp  = [0];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rules
[W,GP]   = quadrature(  noGPs, 'GAUSS', 3 ); % noGPs x noGPs x noGPs quad.
[W1,GP1] = quadrature(  q+1,   'GAUSS', 2 ); % noGPs x noGPs quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

tol      = 1e-5;
iterMax  = 20;  % maximum number of Newton-Raphson iterations

for i = 1:noSteps
    
    disp ('=================================')
    disp (sprintf('   Load step %d',i));
    disp ('=================================')
    disp ('  NR iter : L2-norm residual')
    
    % compute the tangent and internal force
    
    [K,fint] = getTangentStiffnessMatrix3D (u,C);
    
    % compute the external force
    % for simplicity, assume the loaded surface discretized by ONE element

    for gp=1:size(W1,1)
        pt      = GP1(gp,:);
        wt      = W1(gp);
        
        Xi      = parent2ParametricSpace([0 1],pt(1));
        Eta     = parent2ParametricSpace([0 1],pt(2));
        J2      = jacobianPaPaMapping([0 1], [0 1]);
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],q,r,vKnot,wKnot,weights(forcedNodes1));
        
        % updated configuration
        uu      = [u(3*forcedNodes1-2) u(3*forcedNodes1-1) u(3*forcedNodes1)];
        pts     = controlPts(forcedNodes1,:) + uu;
        
        point   = R*pts;
        % basis vectors of the loaded surface
        jacob   = [dRdxi; dRdeta] * pts; % 2x3 matrix
        a1      = jacob(1,:);
        a2      = jacob(2,:);
        a3      = cross(a1,a2); % normal to the surface
        J1      = norm(a3);
                        
        fext(forcedDofsX) = fext(forcedDofsX) + R'*deltaF*pt(2)*a3(1)*J1*J2*wt;
        fext(forcedDofsZ) = fext(forcedDofsZ) + R'*deltaF*pt(2)*a3(3)*J1*J2*wt;
    end
        
    % Newton-Raphson iterations
    error    = 1;
    iiter    = 0;
    
    while error > tol
        
        iiter    = iiter + 1;
        
        % solve for the displacement increment
        
        res      = fext-fint;
        [K,res]  = applyDirichletBCs3D(K,res,udofs,vdofs,wdofs,...
            uFixed,vFixed,wFixed);
        deltaU   = K\res;
        
        % update displacements
        
        u = u + deltaU;
        
        [K,fint] = getTangentStiffnessMatrix3D (u,C);
        
        normf  = norm(fext);
        if normf < 1e-16
            error = norm( fext-fint );
        else
            error = norm( fext-fint ) / normf;
        end
        error = norm(deltaU);
        
        disp (sprintf(' %s %i %s %5.4e ', 'Iter',iiter, ':', error) );
        
        if iiter == iterMax
            error('Newton-Raphson iterations did not converge!');
        end
    end
    
    %% the step has been converged...
    
    rxforce  = [rxforce fint(recordedPtXDof)];
    rzforce  = [rzforce fint(recordedPtZDof)];
    rxdisp   = [rxdisp  u(recordedPtXDof)];
    rzdisp   = [rzdisp  u(recordedPtZDof)];
    
    % write to VTK files
    
    vtuFile = sprintf('../results/%s%d',vtuFileName,i);
    
    writeVTK3D(vtuFile,u,node,elementV);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% write a PVD file that collects all VTU files
% An example is given
% <VTKFile byte_order='LittleEndian' type='Collection' version='0.1'>
% <Collection>
% <DataSet file='cantilever8-0.vtu' groups='' part='0' timestep='0'/>
% <DataSet file='cantilever8-1.vtu' groups='' part='0' timestep='1'/>
% </Collection>
% </VTKFile>

pvdFile = fopen(strcat('../results/',vtuFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:noSteps
    vtuFile = sprintf('%s%d%s',vtuFileName,i,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);
%% plot the load-displacement curve

figure
hold on
plot(-rxdisp,t*rxforce,'bs-','LineWidth',1.2);
plot( rzdisp,t*rxforce,'rs-','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('moment (Nmm)');
%axis([0 8 0 0.27])

%% analytical solution
deltaM = moment/40;
uxx = [];
uzz = [];
mom = [];
for i=0:40
    mo    = i*deltaM;
    mom   = [mom mo];
    theta = mo*a/E0/I;
    uval  = a - a/theta*tan(theta/2)*(1+cos(theta));
    zval  =     a/theta*tan(theta/2)*sin(theta);
    uxx   = [uxx uval];
    uzz   = [uzz zval];
end

figure
hold on
plot( uxx,mom,'bs-','LineWidth',1.2);
plot( uzz,mom,'rs-','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('moment (Nmm)');
legend('horizontal displacement','vertical displacement')
%%

momn=[0];
for i=1:noSteps
    momn = [momn i*moment/noSteps];
end

figure
hold on
plot(-rxdisp,momn,'b-','LineWidth',1.2);
plot( rzdisp,momn,'r-','LineWidth',1.2);
plot( uxx,mom,'bs','LineWidth',1.2);
plot( uzz,mom,'rs','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('moment (Nmm)');

