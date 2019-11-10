%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
% Geometrical nonlinearity is allowed. However strain is small.
% Total Lagrangian formulation is used.
%
% The nodal displacement are stored as [ux1 uy1 uz1 ux2 uy2 uz3...]
%
% Cantilever beam subject to shear force undergoes large displacements.
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

% input file

cantileverBeamShearData

%%

E0          = 1.2e6;  % Young modulus
nu0         = 0.0;  % Poisson ratio



%% COMPUTE ELASTICITY MATRIX

C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);


tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



vtuFileName = 'cantileverBeamShearGNL';

noCtrPts       = noPtsX * noPtsY *noPtsZ;
noDofs         = noCtrPts * 3;
noGPs          = p+1;

% find nodes on which Dirichlet BCs are to be imposed

eps = 1e-14;

leftNodes    =  find(controlPts(:,1)==0);
fixedXNodes  =  leftNodes;
fixedYNodes  =  leftNodes;
fixedZNodes  =  leftNodes;

% forced nodes
EPS=1e-8;
forcedNodes1  =  find(abs(controlPts(:,1)-a)<EPS);
forcedNodes2  =  find(abs(controlPts(:,2)-b)<EPS);
forcedNodes3  =  find(abs(controlPts(:,2))  <EPS);
%forcedNodesP  =  intersect(forcedNodes1,forcedNodes3);
forcedNodesM  =  intersect(forcedNodes1,forcedNodes2);

% recorded point for load-displacement curve plot
recordedPnt    = forcedNodesM(1);
recordedPtYDof = 3*recordedPnt-1; % y-direction dof
forcedDofs     = 3*forcedNodesM-1; % y-direction dof

% force, no of steps

F           = 4;
noSteps     = 20;
deltaF      = F/noSteps;

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

ryforce = [0];
rydisp  = [0];


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
    
    fext(forcedDofs) = fext(forcedDofs) + deltaF;
        
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
    
    ryforce  = [ryforce fint(recordedPtYDof)];    
    rydisp   = [rydisp  u(recordedPtYDof)];    
    
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

% refenrece solution from Sze et al, 2004
% Popular benchmark problems for geometric nonlinear analysis of shells

P= F*[ 0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4;0.45;0.5;0.55;0.6;0.65;0.7;...
    0.75;0.8;0.85;0.9;0.95;1];

w=[0.663;1.309;1.922;2.493;3.015;3.488;3.912;...
    4.292;4.631;4.933;5.202;5.444;5.660;5.855;6.031;6.190;6.335;6.467;6.588;6.698];

figure
hold on
plot(rydisp,ryforce,'b-','LineWidth',1.2);
plot(w,P,'rs','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('End force (N)');
%axis([0 8 0 0.27])



%% write the reference solution to a file 
% (disp, force) format

file = fopen('~/code/jive/bezier/large-displacement/beam3D-ref.dat', 'wt');

for i=1:length(w)
   fprintf(file, '  %2.8f %2.8f %2.8f %2.8f', w(i), P(i),rydisp(i),ryforce(i));
   fprintf(file, '\n');
end
fclose(file);

