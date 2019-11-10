%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
% Geometrical nonlinearity is allowed. However strain is small.
% Total Lagrangian formulation is used.
%
% The nodal displacement are stored as [ux1 uy1 ux2 uy2 ...]
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
addpath ../integration/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q element index elRangeU elRangeV uKnot vKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP  noPtsX noPtsY De

E0          = 100;  % Young modulus
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS

deltaF      = -0.01;
noSteps     = 20;

%% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input file

cantileverBeamData

vtuFileName = 'cantileverBeam';

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;
noGPs          = p+1;

% find boundary nodes for boundary conditions

leftNodes    =  find(controlPts(:,1)==0);
fixedXNodes  =  leftNodes;
fixedYNodes  =  leftNodes;
forcedNodes  =  2*noCtrPts;

% build connectivity ...

generateIGA2DMesh

% build visualization mesh only once

buildVisualizationMesh;

% initialization

u    = zeros(noDofs,1);        % displacement vector
fext = zeros(noDofs,1);        % external force vector


% values, row indices, columns indices of the global K matrix
% vSprGK = zeros(nSprGK,1);
% jSprRw = zeros(nSprGK,1);
% jSprCl = zeros(nSprGK,1);

Ke0    = zeros(size(element,2)*2,size(element,2)*2); % element Ke

%% essential boundary conditions

uFixed     = zeros(size(fixedXNodes))';
vFixed     = zeros(size(fixedYNodes))';

udofs      = 2*fixedXNodes-1;
vdofs      = 2*fixedYNodes;

%% load-displacement data

rforce=[0];
rdisp =[0];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,GP]=gaussianQuad2DNURBS(p+1,q+1);

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
    
    [K,fint] = getTangentStiffnessMatrix (u,De);
    
    % compute the external force
    
    fext(forcedNodes) = fext(forcedNodes) + deltaF;
    
    error    = 1;
    iiter    = 0;
    
    while error > tol
        
        iiter    = iiter + 1;
        
        % solve for the displacement increment
        
        res      = fext-fint;
        [K,res]  = applyDirichletBCs(K,res,udofs,vdofs,uFixed,vFixed);
        deltaU   = K\res;
        
        % update displacements
        
        u = u + deltaU;
        
        [K,fint] = getTangentStiffnessMatrix (u,De);
        
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
    
    rforce = [rforce fint(forcedNodes)];
    rdisp  = [rdisp u(forcedNodes)];
    
    % write to VTK files
    
    vtuFile = sprintf('../results/%s%d',vtuFileName,i);
    
    writeVTK(vtuFile,u,node,elementV);
    
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
plot(-rdisp,-rforce,'bs-','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('force (N)');
%axis([0 7 0 0.2])

%% write the reference solution to a file 
% (disp, force) format

file = fopen('~/code/jive/bezier/large-displacement/matlab.dat', 'wt');

for i=1:length(rdisp)
   fprintf(file, '  %2.8f %2.8f', -rdisp(i), -rforce(i));
   fprintf(file, '\n');
end
fclose(file);
