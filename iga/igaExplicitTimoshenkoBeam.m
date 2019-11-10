%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Transient analysis of the Timoshenko beam.
% Time integration: explicit central difference.
%
% Vinh Phu Nguyen,
% March 2013
% Cardiff University, Wales, UK
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../geometric-nonlinear/

tic

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho C noDofs noCtrPts elConnU elConnV Ke0 Me0 damping I
global noPtsX noPtsY

E0          = 3e7;  % Young modulus N/m^2
nu0         = 0.3;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
rho         = 1; % kg/m^3; 1N = 1 kg * m/s^2
damping     = 0.0;

vtuFileName = '../results/dynamicTimoBeam';

% Newmark integration constants

delta = 0.5;
beta  = 0.25;
dTime = 1e-4;
tMax  = 2;     % in seconds

wf=27;


% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timoshenkoBeamCkData

noGPs  = p+1;
noGPs1 = q+3;

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

eps=1e-12;
fixedNodes  = find(abs(controlPts(:,1))<eps)';
forcedNodes = find(abs(controlPts(:,1)-a)<eps)';

pointA  = intersect (find(abs(controlPts(:,1)-48)<eps),...
                     find(abs(controlPts(:,2))<eps));

%pointA= noCtrPts;

adof = 2*pointA;                 
                 
I = (1/12)*b^3;

% build connectivity ...

generateIGA2DMesh

% build a 1D mesh for the right edge over which
% a traction is applied

bndPoints   = controlPts(forcedNodes,:);

bndMesh = zeros(noElemsV,q+1);

for i=1:noElemsV
    bndMesh(i,:) = forcedNodes(i:i+q);
end

% initialization

u = zeros(noDofs,1);        % displacement vector
v = zeros(noDofs,1);        % velocity vector
a = zeros(noDofs,1);        % acceleration vector

% essential boundary conditions,

udofs      = 2*fixedNodes-1;
vdofs      = 2*fixedNodes;

uFixed     = zeros(size(fixedNodes));
vFixed     = zeros(size(fixedNodes));

dispA = [0];
time  = [0];
f     = [0];

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*2;

Ke0    = zeros(nElDof,nElDof); % element Ke
Me0    = zeros(nElDof,nElDof); % element Me

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule

[W,Q]   = quadrature( noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature
[W1,Q1] = quadrature( noGPs1,'GAUSS', 1 );

buildVisualizationMesh;

%% compute initial acceleration


fext0   = getExternalForce(bndMesh,bndPoints,b,wf,0);
[K,M,Cc] = getMatricesWithLumping();    
    
a        = fext0./M;

a(udofs) = 0;
a(vdofs) = 0;


%% loop on time steps
t = 0;

id = 0;

a0 = a;

for i = 0:100000
    
    disp ('=================================')
    disp (sprintf('   Load step %d',i));
    disp ('=================================')
    disp ('  NR iter : L2-norm residual')
    
    t = i * dTime;
        
    % compute the external force
    
    fext = getExternalForce(bndMesh,bndPoints,b,wf,t);
    
    %fext(adof) = -p0*sin(wf*t);
   
    % solve for the new acceleration
        
    fmod     = fext-K*(u + dTime*v + 0.5*dTime^2*a);              
    a        = fmod./M;
    
    a(udofs) = 0;
    a(vdofs) = 0;
    
    % update displacements and velocity
    
    u   = u + dTime*v + 0.5*dTime^2*a0;
    v   = v + 0.5*dTime*(a + a0);
    
    time  = [time; t];
    dispA = [dispA; u(adof)];
    f     = [f;fext(adof)];
    
    % for every 20 steps write a VTK file
    
%     if rem(i,20) == 0
%         id = id + 1;
%         vtuFile = sprintf('../results/%s%d',vtuFileName,id);
%         
%         writeVTK(vtuFile,u,node,elementV);
%     end
    
    a0 = a;
    
    if t >= tMax
        disp (['  Simulation has been finished in ', num2str(toc), ' seconds.'])
        break;
    end
end

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

for i = 1:id
    vtuFile = sprintf('%s%d%s',vtuFileName,i,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
figure
plot(time,dispA,'r-','LineWidth',1.0);
xlabel('time')
ylabel('Vertical displacement v');
legend('$\Delta t=1e-3, c=0$');
%%

%plot(time,f,'r-','LineWidth',1.1);







