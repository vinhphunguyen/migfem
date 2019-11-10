%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Stress wave in 1D bar modelled by 2D plane stress elements.
% Time integration: implicit Newmark (a-form of Prof. Hughes).
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
addpath ../nurbs-geopdes/inst/

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho C noDofs noCtrPts elConnU elConnV Ke0 Me0 damping noPtsX noPtsY

E0          = 49;  % Young modulus N/m^2
nu0         = 0.0;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
rho         = 0.01; % kg/m^3; 1N = 1 kg * m/s^2
damping     = 0.0;

% Newmark integration constants

delta = 0.5;
beta  = 0.25;
dTime = 1e-3;
tMax  = 2;     % in seconds

f0 = -0.1;

vtuFileName = '../results/stressWave1D';

% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oneDBar

noGPs  = p+1;
noGPs1 = q+3;

noCtrPts       = noPtsX * noPtsY;
noDofs         = noCtrPts * 2;

% find boundary nodes for boundary conditions

eps=1e-12;
fixedNodes  = find(abs(controlPts(:,1))<eps)';
forcedNodes = find(abs(controlPts(:,1)-a)<eps)';


pointA= noCtrPts;

adof = 2*pointA-1;

% build connectivity ...

generateIGA2DMesh

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
stress     = [0];

%% one the assembly
nElNod = size(element,2);
nElDof = nElNod*2;

Ke0    = zeros(nElDof,nElDof); % element Ke
Me0    = zeros(nElDof,nElDof); % element Me

% build visualization mesh
buildVisualizationMesh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule

[W,Q]   = quadrature( noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature
[W1,Q1] = quadrature( noGPs1,'GAUSS', 1 );

%% compute initial acceleration

fext=zeros(noDofs,1);

fext(2*forcedNodes-1) = (1/noPtsY)*f0;

[K,M,Cc] = getMatrices();

[M,fext] = applyDirichletBCs(M,fext,udofs,vdofs,uFixed,vFixed);
a        = M\fext;

% stiffness, mass and damping matrices are constant

[K,M,Cc] = getMatrices();

Kmod    = M + delta*dTime*Cc + beta*dTime^2*K;

%% loop on time steps
t = 0;
id=0;
for i = 0:100000
    
    disp ('=================================')
    disp (sprintf('   Load step %d',i));
    disp ('=================================')
    disp ('  NR iter : L2-norm residual')
    
    t = i * dTime;
    
    % compute predictor velocities and displacements
    
    predDisp = u + dTime*v + dTime^2*0.5*(1-2*beta)*a;
    predVelo =           v + dTime*(1-delta)*a;
    
    % compute the external force
    
    % solve for the new acceleration
    
    fmod     = fext-K*predDisp - Cc*predVelo;
    
    [Kmod,fmod] = applyDirichletBCs(Kmod,fmod,udofs,vdofs,uFixed,vFixed);
    a           = Kmod\fmod;
    
    % update displacements and velocity
    
    u   = predDisp + beta  * dTime^2 * a;
    v   = predVelo + delta * dTime   * a;
    
    % compute some quantities for postprocessing
    
    sigma  = getStressAt(noElems/2,u);
    
    time   = [time; t];
    dispA  = [dispA; u(adof)];
    stress = [stress;sigma(1)];
    
    % for every 20 steps write a VTK file
    
    if rem(i,20) == 0
        id = id + 1;
        vtuFile = sprintf('../results/%s%d',vtuFileName,id);
        
        writeVTK(vtuFile,u,node,elementV);
    end
    
    % stop the simulation
    
    if t >= tMax
        disp ('  Simulation has been finished.')
        break;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
figure
plot(time,dispA,'r-','LineWidth',1.1);
xlabel('time')
ylabel('Vertical displacement v');
legend('$\Delta t=1e-3, c=0$');
%%

figure
plot(time,stress,'r-','LineWidth',1.1);
xlabel('time')
ylabel('Vertical displacement v');
legend('$\Delta t=1e-3, c=0$');
%%



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

%% save some data to MAT file

savefile = 'stresswave1D2.mat';
save(savefile, 'stress', 'time');

% load that file later
S = load(savefile); 


