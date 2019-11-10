%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for two dimensional elasticity problems.
%
% Transient analysis of the Timoshenko beam.
% Time integration: implicit Newmark (a-form of Prof. Hughes).
%
% The applied force is:
%
% f = -1000 * gt /(2*I)*((D*D)/4-x(1,2)^2);
% in which gt=g(t)=sin(omega*t) or any time dependent function.
%
% Function handle @ is used to pass the dynamic function gt.
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

tic 

clc
clear all

global element controlPts p q uKnot vKnot weights index elRangeU elRangeV
global Q W Q1 W1 rho De noDofs noCtrPts elConnU elConnV Ke0 Me0 damping I
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
dTime = 0.005;
tMax  = 2;     % in seconds

wf=27;

% handle to the dynamic function
% choose one of the following

gt       = @(t) (sin(wf*t)); % harmonic function
%gt       = @(t) (t <= 0.5);  % step function


% COMPUTE ELASTICITY MATRIX

De = elasticityMatrix(E0,nu0,stressState);

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


fext0   = getExternalForce(bndMesh,bndPoints,b,gt,0);
[K,M,Cc] = getMatrices();    
    
[M,fext0] = applyDirichletBCs(M,fext0,udofs,vdofs,uFixed,vFixed);
a         = M\fext0;

% stiffness, mass and damping matrices are constant

[K,M,Cc] = getMatrices();

Kmod    = M + delta*dTime*Cc + beta*dTime^2*K;

%% loop on time steps
t = 0;

id = 0;
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
    
    fext = getExternalForce(bndMesh,bndPoints,b,gt,t);
    
    %fext(adof) = -p0*sin(wf*t);
    
    % solve for the new acceleration
        
    fmod     = fext-K*predDisp - Cc*predVelo;    
    
    [Kmod,fmod] = applyDirichletBCs(Kmod,fmod,udofs,vdofs,uFixed,vFixed);    
    a           = Kmod\fmod;
    
    % update displacements and velocity
    
    u   = predDisp + beta  * dTime^2 * a;
    v   = predVelo + delta * dTime   * a;
    

    xi  = 1;
    eta = 0.5;
    
    Ux    = u(1:noCtrPts);
    Uy    = u(noCtrPts+1:noDofs);

    projcoord = nurb2proj(noPtsX*noPtsY, [Ux Uy], weights);
    dim       = size(projcoord,2);
    tem       = SurfacePoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                   projcoord,dim,xi,eta);
    %format long
    numUy = tem(2)/tem(3);
    numUx = tem(1)/tem(3);


    time  = [time; t];
    dispA = [dispA; u(adof)];
    f     = [f;fext(adof)];
    
    % for every 20 steps write a VTK file
    
    if rem(i,20) == 0
        id = id + 1;
        vtuFile = sprintf('../results/%s%d',vtuFileName,id);
        
        writeVTK(vtuFile,u,node,elementV);
    end
    
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







