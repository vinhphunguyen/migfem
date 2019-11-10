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
addpath ../post-processing/
addpath ../fem-functions/
addpath ../analytical-solutions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q element index elRangeU elRangeV uKnot vKnot nElmLK...
    weights noCtrPts noDofs Ke0 controlPts W GP  noPtsX noPtsY C

E0          = 1.2e5;  % Young modulus
nu0         = 0.0;  % Poisson ratio
stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS



%% COMPUTE ELASTICITY MATRIX

C = elasticityMatrix(E0,nu0,stressState);

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
forcedNodes  =  find(abs(controlPts(:,1)-a)<1e-8);

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

forcedDofsX      = 2*forcedNodes-1;
forcedDofsY      = 2*forcedNodes;

%% load-displacement data

rxforce=[];
ryforce=[];
rxdisp =[];
rydisp =[];

b = 1.0;
I      = b*t^3/12;
moment = 2*pi*E0*I/a; % final bending moment to make the strip a circle
force  = 100*2*pi*(t/2)*E0/a;

noSteps     = 60;
deltaM      = moment/noSteps;
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
[W,GP]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature
[W1,GP1] = quadrature(2,   'GAUSS', 1 ); % noGPs x noGPs quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

tol      = 1e-6;
iterMax  = 20;  % maximum number of Newton-Raphson iterations

for i = 1:noSteps
    
    disp ('=================================')
    disp (sprintf('   Load step %d',i));
    disp ('=================================')
    disp ('  NR iter : L2-norm residual')
    
    % compute the tangent and internal force
    
    [K,fint] = getTangentStiffnessMatrix (u,C);
    
    deltaF   = -deltaM/I;
    % compute the external force
    for gp=1:size(W1,1)
        pt      = GP1(gp,:);
        wt      = W1(gp);
        
        Xi      = 0.5 * ( ( 1 - 0 ) * pt + 1 + 0);% coord in parameter space
        J2      = 0.5 * ( 1 - 0 );
        [R dRdxi] = NURBS1DBasisDers(Xi,q,vKnot,weights(forcedNodes));
        [R0 dRdxi0] = NURBS1DBasisDers(0.5,q,vKnot,weights(forcedNodes));
        
        % updated configuration
        uu      = [u(2*forcedNodes-1) u(2*forcedNodes)];
        pts     = controlPts(forcedNodes,:) + uu;
        
        point   = R*pts;
        point0  = R0*pts;
        dp      = point - point0;
        yy      = sign(dp(2))*norm(dp)
        % basis vectors of the loaded surface
        jacob   = dRdxi * pts; % 1x2 matrix
        
        J1      = norm(jacob);
        
        fext(forcedDofsX) = fext(forcedDofsX) - R'*deltaF*yy*jacob(2)*J1*J2*wt;
        fext(forcedDofsY) = fext(forcedDofsY) + R'*deltaF*yy*jacob(1)*J1*J2*wt;
    end
    
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
        
        [K,fint] = getTangentStiffnessMatrix (u,C);
        
        % norm of force should taken reaction forces into account
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
    
    rxforce = [rxforce fint(forcedDofsX(1))];
    ryforce = [ryforce fint(forcedDofsY(1))];
    rydisp  = [rydisp u(forcedDofsY(1))];
    rxdisp  = [rxdisp u(forcedDofsX(1))];
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

%% analytical solution
deltaM = moment/noSteps;
uxx = [];
uzz = [];
mom = [];
for i=1:noSteps
    mo    = i*deltaM;
    mom   = [mom mo];
    theta = mo*a/E0/I;
    uval  = a - a/theta*tan(theta/2)*(1+cos(theta));
    zval  =     a/theta*tan(theta/2)*sin(theta);
    uxx   = [uxx uval];
    uzz   = [uzz zval];
end

momn = 0.5*t*sqrt(rxforce.^2+ryforce.^2);

figure
hold on
plot( -rxdisp,mom,'b*-','LineWidth',1.2);
plot( uxx,mom,'bs','LineWidth',1.2);
xlabel('horizontal displacment (mm)');
ylabel('moment (Nmm)');
%axis([0 8 0 7])

figure
hold on
plot( -rydisp,mom,'r*-','LineWidth',1.2);
plot( uzz,mom,'rs','LineWidth',1.2);
xlabel('vertical displacment (mm)');
ylabel('moment (Nmm)');
%axis([0 8 0 7])
