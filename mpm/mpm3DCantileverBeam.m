% This file implements the Material Point Method
% Three dimensional problems.
% The grid is a structured mesh consisting of 8-node brick elements.
%
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% February 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../nurbs-geopdes/inst/
addpath ../delamination/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../post-processing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E0   = 1000;           % Young's modulus
nu0  = 0.3;            % Poisson ratio
rho = 1000;            % density
kappa = 3-4*nu0;       % Kolosov constant
mu    = E0/2/(1+nu0);  % shear modulus
g     = 10;            % gravity mm/s^2

C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
                                  nu0 1-nu0 nu0;
                                  nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

I = [1 0 0;0 1 0; 0 0 1];

vtkFileName  = 'mpm3DBeam';
vtkFileName1  = '../results/mpm3DBeamGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

Lx = 6;
Ly = 3;
Lz = 1;

noX = 4;
noY = 2;
noZ = 0;

[mesh]=buildGrid3D(Lx,Ly,Lz,noX,noY,noZ);


numx = mesh.noElemsX;      % number of elements along X direction
numy = mesh.noElemsY;      % number of elements along Y direction
numz = mesh.noElemsZ;      % number of elements along Z direction

deltax = Lx/numx;
deltay = Ly/numy;
deltaz = Lz/numz;

element = mesh.element;
node    = mesh.node;

elemCount = size(element,1);
nodeCount = size(node,1);


%%   particle distribution 

Lx = 4;
Ly = 1;
Lz = 1;

noX = 3;
noY = 1;
noZ = 0;

pMesh =buildGrid3D(Lx,Ly,Lz,noX,noY,noZ); % particle mesh

pMesh.node(:,2) = pMesh.node(:,2) + 3/2-Ly/2;

ngp = 2;
[W,Q]=quadrature( ngp, 'GAUSS', 3 ); % two point quadrature

pCount = size(pMesh.element,1)*size(W,1);

Mp  = ones(pCount,1);                 % mass
Vp  = ones(pCount,1);                 % volume
Fp  = ones(pCount,9);                 % gradient deformation
s   = zeros(pCount,6);                % stress
eps = zeros(pCount,6);                % strain
vp  = zeros(pCount,3);                % velocity
xp  = zeros(pCount,3);                % position

% initialise particle position, mass ,volume, velocity
% particle taken as global Gauss points

% particles as Gauss points of the particle mesh
id=1;
for e = 1:size(pMesh.element,1)
    sctr = pMesh.element(e,:);
    pts  = pMesh.node(sctr,:);
    for q=1:size(W)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis('B8',pt);  % element shape functions
        J0=dNdxi'*pts;
        detJ0=det(J0);
        a     = wt*detJ0;
        Vp(id) = a;
        Mp(id) = a*rho;
        xp(id,:) =  N'*pts;
        Fp(id,:) = [1 0 0 0 1 0 0 0 1];
        
        id = id + 1;
    end
end

Vp0 = Vp;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    z = xp(p,3);
    e = floor(x/deltax) + 1 + numx*floor(y/deltay) + numx*numy*floor(z/deltaz) ;
    pElems(p) = e;
end

for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,3);  % nodal momentum vector
niforce   = zeros(nodeCount,3);  % nodal internal force vector
neforce   = zeros(nodeCount,3);  % nodal external force vector

%% Dirichlet boundary 

eps0=1e-9;
leftNode  = find( node(:,1) == 0 );

%% plot mesh, particles

hold on
plot_mesh(node,element,'B8','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot3(xp(:,1),xp(:,2),xp(:,3),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E0/rho);
dtime = 0.001;
time  = 0.9;
t     = 0;


nsteps = floor(time/dtime);

interval = 10;

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 0;

while ( t < time )
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            
            % particle mass and momentum to node
            stress = s(pid,:);
            
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - node(id,:);
                [N,dNdx]=getMPM3D(x,deltax,deltay,deltaz);
                dNIdx = dNdx(1);
                dNIdy = dNdx(2);
                dNIdz = dNdx(3);
                nmass(id)       = nmass(id)       + N*Mp(pid);
                nmomentum(id,:) = nmomentum(id,:) + N*Mp(pid)*vp(pid,:);
                
                niforce(id,1)   = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(6)*dNIdy + stress(5)*dNIdz);
                niforce(id,2)   = niforce(id,2) - Vp(pid)*(stress(6)*dNIdx + stress(2)*dNIdy + stress(4)*dNIdz);
                niforce(id,3)   = niforce(id,3) - Vp(pid)*(stress(5)*dNIdx + stress(4)*dNIdy + stress(3)*dNIdz);
                
                neforce(id,2)   = neforce(id,2) - g * N*Mp(pid);
            end
        end
    end
    
    % debug
    
    
    % update nodal momenta
    
    nmomentum= nmomentum + (niforce + neforce)*dtime;
   
    nmomentum(leftNode,:) = 0;
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            
            Lp = zeros(3,3);
            for i=1:length(esctr)
                id = esctr(i);
                x     = xp(pid,:) - node(id,:);
                [N,dNdx]=getMPM3D(x,deltax,deltay,deltaz);
              
                vI = [0 0 0];
                if nmass(id) > tol
                    vp(pid,:)  = vp(pid,:) + dtime * N*(niforce(id,:) + neforce(id,:)) /nmass(id);
                    xp(pid,:)  = xp(pid,:) + dtime * N*nmomentum(id,:)/nmass(id);
                    vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
                end
                Lp = Lp + vI'*dNdx;         % particle gradient velocity
            end
            
            F       = ( I + Lp*dtime)*reshape(Fp(pid,:),3,3);
            Fp(pid,:)= reshape(F,1,9);
            Vp(pid) = det(F)*Vp0(pid);
            dEps    = dtime * 0.5 * (Lp+Lp');
            dsigma  = C * [dEps(1,1);dEps(2,2);dEps(3,3);2*dEps(2,3);2*dEps(1,3); 2*dEps(1,2)] ;
            s(pid,:)  = s(pid,:) + dsigma';
            eps(pid,:)= eps(pid,:) +  [dEps(1,1) dEps(2,2) dEps(3,3) 2*dEps(2,3) 2*dEps(1,3) 2*dEps(1,2)] ;
            
            k = k + 0.5*(vp(pid,1)^2+vp(pid,2)^2)*Mp(pid);
            u = u + 0.5*Vp(pid)*s(pid,:)*eps(pid,:)';
        end
    end
    
    % store time,velocty for plotting
    
    pos{istep+1} = xp;
    vel{istep+1} = vp;
    
    % update the element particle list
    
    for p=1:pCount
        x = xp(p,1);
        y = xp(p,2);
        z = xp(p,3);
        e = floor(x/deltax) + 1 + numx*floor(y/deltay) + numx*numy*floor(z/deltaz);
        pElems(p) = e;
    end
    
    for e=1:elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep,interval) == 0 )
        vtkFile = sprintf('../results/%s%d',vtkFileName,istep);
        VTKParticles(xp,vtkFile,s);
    end
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
end


%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:istep
    if (  mod(i,interval) == 0 )
        vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
        fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
    end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess3d(node,element,'B8',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY sigmaXX sigmaXX sigmaXX],[Ux Uy Ux]);



figure
set(gca,'FontSize',14)
hold on
plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
xlabel('Time')
ylabel('Energy')
legend('kinetic','strain','total')
%set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
%axis([0 3 0 3])

disp([num2str(toc),'   DONE '])
