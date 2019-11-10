% This file implements the Generalized Interpolation Material Point Method
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-node elements.
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

E0   = 3e5;           % Young's modulus
nu0  = 0.3;            % Poisson ratio
rho = 1000;            % density
kappa = 3-4*nu0;       % Kolosov constant
mu    = E0/2/(1+nu0);  % shear modulus
g     = 100;            % gravity mm/s^2

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E0,nu0,stressState);

I = [1 0;0 1];

vtkFileName  = 'gimp2DBeam';
vtkFileName1  = '../results/gimp2DBeamGrid';

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

Lx = 9;
Ly = 3;

noX = 36;
noY = 12;

[mesh]=buildGrid2D(Lx,Ly,noX,noY);


numx = noX;      % number of elements along X direction
numy = noY;      % number of elements along Y direction


deltax = Lx/numx;
deltay = Ly/numy;


element = mesh.element;
node    = mesh.node;

elemCount = size(element,1);
nodeCount = size(node,1);


%%   particle distribution 

Lx = 8;
Ly = 1;

noX = 32;
noY = 4;


pMesh =buildGrid2D(Lx,Ly,noX,noY); % particle mesh

pMesh.node(:,2) = pMesh.node(:,2) + 3/2-Ly/2;

ngp = 3;
[W,Q]=quadrature( ngp, 'GAUSS', 2 ); % two point quadrature



pCount = size(pMesh.element,1)*size(W,1);

Mp  = ones(pCount,1);                 % mass
Vp  = ones(pCount,1);                 % volume
Fp  = ones(pCount,4);                 % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position
lp  = zeros(pCount,2);                % particle size (used in GIMP basis functions)

lpx0 = deltax/ngp;
lpy0 = deltay/ngp;

lp(:,1) = deltax/ngp;
lp(:,2) = deltay/ngp;

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
        [N,dNdxi]=lagrange_basis('Q4',pt);  % element shape functions
        J0=dNdxi'*pts;
        detJ0=det(J0);
        a     = wt*detJ0;
        Vp(id) = a;
        Mp(id) = a*rho;
        xp(id,:) =  N'*pts;
        Fp(id,:) = [1 0 0 1];
        
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
    e = floor(x/deltax) + 1 + numx*floor(y/deltay);
    pElems(p) = e;
end

for e=1:elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector

%% find GIMP element connectivity

gimpElement = cell(elemCount,1);

for e=1:elemCount
    neighbors = getNeighbors(e, numx, numy);
    temp = [];
    for ie=1:length(neighbors)
        neighborId = neighbors(ie);
        neighborNd = element(neighborId,:);
        temp = [temp neighborNd];
    end
    gimpElement{e} = unique(temp);
end

%% Dirichlet boundary 

eps0=1e-9;
leftNode  = find( node(:,1) == 0 );

%% plot mesh, particles

hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-1; % mass tolerance

c     = sqrt(E0/rho);
dtime = 0.8*deltax/c;
time  = 2;
t     = 0;


nsteps = floor(time/dtime);

interval = 10;

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 0;

    for e=1:elemCount
        esctr = gimpElement{e};    % GIMP (extended) element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        
        if isempty(mpts) continue; end % skip element without particles
         
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            lpx = lp(pid,1);
            lpy = lp(pid,2);
            
            % particle mass and momentum to node
            stress = s(pid,:);
            
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - node(id,:);
               [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
 
                neforce(id,2)   = neforce(id,2) - g * N*Mp(pid);
            end
        end
    end
    
while ( t < time )
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    
    if ( t ~= 0 ) neforce(:)   = 0; end
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = gimpElement{e};    % GIMP (extended) element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            lpx = lp(pid,1);
            lpy = lp(pid,2);
            % particle mass and momentum to node
            stress = s(pid,:);
            
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - node(id,:);
                [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
                dNIdx = dNdx(1);
                dNIdy = dNdx(2);

                nmass(id)       = nmass(id)       + N*Mp(pid);
                nmomentum(id,:) = nmomentum(id,:) + N*Mp(pid)*vp(pid,:);
                
                niforce(id,1)   = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy );
                niforce(id,2)   = niforce(id,2) - Vp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy );
                
                %neforce(id,2)   = neforce(id,2) - g * N*Mp(pid);
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
        esctr = gimpElement{e};    % GIMP (extended) element connectivity
        enode = node(esctr,:);
        mpts  = mpoints{e};
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            lpx = lp(pid,1);
            lpy = lp(pid,2);
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                x     = xp(pid,:) - node(id,:);
                [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
              
                vI = [0 0];
                if nmass(id) > tol
                    vp(pid,:)  = vp(pid,:) + dtime * N*(niforce(id,:) + neforce(id,:)) /nmass(id);
                    xp(pid,:)  = xp(pid,:) + dtime * N*nmomentum(id,:)/nmass(id);
                    vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
                end
                Lp = Lp + vI'*dNdx;         % particle gradient velocity
            end
            
            F       = ( I + Lp*dtime)*reshape(Fp(pid,:),2,2);
            Fp(pid,:)= reshape(F,1,4);
            Vp(pid) = det(F)*Vp0(pid);
            dEps    = dtime * 0.5 * (Lp+Lp');
            dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
            s(pid,:)  = s(pid,:) + dsigma';
            eps(pid,:)= eps(pid,:) +  [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
            
%             lp(pid,1) = lpx0 * F(1,1);
%             lp(pid,2) = lpy0 * F(2,2);
            
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
        e = floor(x/deltax) + 1 + numx*floor(y/deltay);
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

VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY],[Ux Uy]);



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
