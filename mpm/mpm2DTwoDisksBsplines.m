% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of B-splines.
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% February 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../delamination/
addpath ../nurbs-geopdes/inst/
addpath ../meshing/
addpath ../C_files/
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

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

interval     = 100;
vtkFileName  = 'mpm2DTwoDisks';

v   = 0.1;

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])

%%   particle distribution from a mesh
%

meshFile = 'disks.msh';
mesh     = load_gmsh (meshFile);

elemType = 'T3';
numnode  = mesh.nbNod;
numelem  = mesh.nbTriangles;
node1    = mesh.POS(:,1:2);
element1 = mesh.TRIANGLES(1:numelem,1:3);

% check if Jacobian is negative

element1  = tricheck(node1,element1,1);

% one particle per triangle

pCount = numelem;

Mp  = ones(pCount,1);                 % mass
Vp  = ones(pCount,1);                 % volume
Fp  = ones(pCount,4);                 % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position

% initialise particle position, mass ,volume, velocity

for e = 1:numelem
    coord = node1(element1(e,:),:);
    a     = det([coord,[1;1;1]])/2;
    Vp(e) = a;
    Mp(e) = a*rho;
    xp(e,:) = mean(coord);
    
    if xp(e,1) < 0.5
        vp(e,:) = [v v];
    else
        vp(e,:) = [-v -v];
    end
    
    Fp(e,:) = [1 0 0 1];
end
Vp0 = Vp;

%% Computational grid

L = 1;
w = 1;

controlPts          = zeros(4,2,2);

controlPts(1:2,1,1) = [0;0];
controlPts(1:2,2,1) = [L;0];
controlPts(1:2,1,2) = [0;w];
controlPts(1:2,2,2) = [L;w];

controlPts(4,:,:)   = 1;

uKnot = [0 0 1 1];
vKnot = [0 0 1 1];

pnew=1+0; % new order basis in x dir
qnew=1+0; % new order basis in y dir

surf    = nrbmak(controlPts,{uKnot vKnot});   % Bspline surface
surf    = doKRefinementSurface(surf,pnew,qnew,3,3); % k-refinement
igaMesh = buildIGA2DMesh(surf);               % parameter mesh
vMesh   = buildVisualizationMesh2D(surf);     % physical mesh

elemCount = igaMesh.elemCount;                % # of elements
nodeCount = length(igaMesh.weights);          % # of grid nodes
deltax    = L/igaMesh.noElemsU;
deltay    = w/igaMesh.noElemsV;
xMin      = 0;
xMax      = L;
yMin      = 0;
yMax      = w;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    e = floor(x/deltax) + 1 + igaMesh.noElemsU*floor(y/deltay);
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

%% plot B-splines mesh, particles

hold on
plot_mesh(node1,element1,elemType,'r-',0.2);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
nrbkntplot(surf);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
critT = deltax/c;
dtime = 0.001;
time  = 0.001;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    % loop over computational cells or elements
    for e=1:elemCount                
        esctr  = igaMesh.globElems (e,:);     % element connectivity
        pts    = igaMesh.controlPts(esctr,:); % element nodal coords
        mpts   = mpoints{e};                  % particles inside element e        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            x    = xp(pid,1);
            y    = xp(pid,2);
            xi   = (x-xMin)/L;
            et   = (y-yMin)/w;
            
            [N, dNdxi, dNdeta] = NURBS2DBasisDers([xi;et],igaMesh.p,igaMesh.q,...
                igaMesh.uKnot,igaMesh.vKnot,igaMesh.weights');
                        
            jacob   = pts' * [dNdxi' dNdeta'];
            J1      = det(jacob);
            dNdx    = [dNdxi' dNdeta'] * inv(jacob);
                                    
            % particle mass and momentum to node
            stress = s(pid,:);
            
            for i=1:length(esctr)
                id    = esctr(i);
                dNIdx = dNdx(i,1);
                dNIdy = dNdx(i,2);
                nmass(id)       = nmass(id)       + N(i)*Mp(pid);
                nmomentum(id,:) = nmomentum(id,:) + N(i)*Mp(pid)*vp(pid,:);
                niforce(id,1)   = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy);
                niforce(id,2)   = niforce(id,2) - Vp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy);
            end
        end
    end
    
    % debug
    
    % update nodal momenta
    for i=1:nodeCount
        nmomentum(i,:) = nmomentum(i,:) + niforce(i,:)*dtime;
    end
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr  = igaMesh.globElems (e,:);     % element connectivity
        pts    = igaMesh.controlPts(esctr,:); % element nodal coords        
        mpts   = mpoints{e};        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            x    = xp(pid,1);
            y    = xp(pid,2);
            xi   = (x-xMin)/L;
            et   = (y-yMin)/w;
            
            [N, dNdxi, dNdeta] = NURBS2DBasisDers([xi;et],igaMesh.p,igaMesh.q,...
                igaMesh.uKnot,igaMesh.vKnot,igaMesh.weights');
                        
            jacob   = pts' * [dNdxi' dNdeta'];
            J1      = det(jacob);
            dNdx    = [dNdxi' dNdeta'] * inv(jacob);
            
            
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                vI = [0 0];
                if nmass(id) > tol
                    vp(pid,:)  = vp(pid,:) + dtime * N(i)*niforce(id,:)  /nmass(id);
                    xp(pid,:)  = xp(pid,:) + dtime * N(i)*nmomentum(id,:)/nmass(id);
                    vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
                end
                Lp = Lp + vI'*dNdx(i,:);         % particle gradient velocity
            end
            
            F       = ([1 0;0 1] + Lp*dtime)*reshape(Fp(pid,:),2,2);
            Fp(pid,:)= reshape(F,1,4);
            Vp(pid) = det(F)*Vp0(pid);
            dEps    = dtime * 0.5 * (Lp+Lp');
            dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
            s(pid,:)  = s(pid,:) + dsigma';
            eps(pid,:)= eps(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
            
            k = k + 0.5*(vp(pid,1)^2+vp(pid,2)^2)*Mp(pid);
            %             u = u + 0.25/mu*( 0.25*(kappa+1)*(s(pid,1)^2+s(pid,2)^2) ...
            %                   - 2*(s(pid,1)*s(pid,2)-s(pid,3)^2))*Vp(pid);
            u = u + 0.5*Vp(pid)*s(pid,:)*eps(pid,:)';
        end
    end
    
    % store time,velocty for plotting
    
    pos{istep} = xp;
    vel{istep} = vp;
    
    % update the element particle list
    
    for p=1:pCount
        x = xp(p,1);
        y = xp(p,2);
        e = floor(x/deltax) + 1 + igaMesh.noElemsU*floor(y/deltay);
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
    
    if (  mod(istep-1,interval) == 0 )
        xp = pos{istep};        
        vtkFile = sprintf('../results/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,s);
    end
    
        
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end
a=0;
for e=1:elemCount
    pts=mpoints{e};
    a = a + length(pts);
end



%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../results/',vtkFileName,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:nsteps
     if (  mod(i,interval) == 0 )
    vtuFile = sprintf('%s%d%s',vtkFileName,i,'.vtp');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''0'' timestep=''%d''/>\n',vtuFile,i);
     end
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);

% interval = 25;
% m = 1;
% clf;
% for it=1:nsteps
%     if (  mod(it,interval) == 0 )
%         xp = pos{it};
%         vp = vel{it};
%         figure(m)
%         hold on
%         plot_mesh(node,element,'Q4','k-',1.);
%         %scatter(xp(:,1),xp(:,2),40,sqrt(vp(:,1).^2+vp(:,2).^2),'full');
%         plot(xp(:,1),xp(:,2),'k.','markersize',10);
%         %plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
%         axis equal
%         axis off
%         title( sprintf( '%s: %f', 'time', ta(it) ) );
%         %colorbar
%         %print(m, '-djpeg90', ['disks-',num2str(m),'.jpg']);
%         exportfig(gcf,sprintf( '%s %d %s', 'disks', m, '.eps' ),opts)
%         m = m + 1;
%     end
% end


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
axis([0 3 0 3])

disp([num2str(toc),'   DONE '])
