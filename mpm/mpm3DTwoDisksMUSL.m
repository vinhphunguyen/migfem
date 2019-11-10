% This file implements the Material Point Method of Sulsky 1994.
% Three dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% USF formulation and the contact-release algorithm of York.
% Use data structures suitable for contact bodies.
% The shape functions are written in the global coords so that one does not
% have to convert particle position to natural coordinates.
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% March 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/
addpath ../post-processing/
addpath ../delamination/
addpath ../meshing/

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E0   = 1000;        % Young's modulus
nu0  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu0;    % Kolosov constant
mu    = E0/2/(1+nu0);% shear modulus

v   = 0.1;

C=zeros(6,6);
C(1:3,1:3)=E0/(1+nu0)/(1-2*nu0)*[ 1-nu0 nu0 nu0;
    nu0 1-nu0 nu0;
    nu0 nu0 1-nu0];
C(4:6,4:6)=E0/2/(1+nu0)*eye(3);

I = [1 0 0;0 1 0; 0 0 1];

vtkFileName  = 'mpm3DTwoDisk';
vtkFileName1  = '../results/mpm3DTwoDisksGrid';
interval     = 1;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

Lx = 1;
Ly = 1;
Lz = 1;

noX = 4;
noY = 4;
noZ = 4;

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

%% generate material points

[W,Q]=quadrature(  2, 'GAUSS', 3 ); % 2x2 Gaussian quadrature


% body1

center = [0.2 0.2 0.5];
radius = 0.2;


volume = [];
mass   = [];
coord  = [];


for e=1:elemCount                          % start of element loop
    sctr = element(e,:);                   % element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                      % quadrature loop
        pt=Q(q,:);                         % quadrature point
        wt=W(q);                           % quadrature weight
        [N,dNdxi]=lagrange_basis('B8',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        r  = norm(x-center);
        if ( r-radius < 0 )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

particles1.volume = volume;
particles1.volume0 = volume;
particles1.mass   = mass;
particles1.coord  = coord;
particles1.deform = repmat([1 0 0 0 1 0 0 0 1],pCount,1);     % gradient deformation
particles1.stress = zeros(pCount,6);                % stress
particles1.strain = zeros(pCount,6);                % strain
particles1.velo   = zeros(pCount,3) ;               % velocity
particles1.velo(:,1:2) = v;


% body 2
center = [0.8 0.8 0.5];
radius = 0.2;

volume = [];
mass   = [];
coord  = [];

for e=1:elemCount                          % start of element loop
    sctr = element(e,:);                   %  element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                      % quadrature loop
        pt=Q(q,:);                         % quadrature point
        wt=W(q);                           % quadrature weight
        [N,dNdxi]=lagrange_basis('B8',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        r  = norm(x-center);
        if ( r-radius < 0 )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

particles2.volume  = volume;
particles2.volume0 = volume;
particles2.mass   = mass;
particles2.coord  = coord;
particles2.deform = repmat([1 0 0 0 1 0 0 0 1],pCount,1);     % gradient deformation
particles2.stress = zeros(pCount,6);                % stress
particles2.strain = zeros(pCount,6);                % strain
particles2.velo   = zeros(pCount,3) ;               % velocity
particles2.velo(:,1:2) = -v;

body1.particles = particles1;
body2.particles = particles2;

bodies{1} = body1;
bodies{2} = body2;


bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles


for ib=1:length(bodies)
    body      = bodies{ib};
    particles = body.particles;
    elems     = ones(length(particles.volume),1);
    
    for ip=1:length(particles.volume)
        x = particles.coord(ip,1); y = particles.coord(ip,2); z = particles.coord(ip,3);
        e = floor(x/deltax) + 1 + numx*floor(y/deltay) + numx*numy*floor(z/deltaz) ;
        elems(ip) = e;
    end
    
    bodies{ib}.elements = unique(elems);
    bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
    mpoints = cell(elemCount,1);
    for ie=1:elemCount
        id  = find(elems==ie);
        mpoints{ie}=id;
    end
    
    
    bodies{ib}.mpoints  = mpoints;
end


%% plot mesh, particles

coords=[bodies{1}.particles.coord;bodies{2}.particles.coord];
hold on
plot_mesh(node,element,'B8','k-',0.7);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot3(coords(:,1),coords(:,2),coords(:,3),'k.','markersize',20);
axis off

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,3);  % nodal momentum vector
niforce   = zeros(nodeCount,3);  % nodal internal force vector
neforce   = zeros(nodeCount,3);  % nodal external force vector
nvelo     = zeros(nodeCount,3*(bodyCount+1)); % nodal velocities
nvelo0    = zeros(nodeCount,3*(bodyCount+1));;
nacce     = zeros(nodeCount,3);  % nodal acceleration vector

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E0/rho);
dtime = 0.001;
time  = 2,5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    %% reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    
    %% loop over bodies (update nodal momenta without contact)
    
    for ib=1:bodyCount
        body      = bodies{ib};
        particles = body.particles;
        elems     = body.elements;
        mpoints   = body.mpoints;
        % loop over computational cells or elements
        for ie=1:length(elems)
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};       % particles inside element e
            
            % loop over particles
            for p=1:length(mpts)
                pid    = mpts(p);
                xp     = particles.coord(pid,:);
                stress = particles.stress(pid,:);
                Mp     = particles.mass(pid);
                vp     = particles.velo(pid,:);
                Vp     = particles.volume(pid);
                % loop over nodes of current element "ie"
                for i=1:length(esctr)
                    id    = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getMPM3D(x,deltax,deltay,deltaz);
                    
                    dNIdx = dNdx(1);
                    dNIdy = dNdx(2);
                    dNIdz = dNdx(3);
                    nmass(id)       = nmass(id)       + N*Mp;
                    nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
                    
                    niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(6)*dNIdy + stress(5)*dNIdz);
                    niforce(id,2)   = niforce(id,2) - Vp*(stress(6)*dNIdx + stress(2)*dNIdy + stress(4)*dNIdz);
                    niforce(id,3)   = niforce(id,3) - Vp*(stress(5)*dNIdx + stress(4)*dNIdy + stress(3)*dNIdz);
                end
            end
        end
        
        % update nodal momenta
        
        activeNodes=bodies{ib}.nodes;
        nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
        
        %         nvelo(activeNodes,2*ib-1) = nmomentum(activeNodes,1)./nmass(activeNodes);
        %         nvelo(activeNodes,2*ib)   = nmomentum(activeNodes,2)./nmass(activeNodes);
        %         nacce(activeNodes,2*ib-1) = niforce(activeNodes,1)./nmass(activeNodes);
        %         nacce(activeNodes,2*ib)   = niforce(activeNodes,2)./nmass(activeNodes);
    end
    
    %% update particle velocity and position and stresses
    k = 0; u = 0;
    for ib=1:bodyCount
        body      = bodies{ib};
        particles = body.particles;
        elems     = body.elements;
        mpoints   = body.mpoints;
        % loop over computational cells or elements
        for ie=1:length(elems)
            e     = elems(ie);
            esctr = element(e,:);      % element connectivity
            enode = node(esctr,:);     % element node coords
            mpts  = mpoints{e};       % particles inside element e
            % loop over particles
            for p=1:length(mpts)
                pid  = mpts(p);
                xp   = particles.coord(pid,:);
                xp0  = particles.coord(pid,:);
                Mp   = particles.mass(pid);
                vp   = particles.velo(pid,:);
                Vp   = particles.volume(pid);
                
                Lp   = zeros(3,3);
                for i=1:length(esctr)
                    id = esctr(i);
                    x     = xp0 - node(id,:);
                    [N,dNdx]=getMPM3D(x,deltax,deltay,deltaz);
                    
                    vI = [0 0 0];
                    if nmass(id) > tol
                        vp  = vp + dtime * N*(niforce(id,:) + neforce(id,:)) /nmass(id);
                        xp  = xp + dtime * N*nmomentum(id,:)/nmass(id);
                        vI         = nmomentum(id,:)/nmass(id);  % nodal velocity
                    end
                    Lp = Lp + vI'*dNdx;         % particle gradient velocity
                end
                
                bodies{ib}.particles.velo(pid,:) = vp;
                bodies{ib}.particles.coord(pid,:)= xp;
                F       = (I+ Lp*dtime)*reshape(bodies{ib}.particles.deform(pid,:),3,3);
                bodies{ib}.particles.deform(pid,:) = reshape(F,1,9);
                bodies{ib}.particles.volume(pid  ) = det(F)*bodies{ib}.particles.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                dsigma  = C * [dEps(1,1);dEps(2,2);dEps(3,3);2*dEps(2,3);2*dEps(1,3); 2*dEps(1,2)] ;
                bodies{ib}.particles.stress(pid,:)  = bodies{ib}.particles.stress(pid,:) + dsigma';
                bodies{ib}.particles.strain(pid,:)  = bodies{ib}.particles.strain(pid,:) + ...
                   [dEps(1,1) dEps(2,2) dEps(3,3) 2*dEps(2,3) 2*dEps(1,3) 2*dEps(1,2)] ;
                
                vp   = bodies{ib}.particles.velo(pid,:);
                k = k + 0.5*(vp(1)^2+vp(2)^2+vp(3)^2)*Mp;
                u = u + 0.5*Vp*bodies{ib}.particles.stress(pid,:)*bodies{ib}.particles.strain(pid,:)';
            end
        end
    end
    
    % update the element particle list
    
    for ib=1:length(bodies)
        body      = bodies{ib};
        particles = body.particles;
        elems     = ones(length(particles.volume),1);
        
        for ip=1:length(particles.volume)
            x = particles.coord(ip,1); y = particles.coord(ip,2); z = particles.coord(ip,3);
            e = floor(x/deltax) + 1 + numx*floor(y/deltay) + numx*numy*floor(z/deltaz) ;
            elems(ip) = e;
        end
        
        bodies{ib}.elements = unique(elems);
        bodies{ib}.nodes    = unique(element(bodies{ib}.elements,:));
        
        mpoints = cell(elemCount,1);
        for ie=1:elemCount
            id  = find(elems==ie);
            mpoints{ie}=id;
        end
        
        bodies{ib}.mpoints  = mpoints;
    end
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % VTK output
    
    if (  mod(istep-1,interval) == 0 )
        xp = [bodies{1}.particles.coord;bodies{2}.particles.coord];
        s  = [bodies{1}.particles.stress;bodies{2}.particles.stress];
        vtkFile = sprintf('../results/mpm/two-disks-3D/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,s);
    end
    
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
    
    nvelo0 = nvelo;
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


Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
    [sigmaXX sigmaYY sigmaXY],[Ux Uy]);

% plot kinetic, strain energies

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
