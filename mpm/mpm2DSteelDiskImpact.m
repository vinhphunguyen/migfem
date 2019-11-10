% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% USF formulation and the contact-release algorithm of York.
% Use data structures suitable for contact bodies.
% The shape functions are written in the global coords so that one does not
% have to convert particle position to natural coordinates.
%
% Steel disk impacts the plastic target.
% Example taken from Coetze's PhD thesis.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% March 2014.

%%

addpath ../fem_util/
addpath ../fem-functions/
addpath ../post-processing/

%%
clc
clear all
colordef white

eye2x2 = [1 1 0; 
          1 1 0; 
          0 0 0];
I_dev  = eye(3) - 0.5*eye2x2;

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties
%

E2   = 200e3;       % Young's modulus
nu2  = 0.3;         % Poisson ratio
rho2 = 7850e-12;    % density
mu2  = E2/2/(1+nu2);    % shear modulus

E1   = 78.2e3;      % Young's modulus
nu1  = 0.3;         % Poisson ratio
rho1 = 2700e-12;    % density
yield=300;          % yield stress
k1   = 0;           % Hardening modulus (perfect plasticity)
mu1      = E1/2/(1+nu1);    % shear modulus
lambda1  = E1*nu1/((1+nu1)*(1-2*nu1));
kappa1   = lambda1 + mu1;

v0   = 1160*1000;   % initial velocity of the disk

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);


vtkFileName  = 'mpm2DSteelDiskImpact';
vtkFileName1 = '../results/mpm2DSteelDiskImpactGrid';
interval     = 100;

tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

l = 60;

numx2 = 50;      % number of elements along X direction
numy2 = 50;      % number of elements along Y direction

deltax = l/numx2;
deltay = l/numy2;

nnx=numx2+1;
nny=numy2+1;
node=square_node_array([0 0],[l 0],[l l],[0 l],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

elemCount = size(element,1);
nodeCount = nnx*nny;

%% generate material points

[W,Q]=quadrature(  2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature


% the aluminium target

h = 40;


volume = [];
mass   = [];
coord  = [];


for e=1:elemCount                          % start of element loop
    sctr = element(e,:);                   % element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                      % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        if ( x(2) < h )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho1];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

body1.volume  = volume;
body1.volume0 = volume;
body1.mass    = mass;
body1.coord   = coord;
body1.deform  = repmat([1 0 0 1],pCount,1);     % gradient deformation
body1.stress  = zeros(pCount,3);                % stress
body1.strain  = zeros(pCount,3);                % strain
body1.pstrain = zeros(pCount,3);                % plastic strain
body1.alpha   = zeros(pCount,1);                % equivalent plastic strain
body1.velo    = zeros(pCount,2);                % velocity


% body 2 (steel disk)
center = [l/2 40+10];
radius = 9.53/2;

volume = [];
mass   = [];
coord  = [];

for e=1:elemCount                          % start of element loop
    sctr = element(e,:);                   %  element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        r  = norm(x-center);
        if ( r-radius < 0 )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho2];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

body2.volume = volume;
body2.volume0 = volume;
body2.mass   = mass;
body2.coord  = coord;
body2.deform = repmat([1 0 0 1],pCount,1);     % gradient deformation
body2.stress = zeros(pCount,3);                % stress
body2.strain = zeros(pCount,3);                % strain
body2.velo   = zeros(pCount,2);                % velocity
body2.velo(:,2) = -v0;

bodies    = cell(2,1);
bodies{1} = body1;
bodies{2} = body2;
bodyCount = length(bodies);

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles


for ib=1:length(bodies)
    body      = bodies{ib};
    coord     = body.coord;
    elems     = ones(size(coord,1),1);
    
    for ip=1:size(coord,1)
        x = coord(ip,1); y = coord(ip,2);
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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

% find boundary conditions

bottom = find(node(:,2)==0);
left   = find(node(:,1)==0);
right  = find(abs(node(:,1)-60)<1e-12);

fixedNodes = [bottom; left; right];

%% plot mesh, particles

coords=[bodies{1}.coord;bodies{2}.coord];
hold on
plot_mesh(node,element,'Q4','b-',0.8);
%plot_mesh(node,element(bodies{1}.elements,:),'Q4','cy-',2.1);
%plot_mesh(node,element(bodies{2}.elements,:),'Q4','r-',2.1);
plot(coords(:,1),coords(:,2),'ks','markersize',3);
axis off

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2); % nodal velocities
nvelo0    = zeros(nodeCount,2*(bodyCount+1));
nacce     = zeros(nodeCount,2);  % nodal acceleration vector

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-16; % mass tolerance

c      = sqrt(E2/rho2);
critdT = deltax/c;
dtime  = critdT/20;
time   = 50e-6;
t      = 0;


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
    nvelo(:)     = 0;
    
    %% loop over bodies (update nodal momenta without contact)
    
    for ib=1:bodyCount
        body      = bodies{ib};
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
                xp     = body.coord(pid,:);
                stress = body.stress(pid,:);
                Mp     = body.mass(pid);
                vp     = body.velo(pid,:);
                Vp     = body.volume(pid);
                % loop over nodes of current element "ie"
                for i=1:length(esctr)
                    id    = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getMPM2D(x,deltax,deltay);
                    
                    dNIdx = dNdx(1);
                    dNIdy = dNdx(2);
                    nmass(id)       = nmass(id)       + N*Mp;
                    nmomentum(id,:) = nmomentum(id,:) + N*Mp*vp;
                    niforce(id,1)   = niforce(id,1) - Vp*(stress(1)*dNIdx + stress(3)*dNIdy);
                    niforce(id,2)   = niforce(id,2) - Vp*(stress(3)*dNIdx + stress(2)*dNIdy);
                end
            end
        end
        
        % update nodal momenta
        
        activeNodes=bodies{ib}.nodes;
        nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + niforce(activeNodes,:)*dtime;
        
        nmomentum(fixedNodes,:) = 0; % fixed boundary conditions
        
%         nvelo(activeNodes,2*ib-1) = nmomentum(activeNodes,1)./nmass(activeNodes);
%         nvelo(activeNodes,2*ib)   = nmomentum(activeNodes,2)./nmass(activeNodes);
    end
    
    %% update particle velocity 
    
    for ib=1:bodyCount
        body      = bodies{ib};
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
                xp   = body.coord(pid,:);
                Mp   = body.mass(pid);
                vp   = body.velo(pid,:);
                Vp   = body.volume(pid);
   
                for i=1:length(esctr)
                    id = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getMPM2D(x,deltax,deltay);
                    vp  = vp + dtime * N*niforce(id,:)/nmass(id);
                end
                
                bodies{ib}.velo(pid,:) = vp;
                
                for i=1:length(esctr)
                    id = esctr(i);
                    x     = xp - node(id,:);
                    [N,dNdx]=getMPM2D(x,deltax,deltay);
                    nvelo(id,1:2)  = nvelo(id,1:2) + N*Mp*vp;
                end   
            end
        end
    end
    
    activeNodes = [bodies{1}.nodes; bodies{2}.nodes];
    
    nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
    nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
    
    nvelo(fixedNodes,:) = 0;
    
    k = 0; u = 0;
    for ib=1:bodyCount
        body      = bodies{ib};
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
                xp   = body.coord(pid,:);
                xp0  = body.coord(pid,:);
                Mp   = body.mass(pid);
                Vp   = body.volume(pid);
                
                Lp   = zeros(2,2);
                for i=1:length(esctr)
                    id = esctr(i);
                    x     = xp0 - node(id,:);
                    [N,dNdx]=getMPM2D(x,deltax,deltay);
                    vI  = nvelo(id,1:2);
                    xp  = xp + dtime * N*nmomentum(id,:)/nmass(id);   
                    Lp  = Lp  + vI'*dNdx;         % particle gradient velocity
                end
                
                
                bodies{ib}.coord(pid,:)= xp;
                F       = ([1 0;0 1] + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
                bodies{ib}.deform(pid,:) = reshape(F,1,4);
                bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                if ( ib == 2 ) % elastic
                    dsigma  = C2 * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
                    bodies{ib}.stress(pid,:)  = bodies{ib}.stress(pid,:) + dsigma';
                else           % plastic
                     
                    epsilon   =  bodies{ib}.strain (pid,:)';   % column vector
                    epsilonp0 =  bodies{ib}.pstrain(pid,:)';
                    alpha0    =  bodies{ib}.alpha  (pid);
                    
                    [sigma,alpha,epsilonp] = updateVonMisesMaterial ( epsilon,epsilonp0,alpha0,mu1,kappa1,k1,yield );
                    
                    bodies{ib}.stress (pid,:) = sigma;
                    bodies{ib}.pstrain(pid,:) = epsilonp;
                    bodies{ib}.alpha  (pid,:) = alpha;
                end
                
                
                vp   = bodies{ib}.velo(pid,:);
                k = k + 0.5*(vp(1)^2+vp(2)^2)*Mp;
                u = u + 0.5*Vp*bodies{ib}.stress(pid,:)*bodies{ib}.strain(pid,:)';
            end
        end
    end
    
    % update the element particle list
    
    for ib=1:length(bodies)
        body      = bodies{ib};
        coord     = body.coord;
        elems     = ones(size(coord,1),1);
        
        for ip=1:length(elems)
            x = coord(ip,1); y = coord(ip,2);
            e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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
        xp = [bodies{1}.coord;bodies{2}.coord];
        s2 = [bodies{2}.stress zeros(size(bodies{2}.stress,1),1)];
        s1 = [bodies{1}.stress zeros(size(bodies{1}.stress,1),1)];
        for i=1:size(bodies{1}.stress,1)
            devSig = I_dev*s1(i,1:3)';
            s1(i,4)= sqrt(  devSig(1)^2 + devSig(2)^2 + 2*devSig(3)^2  );
        end
        data.stress  = [s1;s2];
        data.pstrain = [bodies{1}.alpha;zeros(size(bodies{2}.stress,1),1)];
        vtkFile = sprintf('../results/mpm/steel-disk/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
        pos{istep} = xp;
    end
    
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
    
    nvelo0 = nvelo;
end


%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

figure
coords=pos{1297};
hold on
plot_mesh(node,element,'Q4','b-',1.);
plot(coords(:,1),coords(:,2),'ks','markersize',3);
axis off


figure
hold on
plot_mesh(node,element,'Q4','b-',1.);
scatter(coords(:,1),coords(:,2),50,[bodies{1}.particles.stress(:,2);bodies{2}.particles.stress(:,2)],'filled');
colorbar
axis off



pvdFile = fopen(strcat('../results/mpm/steel-disk/',vtkFileName,'.pvd'), 'wt');

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
