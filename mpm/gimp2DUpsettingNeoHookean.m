% This file implements the Generalized Interpolation Material Point Method
% (GIMP). Only uGIMP is implemented where the particle domain is not updated.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
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
addpath ../post-processing/


%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


I  = [1 0;0 1];

%% Material properties
%

% platen
E2   = 200e3;       % Young's modulus
nu2  = 0.3;         % Poisson ratio
rho2 = 7850e-12;    % density
mu2  = E2/2/(1+nu2);    % shear modulus

% billet

E1   = 200e3;      % Young's modulus
nu1  = 0.3;         % Poisson ratio
rho1 = 7800e-12;    % density
yield=700;          % yield stress
k1   = 300;           % Hardening modulus
mu1      = E1/2/(1+nu1);    % shear modulus
lambda1  = E1*nu1/((1+nu1)*(1-2*nu1));
kappa1   = lambda1 + mu1;

v0   = 1*1000;   %velocity of the platen

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C1 = elasticityMatrix(E1,nu1,stressState);
C2 = elasticityMatrix(E2,nu2,stressState);


vtkFileName  = 'mpm2DUpSetting';
interval     = 1;


tic;

%% Computational grid

w0    = 15;
numy0 = 15;

l = 30;
w = w0 + 2*w0/numy0;

numx2 = 30;      % number of elements along X direction
numy2 = numy0+2; % number of elements along Y direction

deltax = l/numx2;
deltay = w/numy2;

nnx=numx2+1;
nny=numy2+1;
node=square_node_array([0 0],[l 0],[l w],[0 w],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

elemCount = size(element,1);
nodeCount = nnx*nny;

%% generate material points

noParticleX = 2;
noParticleY = 2;

[W,Q]=quadrature(  4, 'GAUSS', 2 ); % 2x2 Gaussian quadrature


% cylinder

hx = 10;
hy = 15;


volume = [];
mass   = [];
coord  = [];


for e=1:elemCount                 % start of element loop
    sctr = element(e,:);          %  element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        
        if ( x(1) < hx ) && ( x(2) < hy )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho1];
            coord   = [coord;x];
        end
    end
end

coord(:,1) = coord(:,1) + deltax;
coord(:,2) = coord(:,2) + deltay;

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

% platen

hx = 20;
hy = 15+2*deltay;

volume = [];
mass   = [];
coord  = [];

for e=1:elemCount                          % start of element loop
    sctr = element(e,:);                   %  element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                      % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        x(1) = x(1) + deltax;
        x(2) = x(2) + deltay;
        if ( x(2) < hy ) && ( x(2) > 15+deltay ) && ( x(1) < l-deltax )
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

%% node quantities

nmass     = zeros(nodeCount,1);  % nodal mass vector
nmomentum = zeros(nodeCount,2);  % nodal momentum vector
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector
nvelo     = zeros(nodeCount,2); % nodal velocities

%% find boundary conditions

eps=1e-12;
bottom = find(abs(node(:,2)-deltay)<eps);
left   = find(abs(node(:,1)-deltax)<eps);

fixedBoth = bottom;
fixedX    = left;

%% plot mesh, particles

figure
hold on
plot_mesh(node,element,'Q4','k-',1.);
xp1 = bodies{1}.coord;
xp2 = bodies{2}.coord;
plot(xp1(:,1),xp1(:,2),'k.','markersize',15);
plot(xp2(:,1),xp2(:,2),'r.','markersize',15);
plot(node(bottom,1),node(bottom,2),'cy.','markersize',15);
plot(node(left,1),node(left,2),'b.','markersize',15);
%axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E1/rho1);
dtime = deltax/c*0.07;
time  = 7e-3;
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
    nvelo(:)     = 0;
    
    %% loop over bodies (update nodal momenta without contact)
    
    for ib=1:1
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
    end
    
    
    nmomentum(fixedBoth,2) = 0; % fixed boundary conditions
    nmomentum(fixedX,1)    = 0;
    
    niforce(fixedBoth,2) = 0; % fixed boundary conditions
    niforce(fixedX,1)    = 0;
    
    nmomentum(bodies{2}.nodes,2) = nmass(bodies{2}.nodes)*(-v0);
    nmomentum(bodies{2}.nodes,1) = 0;
    
    %% update particle velocity
    
    for ib=1:1
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
    
    %     nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
    %     nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
    
    nvelo(fixedBoth,2) = 0;
    nvelo(fixedX,   1) = 0;
    
    nvelo(bodies{2}.nodes,2) = -v0;
    nvelo(bodies{2}.nodes,1) = 0;
    
    k = 0; u = 0;
    for ib=1:1
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
                    xp  = xp + dtime * vI*N;
                    Lp  = Lp  + vI'*dNdx;         % particle gradient velocity
                end
                
                
                bodies{ib}.coord(pid,:)= xp;
                F       = (I + Lp*dtime)*reshape(bodies{ib}.deform(pid,:),2,2);
                bodies{ib}.deform(pid,:) = reshape(F,1,4);
                bodies{ib}.volume(pid  ) = det(F)*bodies{ib}.volume0(pid);
                dEps    = dtime * 0.5 * (Lp+Lp');
                bodies{ib}.strain(pid,:)  = bodies{ib}.strain(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];
                
                
                
                epsilon   =  bodies{ib}.strain (pid,:)';   % column vector
                epsilonp0 =  bodies{ib}.pstrain(pid,:)';
                alpha0    =  bodies{ib}.alpha  (pid);
                
                [sigma,alpha,epsilonp] = updateVonMisesMaterial ( epsilon,epsilonp0,alpha0,mu1,kappa1,k1,yield );
                
                bodies{ib}.stress (pid,:) = sigma;
                bodies{ib}.pstrain(pid,:) = epsilonp;
                bodies{ib}.alpha  (pid,:) = alpha;
                
                
                vp   = bodies{ib}.velo(pid,:);
                k = k + 0.5*(vp(1)^2+vp(2)^2)*Mp;
                u = u + 0.5*Vp*bodies{ib}.stress(pid,:)*bodies{ib}.strain(pid,:)';
            end
        end
    end
    
    bodies{2}.coord = bodies{2}.coord + dtime* bodies{2}.velo;
    
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
        data.stress  = [s1;s2];
        data.pstrain = [bodies{1}.alpha;zeros(size(bodies{2}.stress,1),1)];
        vtkFile = sprintf('../results/mpm/upsetting/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,data);
        pos{istep} = xp;
    end
    
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
    
    nvelo0 = nvelo;
end

disp([num2str(toc),'   DONE '])
