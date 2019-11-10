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

v   = 0.1;

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

vtkFileName  = 'gimp2DTwoDisks';
interval     = 50;

tic;

%% Computational grid

l = 1.1;

numx2 = 22;      % number of elements along X direction
numy2 = 22;      % number of elements along Y direction

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

noParticleX = 2;
noParticleY = 2;

[W,Q]=quadrature(  2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature


% body1

center = [0.2+deltax 0.2+deltay];
radius = 0.2;


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
        r  = norm(x-center);
        if ( r-radius < 0 )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

velo1   = ones(pCount,2)*v;               % velocity
velo2   = -ones(pCount,2)*v;               % velocity

% body 2
center = [0.8+deltax 0.8+deltay];
radius = 0.2;


for e=1:elemCount                          % start of element loop
    sctr = element(e,:);          %  element scatter vector
    pts  = node(sctr,:);
    
    for q=1:size(W,1)                      % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
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

pCount = pCount*2;

Mp  = mass;                           % mass
Vp  = volume;                         % volume
Fp  = repmat([1 0 0 1],pCount,1);     % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = [velo1;velo2];                  % velocity
xp  = coord;
Vp0 = Vp;

lpx = deltax/noParticleX;
lpy = deltay/noParticleY;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(elemCount,1);

for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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
    neighbors      = getNeighbors(e, numx2, numy2);
    neighborNodes  = element(neighbors,:);
    gimpElement{e} = unique(neighborNodes);
end

%% plot mesh, particles

hold on
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',20);
plot([deltax 1.05 1.05 deltax deltax],[deltax deltax 1.05 1.05 deltax],'k-','LineWidth',4)
axis off



ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;
time  = 3.3;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

global node

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = gimpElement{e};    % GIMP (extended) element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            
            % particle mass and momentum to node
            stress = s(pid,:);
            [N,dNdx] = Shape_Function_Global_GMPM(esctr,xp(pid,:),[lpx lpy],deltax,deltay);  
            
            for i=1:length(esctr)
                id    = esctr(i);
                x     = xp(pid,:) - node(id,:);
               % [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
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
    
    nmomentum = nmomentum + niforce*dtime;
        
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = gimpElement{e};
        enode = node(esctr,:);
        mpts  = mpoints{e};
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            [N,dNdx] = Shape_Function_Global_GMPM(esctr,xp(pid,:),[lpx lpy],deltax,deltay);  
            sum(dNdx);
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                x     = xp(pid,:) - node(id,:);
               % [N,dNdx]=getGIMP2D(x,deltax,deltay,lpx,lpy);
          
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
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
        pElems(p) = e;
    end
    
    for e=1:elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    if (  mod(istep-1,interval) == 0 )
        xp = pos{istep};        
        vtkFile = sprintf('../results/gimp/two-disks/%s%d',vtkFileName,istep-1);
        stress = [s sum(s,2)/3];
        VTKParticles(xp,vtkFile,stress);
    end
    
    
    % store time,velocty for plotting
    
    ta = [ta;t];
    ka = [ka;k];
    sa = [sa;u];
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
end


%% post processing

disp([num2str(toc),'   POST-PROCESSING '])


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
