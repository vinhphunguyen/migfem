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

rho = 1000;
g   = 9.81;
K_f     = 2e6;
lambda  = 0.001;
gamma   = 7;

CFL = 0.8;


vtkFileName  = 'dam-break';
interval     = 50;

tic;

%% Computational grid
l  = 0.057;
Lx = 5*l;
Ly = 3*l;

noX0 = 70;      % number of elements along X direction
noY0 = 30;      % number of elements along Y direction

ghostCell = 0;

[mesh]=buildGrid2D(Lx,Ly,noX0,noY0, ghostCell);

%% generate material points

noParticle = 2;

[W,Q]=quadrature( noParticle, 'GAUSS', 2 ); % 2x2 Gaussian quadrature

volume = [];
mass   = [];
coord  = [];


for e=1:mesh.elemCount                 % start of element loop
    sctr = mesh.element(e,:);          %  element scatter vector
    pts  = mesh.node(sctr,:);
    
    for q=1:size(W,1)                           % quadrature loop
        pt=Q(q,:);                              % quadrature point
        wt=W(q);                                % quadrature weight
        [N,dNdxi]=lagrange_basis('Q4',pt);
        J0 = pts'*dNdxi;
        x  = N'*pts;
        if ( x(1) < l  ) && ( x(2) < 2*l )
            volume  = [volume;wt*det(J0)];
            mass    = [mass; wt*det(J0)*rho];
            coord   = [coord;x];
        end
    end
end

pCount = length(volume);

Mp  = mass;                           % mass
Vp  = volume;                         % volume
Fp  = repmat([1 0 0 1],pCount,1);     % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = coord;
rhop  = zeros(pCount,1);                % density
IntE  = zeros(pCount,1);                % density
Vp0 = Vp;
rhop(:) = rho;

%% find elements to which particles belong to
% two data structures are used
% 1. particle -> element
% 2. element  -> particles

pElems  = ones(pCount,1);
mpoints = cell(mesh.elemCount,1);

for p=1:pCount
    x = xp(p,1);
    y = xp(p,2);
    e = floor(x/mesh.deltax) + 1 + noX0*floor(y/mesh.deltay);
    pElems(p) = e;
end

for e=1:mesh.elemCount
    id  = find(pElems==e);
    mpoints{e}=id;
end

activeElems = unique(pElems);
activeNodes = unique(mesh.element(activeElems,:));

%% node quantities

nmass     = zeros(mesh.nodeCount,1);  % nodal mass vector
nmomentum = zeros(mesh.nodeCount,2);  % nodal momentum vector
niforce   = zeros(mesh.nodeCount,2);  % nodal internal force vector
neforce   = zeros(mesh.nodeCount,2);  % nodal external force vector
nvelo     = zeros(mesh.nodeCount,2);  % nodal external force vector

%% plot mesh, particles

hold on
plot_mesh(mesh.node,mesh.element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',20);
plot(mesh.node(mesh.lNodes,1),mesh.node(mesh.lNodes,2),'r*','markersize',20);
plot(mesh.node(mesh.rNodes,1),mesh.node(mesh.rNodes,2),'r*','markersize',20);
plot(mesh.node(mesh.bNodes,1),mesh.node(mesh.bNodes,2),'r*','markersize',20);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-12; % mass tolerance

c     = sqrt(K_f/rho);
maxV  = mesh.deltax/c;
dtime = CFL*maxV;
time  = 0.5;
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

while ( t < time )
    disp(['time step ',num2str(t)])
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;
    neforce(:)   = 0;
    nvelo(:)     = 0;
    % loop over computational cells or elements
    for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = mesh.element(e,:);
        enode = mesh.node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            % particle mass and momentum to node
            stress = s(pid,:);
            pt(1)       = (2*xp(pid,1)-(enode(1,1)+enode(2,1)))*mesh.dxInv;
            pt(2)       = (2*xp(pid,2)-(enode(2,2)+enode(3,2)))*mesh.dyInv;
            
            [N,dNdxi]   = lagrange_basis('Q4',pt);	% element shape functions
            J0          = enode'*dNdxi;              % element Jacobian matrix
            dNdx        = dNdxi*inv(J0);
            for i=1:length(esctr)
                id    = esctr(i);
                dNIdx = dNdx(i,1);
                dNIdy = dNdx(i,2);
                nmass(id)       = nmass(id)       + N(i)*Mp(pid);
                nmomentum(id,:) = nmomentum(id,:) + N(i)*Mp(pid)*vp(pid,:);
                niforce(id,1)   = niforce(id,1) - Vp(pid)*(stress(1)*dNIdx + stress(3)*dNIdy);
                niforce(id,2)   = niforce(id,2) - Vp(pid)*(stress(3)*dNIdx + stress(2)*dNIdy);
                neforce(id,2)  =  neforce(id,2) - N(i)*Mp(pid)*g;
            end
        end
    end
    
    % update nodal momenta
    nforce = niforce + neforce;
    nmomentum(activeNodes,:) = nmomentum(activeNodes,:) + nforce(activeNodes,:)*dtime;
    
    % Dirichlet boundary conditions
    nmomentum(mesh.lNodes,1) = 0; nforce(mesh.lNodes,1) = 0;
    nmomentum(mesh.rNodes,1) = 0; nforce(mesh.rNodes,1) = 0;
    nmomentum(mesh.bNodes,2) = 0; nforce(mesh.bNodes,2) = 0;
    
    %% update particle velocity
    
    for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = mesh.element(e,:);
        enode = mesh.node(esctr,:);
        mpts  = mpoints{e};
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))*mesh.dxInv;
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))*mesh.dyInv;
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            
            for i=1:length(esctr)
                id = esctr(i);
                vp(pid,:)  = vp(pid,:) + dtime * N(i)*nforce(id,:)/nmass(id);
            end
            
            for i=1:length(esctr)
                id = esctr(i);
                nvelo(id,:)  = nvelo(id,:) + N(i)*vp(pid,:)*Mp(pid);
            end
        end
    end
    
    nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
    nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
    
    nvelo(mesh.lNodes,1) = 0;
    nvelo(mesh.rNodes,1) = 0;
    nvelo(mesh.bNodes,2) = 0;
    
    % update particle  position and stresses
    k = 0;
    u = 0;
     for ie=1:length(activeElems)
        e     = activeElems(ie);
        esctr = mesh.element(e,:);
        enode = mesh.node(esctr,:);
        mpts  = mpoints{e};
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);
                vI = nvelo(id,:);
                xp(pid,:)  = xp(pid,:) + dtime * N(i)*nmomentum(id,:)/nmass(id);
                Lp = Lp + vI'*dNdx(i,:);         % particle gradient velocity
            end
            
            F       = ([1 0;0 1] + Lp*dtime)*reshape(Fp(pid,:),2,2);
            Fp(pid,:)= reshape(F,1,4);
            Vp(pid) = det(F)*Vp0(pid);
            deps    = dtime * 0.5 * (Lp+Lp');
            eps(pid,:)= eps(pid,:) + [deps(1,1) deps(2,2) 2*deps(1,2)];
            
            temp        = deps(1,1) + deps(2,2);
            rhop(pid)   = rhop(pid) / (1+temp);
            IntE(pid)   = IntE(pid) + s(pid,:)*[deps(1,1); deps(2,2); 2*deps(1,2)]/rhop(pid);
            p           = (gamma-1)*rhop(pid)*IntE(pid);
            s(pid,:) = -p*[1; 1; 0] + 2 * lambda * [deps(1,1); deps(2,2); 2*deps(1,2)] ...
                - 2/3 * lambda * temp*[1; 1; 0];
            
            k = k + 0.5*(vp(pid,1)^2+vp(pid,2)^2)*Mp(pid);
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
        e = floor(x/mesh.deltax) + 1 + noX0*floor(y/mesh.deltay);
        pElems(p) = e;
    end
    
    for e=1:mesh.elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    activeElems = unique(pElems);
    activeNodes = unique(mesh.element(activeElems,:));
    
    if (  mod(istep-1,interval) == 0 )
        xp = pos{istep};
        vtkFile = sprintf('../results/mpm/dam-break/%s%d',vtkFileName,istep-1);
        
        stress = [s sum(s,2)/3];
        data.stress = stress; data.pstrain=[];data.velo=vp;
        VTKParticles(xp,vtkFile,data);
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

hold on
plot_mesh(mesh.node,mesh.element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',20);
axis off

disp([num2str(toc),'   DONE '])
