% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
% 
% Modified Update Stress Last (MUSL): double mapping technique.
%
% Two elastic disks come into contact.
% Example taken from Sulsky et al. 1994 paper.
%
% Vinh Phu Nguyen
% Cardiff University, Wales, UK
% June 2013.

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


%% Material properties 
% 

E   = 1000;        % Young's modulus
nu  = 0.3;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v   = 0.1;

interval     = 100;
vtkFileName  = 'mpm2DTwoDisks';

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

l = 1;

numx2 = 20;      % number of elements along X direction
numy2 = 20;      % number of elements along Y direction

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
nvelo     = zeros(nodeCount,2);  % nodal velocity vector (used for MUSL only)
niforce   = zeros(nodeCount,2);  % nodal internal force vector
neforce   = zeros(nodeCount,2);  % nodal external force vector

%% plot mesh, particles

hold on
plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
axis off

ta = [];           % time
ka = [];           % kinetic energy 
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.001;           % time increment
time  = 3.;             % total simulation time 
t     = 0;


nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;

activeNodes = unique(element(unique(pElems),:));

while ( t < time )
    disp(['time step ',num2str(t)])
    %% reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    niforce(:)   = 0;    
    nvelo(:)     = 0;
    %% loop over computational cells or elements
    for e=1:elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords        
        mpts  = mpoints{e};        % particles inside element e 
        
        % loop over particles        
        for p=1:length(mpts)
            pid  = mpts(p);
            % compute shape functions and derivatives for all nodes at
            % particle "pid"
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            
            stress = s(pid,:);
            pmass  = Mp(pid);
            pvelo  = vp(pid,:);
            pvol   = Vp(pid);
            % particle mass and momentum to node
            for i=1:length(esctr)
                id    = esctr(i);
                NI    = N(i);
                dNIdx = dNdx(i,1);
                dNIdy = dNdx(i,2);
                nmass(id)       = nmass(id)       + NI*pmass;
                nmomentum(id,:) = nmomentum(id,:) + NI*pmass*pvelo;
                niforce(id,1)   = niforce(id,1) - pvol*(stress(1)*dNIdx + stress(3)*dNIdy);
                niforce(id,2)   = niforce(id,2) - pvol*(stress(3)*dNIdx + stress(2)*dNIdy);
            end                                                            
        end
    end
    
    % debug
   
    %% update nodal momenta
    
    nmomentum = nmomentum + niforce*dtime;
   
    %% update particle velocity

    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
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
            
            for i=1:length(esctr)
                id = esctr(i);     
                vp(pid,:)  = vp(pid,:) + dtime * N(i)*niforce(id,:)/nmass(id);                                           
            end        
            
            for i=1:length(esctr)
                id = esctr(i);
                nvelo(id,:)  = nvelo(id,:) + N(i)*vp(pid,:)*Mp(pid);
            end         
        end
    end
    
    nvelo(activeNodes,1) = nvelo(activeNodes,1) ./ nmass(activeNodes);
    nvelo(activeNodes,2) = nvelo(activeNodes,2) ./ nmass(activeNodes);
    
    % update particle  position and stresses
    k = 0;
    u = 0;
    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
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
            dEps    = dtime * 0.5 * (Lp+Lp');    
            dsigma  = C * [dEps(1,1);dEps(2,2);2*dEps(1,2)] ;
            s(pid,:)  = s(pid,:) + dsigma';     
            eps(pid,:)= eps(pid,:) + [dEps(1,1) dEps(2,2) 2*dEps(1,2)];     
            
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
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
        pElems(p) = e;
    end
    
    for e=1:elemCount
        id  = find(pElems==e);
        mpoints{e}=id;
    end
    
    activeNodes = unique(element(unique(pElems),:));
    
    % store time,velocty for plotting
    
    ta = [ta;t];   
    ka = [ka;k];
    sa = [sa;u];
  
    % VTK output
    
    if (  mod(istep-1,interval) == 0 )
        xp = pos{istep};        
        vtkFile = sprintf('../results/mpm/two-disks-elastic/%s%d',vtkFileName,istep-1);
        VTKParticles(xp,vtkFile,s);
    end
    
        
    % advance to the next time step
    
    t     = t + dtime;
    istep = istep + 1;
end

    

%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

pvdFile = fopen(strcat('../results/mpm/two-disks-elastic/',vtkFileName,'.pvd'), 'wt');

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
