% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
%
% Elastic ball impacts on rigid walls.
% Example taken from Sulsky et al. 1995 paper.
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


%% Material properties
%

E1   = 1000;            % Young's modulus of the left ball
E2   = 1000;            % Young's modulus of the right ball
nu   = 0.3;             % Poisson ratio
rho1 = 1000;            % density
rho2 = 1000;            % density

% kappa = 3-4*nu;        % Kolosov constant
% mu    = E/2/(1+nu);    % shear modulus

v0   = 0.1;              % initial particle velocity

stressState ='PLANE_STRAIN'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C1       = elasticityMatrix(E1,nu,stressState);
C2       = elasticityMatrix(E2,nu,stressState);
%D           = inv(C);

vtkFileName  = 'mpm2DHollowDisks';
vtkFileName1  = '../results/mpm2DHollowDisksGrid';


tic;

disp([num2str(toc),'   INITIALISATION '])

%% Computational grid

l0 = 260;
l      = l0;
w      = 100;

numx2 = 60;      % number of elements along X direction
numy2 = 30;      % number of elements along Y direction

deltax = l0/numx2;
deltay = w /numy2;

nnx=numx2+1;
nny=numy2+1;
node=square_node_array([0 0],[l 0],[l w],[0 w],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

elemCount = size(element,1);
nodeCount = nnx*nny;

xMax = max(node(:,1));
yMax = max(node(:,2));

%%   particle distribution from a mesh
%

meshFile = 'disk.msh';
mesh     = load_gmsh (meshFile);

elemType = 'T3';
numnode  = mesh.nbNod;
numelem  = mesh.nbTriangles;
node1    = mesh.POS(:,1:2);
element1 = mesh.TRIANGLES(1:numelem,1:3);

% move the disk to correct position

node1(:,1) = node1(:,1) + l/3;
node1(:,2) = node1(:,2) + w/2;

% check if Jacobian is negative

element1  = tricheck(node1,element1,1);

% "ngp" particles per triangle

ngp = 6;
[W,Q]=quadrature( ngp, 'TRIANGULAR', 2 ); % two point quadrature

pCount = 2* numelem*size(W,1);

Mp  = ones(pCount,1);                 % mass
Vp  = ones(pCount,1);                 % volume
Fp  = ones(pCount,4);                 % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position

% initialise particle position, mass ,volume, velocity
% particle taken as global Gauss points

%% ball particles
id=1;
for e = 1:numelem
    sctr = element1(e,:);
    pts  = node1(sctr,:);
    for q=1:size(W)
        pt=Q(q,:);
        wt=W(q);
        [N,dNdxi]=lagrange_basis('T3',pt);  % element shape functions
        J0=dNdxi'*pts;
        detJ0=det(J0);
        a     = wt*detJ0;
        Vp(id) = a;
        Mp(id) = a*rho1;
        xp(id,:) =  N'*pts;
        vp(id,1) = v0;
        Fp(id,:) = [1 0 0 1];
        
        id = id + 1;
    end
end

noDiskParticles = id-1;

for q=1:pCount/2
        Vp(id) = Vp(q);
        Mp(id) = Mp(q)*rho2/rho1;
        xp(id,1) =  xp(q,1) + l/3;
        xp(id,2) =  xp(q,2);
        vp(id,1) = -v0;
        Fp(id,:) = [1 0 0 1];
        
        id = id + 1;
    end

Vp0 = Vp;

% boundary conditions
% eps0=1e-9;
% leftNode  = find( abs(node(:,1)-deltax)     < eps0 );
% rightNode = find( abs(node(:,1)-(l-deltax)) < eps0 );
% 
% fixedNodes = [leftNode;rightNode];

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

%% plot mesh, particles

figure(1)
hold on
plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.2);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
%plot(node(fixedNodes,1),node(fixedNodes,2),'r.','markersize',20);
axis off

ta = [];           % time
ka = [];           % kinetic energy
sa = [];           % strain energy

%% Solver

disp([num2str(toc),'   SOLVING '])

tol   = 1e-6; % mass tolerance

dtime = 0.1;
time  = 120;
t     = 0;
interval     = 10;

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
    % loop over computational cells or elements
    for e=1:elemCount
        esctr = element(e,:);      % element connectivity
        enode = node(esctr,:);     % element node coords
        mpts  = mpoints{e};        % particles inside element e
        
        if isempty(mpts) continue; end % skip element without particles
        
        % loop over particles
        for p=1:length(mpts)
            pid  = mpts(p);
            pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
            pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));
            
            [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
            J0       = enode'*dNdxi;             % element Jacobian matrix
            invJ0    = inv(J0);
            dNdx     = dNdxi*invJ0;
            
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
    
    if ~isempty( find(nmass<0))
        disp('NEGATIVE NODAL MASS')
    end
    
%     nmomentum(fixedNodes,:) = 0;
%     niforce  (fixedNodes,:) = 0;
    
    % update nodal momenta
%     for i=1:nodeCount
%         nmomentum(i,:) = nmomentum(i,:) + niforce(i,:)*dtime;
%     end
    
    nmomentum = nmomentum + niforce*dtime;
   % nmomentum(fixedNodes,:) = 0;
    
    % update particle velocity and position and stresses
    k = 0;
    u = 0;
    
    for e=1:elemCount
        esctr = element(e,:);
        enode = node(esctr,:);
        mpts  = mpoints{e};
        
        if isempty(mpts) continue; end % skip element without particles
         
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
                vI = [0 0];
                % if nmass(id) > tol
                    nmassIdInv = 1/nmass(id);     
                    vp(pid,:)  = vp(pid,:) + dtime * N(i)*niforce  (id,:) * nmassIdInv;
                    xp(pid,:)  = xp(pid,:) + dtime * N(i)*nmomentum(id,:) * nmassIdInv;
                    vI         = nmomentum(id,:)*nmassIdInv;  % nodal velocity
                %end
                Lp = Lp + vI'*dNdx(i,:);         % particle gradient velocity
            end
            
            if (pid <= noDiskParticles) 
                C = C1;
            else
                C = C2;
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
        % debug
        if ( x > xMax ) disp('OUT OF CELL'); end
        e = floor(x/deltax) + 1 + numx2*floor(y/deltay);
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
        xp = pos{istep};
        vtkFile = sprintf('../results/%s%d',vtkFileName,istep);
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

Ux= zeros(size(node,1),1);
Uy= zeros(size(node,1),1);
sigmaXX = zeros(size(node,1),1);
sigmaYY = zeros(size(node,1),1);
sigmaXY = zeros(size(node,1),1);

VTKPostProcess(node,element,2,'Quad4',vtkFileName1,...
             [sigmaXX sigmaYY sigmaXY],[Ux Uy]);


figure(2)
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
