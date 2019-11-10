% This file implements the Material Point Method of Sulsky 1994.
% Two dimensional problems.
% The grid is a structured mesh consisting of 4-noded bilinear elements (Q4).
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

%%
clc
clear all
colordef white

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,'splinecurve.eps',opts)


%% Material properties 
% 

E   = 1e5;        % Young's modulus
nu  = 0.0;         % Poisson ratio
rho = 1000;        % density
kappa = 3-4*nu;    % Kolosov constant
mu    = E/2/(1+nu);% shear modulus

v   = 0.1;

stressState ='PLANE_STRESS'; % either 'PLANE_STRAIN' or "PLANE_STRESS
C = elasticityMatrix(E,nu,stressState);
D = inv(C);

tic;

disp([num2str(toc),'   INITIALISATION '])
%
%%   particle distribution from a mesh
%
%{
meshFile = 'disks.msh';
mesh     = load_gmsh (meshFile);
elemType = 'T3';
numnode  = mesh.nbNod;
numelem  = mesh.nbTriangles;
node1    = mesh.POS(:,1:2);
element1 = mesh.TRIANGLES(1:numelem,1:3);
%}

% Generating rectangular particles
Dmp = 1;    % Dimensions of the solid body 
Lmp = 5;
npx = 2;    % Number of particle in x direction
npy = 10;   % Number of particle in y direction

% Four corner of the material points
pt1 = [0 0] ;
pt2 = [Dmp 0] ;
pt3 = [Dmp Lmp] ;
pt4 = [0 Lmp] ;
nnx = npx + 1;
nny = npy + 1;
elemType = 'Q4' ;
[node1,element1] = meshRectangularRegion(...
    pt1, pt2, pt3, pt4, nnx,nny,elemType);
numnode = size(node1,1);
numelem = size(element1,1);

% check if Jacobian is negative

%element1  = tricheck(node1,element1,1);

% one particle per triangle

pCount = numelem;

Mp  = ones(pCount,1);                 % mass
Vp  = ones(pCount,1);                 % volume
Fp  = ones(pCount,4);                 % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = zeros(pCount,2);                % velocity
xp  = zeros(pCount,2);                % position
dimp = zeros(pCount,2);               % X and Y dimensions of particles
dimp(:,1) = Dmp/(npx);
dimp(:,2) = Lmp/(npy);
% initialise particle position, mass ,volume, velocity

for e = 1:numelem
    coord = node1(element1(e,:),:);
    a     = (coord(2,1)-coord(1,1))*(coord(3,2)-coord(1,2));
    Vp(e) = a;
    Mp(e) = a*rho;
    xp(e,:) = mean(coord);
    %{
    if xp(e,1) < 0.5
        vp(e,:) = [v v];
    else
        vp(e,:) = [-v -v];
    end
    %}
    Fp(e,:) = [1 0 0 1];
end
Vp0 = Vp;
xp0 = xp;
dimp0 = dimp;
%}
%% Computational grid

D = 1;
L = 5; 

numx2 = 1;      % number of elements along X direction
numy2 = 5;      % number of elements along Y direction

deltax = D/numx2;
deltay = L/numy2;

nnx=numx2+1;
nny=numy2+1;
global node element
node=square_node_array([0 0],[D 0],[D L],[0 L],nnx,nny);
inc_u=1;
inc_v=nnx;
node_pattern=[ 1 2 nnx+2 nnx+1 ];
element =make_elem(node_pattern,numx2,numy2,inc_u,inc_v);

elemCount = size(element,1);
nodeCount = nnx*nny;

% Boundary nodes
% define essential boundaries
uln = nnx*(nny-1)+1;       % upper left node number
urn = nnx*nny;             % upper right node number
lrn = nnx;                 % lower right node number
lln = 1;                   % lower left node number
cln = nnx*(nny-1)/2+1;     % node number at (0,0)

topEdge  = [ uln:1:(urn-1); (uln+1):1:urn ]';
botEdge  = [ lln:1:(lrn-1); (lln+1):1:lrn ]';
rightedge= [nnx:nnx:(nnx*(nny-1));(2*nnx):nnx:(nnx*nny)]';
leftedge = [1:nnx:uln-nnx ; (nnx+1):nnx:uln]' ; 
leftedge1 = [136 151 167; 151 167 182]';
% GET NODES ON DIRICHLET BOUNDARY AND ESSENTIAL BOUNDARY
botNodes   = unique(botEdge);
topNodes   = unique(topEdge);
rightNodes = unique(rightedge);
leftNodes  = unique(leftedge) ;

ydispNodes = [botNodes];
xdispNodes = [leftNodes;rightNodes];

% these two matrix are used to find the neighboring elements of the
% particle's element in GIMPM
global elementpos elementmatrix
for i=1:size(element,1)
    row = floor (i/(nnx-1)) ;
    p = i-(row*(nnx-1)); 
    if p==0
       row = row;
    else
       row = row+1; 
    end
    col = (i-((row-1)*(nnx-1))) ;
    elementpos(i,:) = [row,col];
    elementmatrix(row,col) = i ;
end
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
ntforce   = zeros(nodeCount,2);  % nodal total force vector

%% plot mesh, particles
figure
hold on
%plot_mesh(node1,element1,elemType,'r-',0.2);
plot_mesh(node,element,'Q4','k-',1.);
%plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
n3 = plot (node(xdispNodes,1),node(xdispNodes,2),'bo');
n4 = plot (node(ydispNodes,1),node(ydispNodes,2),'b^');
axis off

%% Solver

shapefun = 'GMPM';  % MPM or GMPM
disp(['Shape Function: ',  shapefun]);
disp([num2str(toc),'   SOLVING '])

tol   = 1e-8; % mass tolerance

c     = sqrt(E/rho);
dtime = 0.01;
time  = 0.01;
t     = 0;
grav = [0 -2];   % Gravity Force
damping = 0.9;   % Damping ratio

nsteps = floor(time/dtime);

pos   = cell(nsteps,1);
vel   = cell(nsteps,1);
istep = 1;
ta = [];           % time
ka = [];           % kinetic energy 
sa = [];           % strain energy
while ( t < time )
    % reset grid data
    nmass(:)     = 0;
    nmomentum(:) = 0;
    ntforce(:)   = 0;
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

            switch shapefun 
                case 'MPM'
                pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
                pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));

                [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
                J0       = enode'*dNdxi;             % element Jacobian matrix
                invJ0    = inv(J0);
                dNdx     = dNdxi*invJ0;
                case 'GMPM'
                % All the surrounding elements' nodes are considered 
                [nodes] = assembly_GMPM(e,nnx,nny);
                esctr = nodes; % Nodes which might be involved with this particle 
                [dNdx,N] = Shape_Function_Global_GMPM(nodes,xp(pid,:),dimp(pid,:),deltax,deltay,shapefun);  
            end
            
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
                neforce(id,2)   = neforce(id,2) + N(i)*Mp(pid)*grav(2);
            end                                                            
        end
    end
    
    % Local damping 
    ntforce = neforce + niforce;
    ntforce = ntforce - damping*abs(ntforce).*sign(nmomentum);
    
    % update nodal momenta
    for i=1:nodeCount
        nmomentum(i,:) = nmomentum(i,:) + ntforce(i,:)*dtime;
    end
    
    % boundary condition
    nmomentum(ydispNodes,2) = 0; nmomentum(xdispNodes,1) = 0;
    ntforce(ydispNodes,2) = 0;   ntforce(xdispNodes,1) = 0;
    
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
            
            switch shapefun 
                case 'MPM'
                pt(1)= (2*xp(pid,1)-(enode(1,1)+enode(2,1)))/(enode(2,1)-enode(1,1));
                pt(2)= (2*xp(pid,2)-(enode(2,2)+enode(3,2)))/(enode(3,2)-enode(2,2));

                [N,dNdxi]=lagrange_basis('Q4',pt);   % element shape functions
                J0       = enode'*dNdxi;             % element Jacobian matrix
                invJ0    = inv(J0);
                dNdx     = dNdxi*invJ0;
                case 'GMPM'
                [nodes] = assembly_GMPM(e,nnx,nny);
                esctr = nodes; % Nodes which might be involved with this particle 
                [dNdx,N] = Shape_Function_Global_GMPM(nodes,xp(pid,:),dimp(pid,:),deltax,deltay,shapefun);   
            end
            
            Lp = zeros(2,2);
            for i=1:length(esctr)
                id = esctr(i);     
                vI = [0 0];
                if nmass(id) > tol
                    vp(pid,:)  = vp(pid,:) + dtime * N(i)*ntforce(id,:)  /nmass(id);
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
            
            % Updating the dimensions of particles
            %dimp(pid,:) = ( [F(1,1),F(2,2)].*dimp0(pid,:) );
            dimp(pid,:) = dimp(pid,:) + ([dEps(1,1),dEps(2,2)] .* dimp(pid,:));
            
            
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
    
    % store time,velocty for plotting
    
    ta = [ta;t];   
    ka = [ka;k];
    sa = [sa;u];
    
    % advance to the next time step
    
    t = t + dtime;
    istep = istep + 1;
end
a=0;
    for e=1:elemCount        
        pts=mpoints{e};
        a = a + length(pts);
    end
    
%% post processing

disp([num2str(toc),'   POST-PROCESSING '])

%%% Analytical solution for settlement of an elastic column
%
hold on
scatter(xp(:,1),xp(:,2),100,s(:,2),'fill')
colorbar

counter = 0;
for i=0:.2:Lmp
   counter = counter + 1 ;
    sigm_analy(counter,1) = (Lmp-i)*grav(2)*rho;
    sigm_analy(counter,2) = i;
    
    deform_analy(counter,1) = (rho*grav(2)/E)*(Lmp*i-((i^2)/2));
    deform_analy(counter,2) = i;
end 

% Error
for i=1:pCount
    def = ((rho*grav(2)/E)*(Lmp*xp0(i,2)-((xp0(i,2)^2)/2)));
   deform_error(i,1) = 100*(def - (xp(i,2)-xp0(i,2))) / def; 
   deform_error(i,2) = xp0(i,2);
end

figure
hold on
plot(sigm_analy(:,1),sigm_analy(:,2))
hold on 
plot(s(:,2),xp0(:,2),'r*');
xlabel('Stress (Pa)');
ylabel('Position');
title( sprintf( '2D %s Verical Stress (--Analytical   *Numerical)', shapefun ) );

figure
hold on
plot(deform_analy(:,1),deform_analy(:,2))
hold on 
plot((xp(:,2)-xp0(:,2)),xp0(:,2),'r*');
xlabel('Displacement (m)');
ylabel('Position');
title( sprintf( '2D %s Deformation (--Analytical   *Numerical)', shapefun ) );

figure
plot(deform_error(:,1),deform_error(:,2),'ob')
xlabel('Error (%)');
ylabel('Position');
title( sprintf( '2D %s Deformation Error', shapefun ) );

%
% interval = 25;
% m = 1;
% clf;
% for it=1:nsteps   
%     if (  mod(it,interval) == 0 )
%     xp = pos{it};
%     vp = vel{it};
%     figure(m)
%     hold on
%     plot_mesh(node,element,'Q4','k-',1.1);
%     %scatter(xp(:,1),xp(:,2),40,sqrt(vp(:,1).^2+vp(:,2).^2),'full');
%     plot(xp(:,1),xp(:,2),'k.','markersize',10);
%     %plot_mesh(node,element(pElems,:),'Q4','cy-',2.1);
%     axis equal
%     axis([0 1 0 1])
%     title( sprintf( '%s: %f', 'time', ta(it) ) );
%     colorbar
%     print(m, '-djpeg90', ['disks-',num2str(m),'.jpg']);
%     exportfig(gcf,sprintf( '%s %d %s', 'disks', m, '.eps' ),opts)
%     m = m + 1;
%     end
% end


% figure 
% set(gca,'FontSize',14)
% hold on
% plot(ta(1:end),ka(1:end),'b-','LineWidth',1.6);
% plot(ta(1:end),sa(1:end),'r--','LineWidth',2);
% plot(ta(1:end),ka(1:end)+sa(1:end),'g-','LineWidth',2.1);
% xlabel('Time')
% ylabel('Energy')
% legend('kinetic','strain','total')
% %set(gca,'XTick',[0 0.5 1.0 1.5 2.0])
% axis([0 3.5 0 3])
%}
disp([num2str(toc),'   DONE '])
