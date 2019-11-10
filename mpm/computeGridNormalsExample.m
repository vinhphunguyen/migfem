% This file illustrates how to compute the normals at grid node
% in contact MPM simulations.
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

tic;

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

%% generate material points

noParticleX = 2;
noParticleY = 2;

[W,Q]=quadrature(  2, 'GAUSS', 2 ); % 2x2 Gaussian quadrature


% body1

center = [0.5 0.5];
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


Mp  = mass;                           % mass
Vp  = volume;                         % volume
Fp  = repmat([1 0 0 1],pCount,1);     % gradient deformation
s   = zeros(pCount,3);                % stress
eps = zeros(pCount,3);                % strain
vp  = [velo1];                        % velocity
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

%% cell density computation

cellDensity = zeros(elemCount,1);

for ie=1:elemCount
    neighbors = getNeighbors(ie, numx2, numy2);
    center    = 1/4*sum( node(element(ie,:),:) );
    for in=1:length(neighbors)
        elemId = neighbors(in);
        mpts  = mpoints{elemId};
        for p=1:length(mpts)
            pid  = mpts(p);
            x    = xp(pid,:);
            [phi,dphi]=getQuadraticBspline2D(x-center,deltax,deltay);
            cellDensity(ie) = cellDensity(ie) + phi*Mp(pid);
        end        
    end
    cellDensity(ie) = cellDensity(ie)/(deltax*deltay);
end

%% only compute normals for boundary elements
%% use level sets to define them

r  = 0.2; % radius
xc = 0.5; % x coord of center
yc = 0.5; % y coord of center
ls = zeros(length(node),1);

for i = 1 : length(node)
    x      = node(i,1);
    y      = node(i,2);
    d      = sqrt((x-xc)^2+(y-yc)^2);
    ls(i)  = d - r; % level set
end


bndElems = [];

% loop over elements

for iel = 1 : length(element)
    sctr    = element(iel,:);
    phi     = ls(sctr);
    mpts  = mpoints{iel};
    if    ( max(phi)*min(phi) < 0 ) && ~isempty(mpts)
        bndElems = [bndElems;iel];
    end
end

bndNodes = unique(element(bndElems,:));

%% Now, compute normals at nodes of boundary elements

normals = zeros(length(node),2);

for iel = 1 : length(bndElems)
    eId    = bndElems(iel);
    sctr   = element(eId,:);
    dens   = cellDensity(eId);
    center = 1/4*sum( node(element(eId,:),:) );
    for in=1:4
        nId = sctr(in);
        levelS = ls(nId);
        %if (levelS >= 0 )
        xI  = node(nId,:);
        [N,dNdx]=getMPM2D(center-xI,deltax,deltay);
        normals(nId,:) = normals(nId,:) + dNdx*dens;
        %end
    end
end

%% plot mesh, particles

hold on
plot_mesh(node,element,'Q4','k-',1.);
plot_mesh(node,element(bndElems,:),'Q4','cy-',2.1);
plot(xp(:,1),xp(:,2),'k.','markersize',10);
for i=1:length(bndNodes)
    nid = bndNodes(i);
    xI  = node(nid,:);
    nI  = normals(nid,:);
    nI = nI / norm(nI);
    le = 0.2;
    plot([xI(1) xI(1)+le*nI(1)],[xI(2) xI(2)+le*nI(2)],'r-','LineWidth',2);
end
axis off

%%


disp([num2str(toc),'   DONE '])
