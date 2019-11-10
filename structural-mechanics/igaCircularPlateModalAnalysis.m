%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for Kirchoff plate problems.
%
% Rotation-free thin plates. Fully clamped or simply supported
% circular plates.
%
% Modal analysis. Example given in paper of Cottrel, Reali, Hughes
% on IGA and structural vibration.
%
% Vinh Phu Nguyen,
% Cardiff University
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/

clc
clear all

global p q

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%  CAD input: control points and knot vectors
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

circularPlateData

% constitutive matrix

D  = E*t^3/(12*(1-nu^2));
C  = D*[1 nu 0;nu 1 0; 0 0 0.5*(1-nu)];

rho = 2.320;

I0  = rho*t;
I2  = rho*t^3/12;

numberOfModes     = 12;
numberOfModesPlot = 12;

% build connectivity ...

generateIGA2DMesh

noCtrPts       = noPtsX   * noPtsY;
noDofs         = noCtrPts * 1;

% find boundary nodes for boundary conditions

EPS = 1e-8;
bottomNodes  =  1:noPtsX;
topNodes     =  noPtsX*(noPtsY-1)+1:noCtrPts;
leftNodes    =  1:noPtsX:noPtsX*(noPtsY-1)+1;
rightNodes   =  noPtsX:noPtsX:noCtrPts;

fixedNodes  =  unique([bottomNodes;rightNodes;topNodes;leftNodes]);

if clamped
    nextToBotNodes = noPtsX+2:2*noPtsX-1;
    nextToRgtNodes = 2*noPtsX-1:noPtsX:noPtsX*(noPtsY-1)-1;
    nextToTopNodes = noPtsX*(noPtsY-2)+2:noPtsX*(noPtsY-1)-1;
    nextToLefNodes = noPtsX+2:noPtsX:noPtsX*(noPtsY-2)+2;
    
    nextNodes      = unique([nextToBotNodes';nextToRgtNodes';...
        nextToTopNodes';nextToLefNodes']);
    
    fixedNodes     = [fixedNodes; nextNodes(:)];
end

plot(controlPts(fixedNodes,1),controlPts(fixedNodes,2),...
    'bs','MarkerEdgeColor','r','MarkerSize',14);


% initialization

K = sparse(noDofs,noDofs);  % global stiffness matrix
M = sparse(noDofs,noDofs);  % global mass matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gauss quadrature rule
noGPs = p+1;
noGpEle = noGPs^2;
[W,Q]=quadrature(  noGPs, 'GAUSS', 2 ); % noGPs x noGPs point quadrature

% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM']);

% Loop over elements (knot spans)

for e=1:noElems
    idu    = index(e,1);
    idv    = index(e,2);
    xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
    etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = element(e,:);          %  element scatter vector
    nn     = length(sctr);
    pts    = controlPts(sctr,:);
    
    % loop over Gauss points
    
    for gp=1:size(W,1)
        pt      = Q(gp,:);
        wt      = W(gp);
        
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        
        % shape functions, first and second derivatives w.r.t natural coords
        
        [R dRdxi dRdeta dR2dxi dR2det dR2dxe] = ...
            NURBS2DBasis2ndDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        
        % compute the jacobian of physical and parameter domain mapping
        % then the derivative w.r.t spatial physical coordinates
        
        jacob  = [dRdxi; dRdeta]          * pts; % 2x2 matrix
        jacob2 = [dR2dxi; dR2det; dR2dxe] * pts; % 3x2 matrix
        
        J1    = det(jacob);
        
        dxdxi = jacob(1,1); dydxi = jacob(1,2);
        dxdet = jacob(2,1); dydet = jacob(2,2);
        
        j33   = [dxdxi^2     dydxi^2     2*dxdxi*dydxi;
            dxdet^2     dydet^2     2*dxdet*dydet;
            dxdxi*dxdet dydxi*dydet dxdxi*dydet+dxdet*dydxi];
        
        % Jacobian inverse and spatial 1st and 2nd derivatives
        
        invJacob   = inv(jacob);
        dRdx       = invJacob*[dRdxi;dRdeta];
        dR2dx      = inv(j33)*([dR2dxi; dR2det; dR2dxe]-jacob2*dRdx);
        
        % B matrix
        
        B          = dR2dx;
        B(3,:)     = B(3,:)*2;
        
        % compute elementary stiffness matrix and
        % assemble it to the global matrix
        wip = J1 * J2 * wt;
        K(sctr,sctr) = K(sctr,sctr) + B' * C * B * wip;
        M(sctr,sctr) = M(sctr,sctr) + (I0* R' * R + I2 * dRdx' * dRdx )* wip;
    end
end

% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM']);

activeDof=setdiff([1:noCtrPts]',[fixedNodes]);

[modeShape,freq]=eigs(K(activeDof,activeDof),M(activeDof,activeDof),...
    numberOfModes,0);

freq=diag(freq);%/(2*pi);  % frequency in kHz
freq=sqrt(freq);
[freq,ind]=sort(freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%% Visualization using a Q4 visualization mesh

buildVisualizationMesh;

disp([num2str(toc),'  PLOT MODE SHAPE']);

for mm=1:numberOfModesPlot
    m=ind(mm);
    disp(['   MODE: ',num2str(m),' ',num2str(freq(m))])
    % PLOT MODE SHAPE
    figure(m); 
    U            = zeros(noCtrPts,1);
    U(activeDof) = modeShape(:,m);
    scaleFactor  = 20/max(abs(U));
    
    displ   = zeros(noElems,size(elementV,2));
    
    for e=1:noElems
        idu    = index(e,1);
        idv    = index(e,2);
        xiE    = elRangeU(idu,:); % [xi_i,xi_i+1]
        etaE   = elRangeV(idv,:); % [eta_j,eta_j+1]
        
        sctr   = element(e,:);         %  element scatter vector
        pts    = controlPts(sctr,:);
        
        uspan = FindSpan(noPtsX-1,p,xiE(1), uKnot);
        vspan = FindSpan(noPtsY-1,q,etaE(1),vKnot);
        
        % loop over Gauss points
        
        gp = 1;
        for iv=1:2
            if (iv==2)
                xiE = sort(xiE,'descend');
            end
            for iu=1:2
                Xi  = xiE(iu);
                Eta = etaE(iv);
                [N dRdxi dRdeta] = NURBS2DBasisDersSpecial([Xi; Eta],...
                    p,q,uKnot,vKnot,weights',[uspan;vspan]);
                
                displ(e,gp)    = N * U(sctr);
                gp = gp +1;
            end
        end
    end
    
    X = zeros(4,noElemsV);
    Y = zeros(4,noElemsV);
    Z = displ';
    
    for i = 1:size(elementV,1)
        sctr   = elementV(i,:);
        X(:,i) = node(sctr,1);
        Y(:,i) = node(sctr,2);
    end
    
    fill3(X,Y,Z,Z);    
    title(['MODE ',num2str(m),', FREQUENCY = ',num2str(freq(m))])
    view(37,36)
    axis off
    %print(m, '-djpeg90', ['afm_mode_',num2str(m),'.jpg']);
end

% run close all to close all figures
% close all;

%% ANIMATE MODE
% Inspired by the FEM code of Jack Chessa.

nCycles = 5;     % number of cycles to animate
fpc     = 10;    % frames per cycle
fact    = sin(linspace(0,2*pi,fpc));
m       = input('What mode would you like to animate (type 0 to exit) ');

% while ( m~=0 )
%     U           = zeros(noCtrPts,1);
%     U(activeDof)= modeShape(:,m);
%     wt=20/max(abs(U));
%     for i=1:fpc
%         scaleFactor=fact(i)*wt;
%         figure(length(freq+1));
%         clf;
%         plot_field(node+[U(1:numnode) U(numnode+1:2*numnode)
%             U(2*numnode+1:3*numnode)]*scaleFactor,topSurface,elemType{topID},...
%             ones(3*numnode,1));
%         hold on
%         plot_mesh(node+[U(1:numnode) U(numnode+1:2*numnode)
%             U(2*numnode+1:3*numnode)]*scaleFactor,topSurface,elemType{topID},?k-?);
%         plot_mesh(node,topSurface,elemType{topID},?w-?);
%         hold on
%         view(37,36)
%         axis([70 240 30 160 -10 10])
%         title([?MODE ?,num2str(m),?, FREQUENCY = ?,num2str(freq(m)),? [kHz]?])
%         axis off
%         film(i)=getframe;
%     end
%     movie(film,nCycles);
%     m=input(?What mode would you like to animate (type 0 to exit) ?);
%     if ( m > length(freq) )
%         disp([?mode must be less than ?,num2str(length(freq))])
%     end
% end

opts = struct('Color','rgb','Bounds','tight','FontMode','fixed','FontSize',20);
%exportfig(gcf,fileName,opts)


% Exact solution check

omega1 = 1.015^2*pi^2/r^2*sqrt(D/rho/t)
omega2 = 1.468^2*pi^2/r^2*sqrt(D/rho/t)
omega3 = 2.007^2*pi^2/r^2*sqrt(D/rho/t)







