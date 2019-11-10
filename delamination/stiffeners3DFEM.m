%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Isogeometric analysis for three dimensional elasticity problems.
% Using the so-called Bezier extraction operator.
%
% Cuvred composite panel with two stiffeners. 3D version.
%
% Vinh Phu Nguyen,
% Cardiff University, UK
% April 2013.
% nvinhphu@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../fem_util/
addpath ../C_files/
addpath ../data/
addpath ../meshing/
addpath ../integration/
addpath ../post-processing/
addpath ../fem-functions/
addpath ../meshing/
addpath ../nurbs-util/
addpath ../nurbs-geopdes/inst/
addpath ../iga/

clc

% input file
stiffener2;

% Dirichlet nodes

fixedXNodes  =  data.xnodes'; % transpose to make it a row vector
fixedYNodes  =  data.ynodes';
fixedZNodes  =  data.znodes';

fixedXNodes = [fixedXNodes fnodes'];

udofs = 3*fixedXNodes-2;    % global indecies  of the fixed x disps
vdofs = 3*fixedYNodes-1;    % global indecies  of the fixed y disps
wdofs = 3*fixedZNodes;      % global indecies  of the fixed z disps

uFixed = [zeros(size(data.xnodes')) ones(size(fnodes'))];
vFixed = zeros(size(fixedYNodes));
wFixed = zeros(size(fixedZNodes));

%% external force

fbar= 1;
b0  = -10;

%% initialization
nDim       = 3;
patchCount = length(data.mesh);
noDofs     = data.pntCount * nDim;

K = sparse(noDofs,noDofs);  % global stiffness matrix
f = zeros(noDofs,1);        % external force vector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assembling system of equation
% Stiffness matrix and external force vector

disp([num2str(toc),'  ASSEMBLING THE SYSTEM'])

% Loop over patches
for ip=1:patchCount
    mesh      = data.mesh{ip};
    globElems = mesh.globElems;
    locElems  = mesh.locElems;
    C         = mesh.C;
    weights   = mesh.weights;
    p         = mesh.p;
    q         = mesh.q;
    r         = mesh.r;
    controlPts= mesh.controlPts;
%    elementSet=elemSet{ip};
    
    % Gauss quadrature rule
    [W,Q] = gaussianQuadNURBS(p+1,q+1,r+1);
    
    %% Pre-compute Bernstein basis and derivatives for ONE Bezier element
    
    noBasis = (p+1)*(q+1)*(r+1);
    noGpEle = noBasis;
    
    shapes  = zeros(noGpEle,noBasis);
    derivs  = zeros(noGpEle,noBasis,3);
    
    for gp=1:size(W,1)
        [shapes(gp,:) derivs(gp,:,:)] = ...
            getShapeGradBernstein3D(p,q,r,Q(gp,1),Q(gp,2),Q(gp,3));
    end
    
    % Loop over elements
    for e=1:size(globElems,1)
        sctrg   = globElems(e,:);         %  global element scatter vector
        sctrl   = locElems (e,:);         %  local element scatter vector
        nn      = length(sctrg);
        
        sctrB(1:3:3*nn)    = 3*sctrg-2;
        sctrB(2:3:3*nn)    = 3*sctrg-1;
        sctrB(3:3:3*nn)    = 3*sctrg-0;
        
        
        Ce     = C(:,:,e);              % element Bezier extraction operator
        we     = diag(weights(sctrl));  % element weights
        pts    = controlPts(sctrl,:);   % element nodes
        Wb     = Ce'*weights(sctrl);    % element Bezier weights
        
        mat = materials{1};      
%         if     ismember(e,elementSet{1})
%             mat = materials{1};
%         elseif ismember(e,elementSet{2})
%             mat = materials{2};
%         elseif ismember(e,elementSet{3})
%             mat = materials{3};
%         elseif ismember(e,elementSet{4})
%             mat = materials{4};
%         end
        
        % loop over Gauss points
        for gp=1:size(W,1)
            pt      = Q(gp,:);
            wt      = W(gp);
            
            %% Bernstein basis and derivatives at GP gp
            Be      = shapes(gp,:)';
            dBedxi  = reshape(derivs(gp,:,:),noBasis,3);
            
            %% Bezier weight functions (denomenator of NURBS)
            wb        = dot(Be,Wb);            % Be(I)*Wb(I)
            dwbdxi(1) = dot(dBedxi(:,1),Wb);   % Be(I)_{,xi} * Wb(I)
            dwbdxi(2) = dot(dBedxi(:,2),Wb);   % Be(I)_{,et} * Wb(I)
            dwbdxi(3) = dot(dBedxi(:,3),Wb);   % Be(I)_{,et} * Wb(I)
            %% Shape function and derivatives
            R          = we*Ce*Be/wb;
            dRdxi(:,1) = we*Ce*(dBedxi(:,1)/wb-dwbdxi(1)*Be/(wb*wb));
            dRdxi(:,2) = we*Ce*(dBedxi(:,2)/wb-dwbdxi(2)*Be/(wb*wb));
            dRdxi(:,3) = we*Ce*(dBedxi(:,3)/wb-dwbdxi(3)*Be/(wb*wb));
            
            %% Jacobian matrix
            dxdxi = pts'*dRdxi;
            
            dxidx = inv(dxdxi);
            dRdx  = dRdxi*dxidx;
            detJ  = det(dxdxi);
            
            % B matrix
            B = strainDispMatrix3d(nn,dRdx);
            
            % compute elementary stiffness matrix and
            % assemble it to the global matrix
            
            K(sctrB,sctrB) = K(sctrB,sctrB) + B' * mat.stiffMat * B * detJ * wt;
            
            %body force acting upon the skin solid
%                         if ip == 1
%                             sctry    = 3*sctrg-1;
%                             f(sctry) = f(sctry) + R * b0 * detJ * wt;
%                         end
        end
    end
end

%% contribution of interface elements
% rigid links
w     = 1e6;
penaltyStiffness = w*[1 -1;-1 1];

for j=1:size(aa,1)
    %for i=1:size(aa,1)
        sctr  = [aa(j) bb(j)];
        sctrx = 3*sctr-2;
        sctry = 3*sctr-1;
        sctrz = 3*sctr;
        K(sctrx,sctrx) = K(sctrx,sctrx) + penaltyStiffness;
        K(sctry,sctry) = K(sctry,sctry) + penaltyStiffness;
        K(sctrz,sctrz) = K(sctrz,sctrz) + penaltyStiffness;
    %end
end

%% Computing external force

p       = data.mesh{1}.r;
q       = data.mesh{1}.q;
uKnot   = data.mesh{1}.wKnot;
vKnot   = data.mesh{1}.vKnot;
weights = ones(data.mesh{1}.noPtsZ*data.mesh{1}.noPtsY);
controlPts = data.mesh{1}.controlPts;
[W,Q] = gaussianQuad2DNURBS(q+1,r+1);

for e=1:size(bndElems1,1)
    idu    = index1(e,1);
    idv    = index1(e,2);
    xiE    = data.mesh{1}.rangeW(idu,:); % [xi_i,xi_i+1]
    etaE   = data.mesh{1}.rangeV(idv,:); % [eta_j,eta_j+1]
    
    sctr   = bndElems1(e,:);          %  element scatter vector
    nn     = length(sctr);
    sctrFx = 3*sctr-2;                % force in x-direction
    sctrFy = 3*sctr-1;                % force in y-direction
    sctrFz = 3*sctr;                  % force in z-direction
    pts    = controlPts(sctr,:);            
    for gp=1:size(W,1)% loop over Gauss points
        pt      = Q(gp,:);
        wt      = W(gp);
        % compute coords in parameter space
        Xi      = parent2ParametricSpace(xiE,pt(1));
        Eta     = parent2ParametricSpace(etaE,pt(2));
        J2      = jacobianPaPaMapping(xiE,etaE);
        % compute derivatives of shape functions
        [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
        jacob      = [dRdxi; dRdeta] * pts;
        a1    = jacob(1,:);
        a2    = jacob(2,:);
        a3    = cross(a1,a2); 
        norma = norm(a3);     a3=a3/norma ;  
        J1    = norma;        
        J          = J1 * J2;
        %f(sctrFx)  = f(sctrFx) + R' * fbar * a3(1) * J * wt;
        %f(sctrFy)  = f(sctrFy) + R' * fbar * a3(2) * J * wt;
        %f(sctrFz)  = f(sctrFz) + R' * fbar * a3(3) * J * wt;
    end
end

% for e=1:size(bndElems2,1)
%     idu    = index2(e,1);
%     idv    = index2(e,2);
%     xiE    = data.mesh{1}.rangeW(idu,:); % [xi_i,xi_i+1]
%     etaE   = data.mesh{1}.rangeV(idv,:); % [eta_j,eta_j+1]
%     
%     sctr   = bndElems2(e,:);          %  element scatter vector
%     nn     = length(sctr);
%     sctrF  = 3*sctr-2;                % force in x-direction
%     pts    = controlPts(sctr,:);
%     
%     % loop over Gauss points
%     
%     for gp=1:size(W,1)
%         pt      = Q(gp,:);
%         wt      = W(gp);
%         % compute coords in parameter space
%         Xi      = parent2ParametricSpace(xiE,pt(1));
%         Eta     = parent2ParametricSpace(etaE,pt(2));
%         J2      = jacobianPaPaMapping(xiE,etaE);
%         % compute derivatives of shape functions
%         [R dRdxi dRdeta] = NURBS2DBasisDers([Xi; Eta],p,q,uKnot,vKnot,weights');
%         jacob      = [dRdxi; dRdeta] * pts;
%         a1    = jacob(1,:);
%         a2    = jacob(2,:);
%         a3    = cross(a1,a2);                 
%         J1    = norma;
%         J        = J1 * J2;
%         f(sctrF) = f(sctrF) - R' * fbar * J * wt;
%     end
% end

%% Apply BCs
disp([num2str(toc),'  APPLYING BOUNDARY  CONDITIONS'])

bcwt=mean(diag(K)); % a measure of the average  size of an element in K
% used to keep the  conditioning of the K matrix

f=f-K(:,udofs)*uFixed';  % modify the  force vector
f=f-K(:,vdofs)*vFixed';
f=f-K(:,wdofs)*wFixed';

f(udofs) = bcwt*uFixed;
f(vdofs) = bcwt*vFixed;
f(wdofs) = bcwt*wFixed;

K(udofs,:)=0;  % zero out the rows and  columns of the K matrix
K(vdofs,:)=0;
K(wdofs,:)=0;
K(:,udofs)=0;
K(:,vdofs)=0;
K(:,wdofs)=0;

K(udofs,udofs)=bcwt*speye(length(udofs));  % put ones*bcwt on the diagonal
K(vdofs,vdofs)=bcwt*speye(length(vdofs));
K(wdofs,wdofs)=bcwt*speye(length(wdofs));

%% SOLVE SYSTEM
disp([num2str(toc),'  SOLVING THE SYSTEM'])
U=K\f;

umatrix=[U(1:3:end) U(2:3:end) U(3:3:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% POST-PROCESSING
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([num2str(toc),'  POST-PROCESSING'])

vtuFile0 = '../results/two-stiffeners';

%elemSet=[];
% Loop over patches
for ip=1:patchCount
    vtuFile = strcat(vtuFile0,num2str(ip));
    figure; hold on;
    ok      = plotStress3DForPatch(data,ip,vtuFile,U,elemSet,materials);
end


pvdFile = fopen(strcat('../results/',vtuFile0,'.pvd'), 'wt');

fprintf(pvdFile,'<VTKFile byte_order="LittleEndian" type="Collection" version="0.1">\n');
fprintf(pvdFile,'<Collection>\n');

for i = 1:patchCount
    vtuFile = sprintf('%s%d%s',vtuFile0,i,'.vtu');
    fprintf(pvdFile,'<DataSet file=''%s'' groups='''' part=''%d'' timestep=''1''/>\n',vtuFile,i);
end

fprintf(pvdFile,'</Collection>\n');
fprintf(pvdFile,'</VTKFile>\n');

fclose(pvdFile);







