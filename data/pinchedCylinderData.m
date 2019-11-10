% Data for the pinched cylinder problem
% Quadratic NURBS for 3 directions
% Vinh Phu Nguyen
% Delft University of Technology/Johns Hopkins University

tanPi8 = tan(pi/8);
R      = 300;
L      = 300;
t      = 3;
Ri     = R -t;
Rii    = 0.5*(R+Ri);

controlPts=[R 0 0;R R*tanPi8 0; R*tanPi8 R 0; 0 R 0;
    Rii 0 0;Rii Rii*tanPi8 0; Rii*tanPi8 Rii 0; 0 Rii 0;
    Ri 0 0;Ri Ri*tanPi8 0; Ri*tanPi8 Ri 0; 0 Ri 0;
    ...
    R 0 0.5*L;R R*tanPi8 0.5*L; R*tanPi8 R 0.5*L; 0 R 0.5*L;
    Rii 0 0.5*L;Rii Rii*tanPi8 0.5*L; Rii*tanPi8 Rii 0.5*L; 0 Rii 0.5*L;
    Ri 0 0.5*L;Ri Ri*tanPi8 0.5*L; Ri*tanPi8 Ri 0.5*L; 0 Ri 0.5*L;
    ...
    R 0 L;R R*tanPi8 L; R*tanPi8 R L; 0 R L;
    Rii L L;Rii Rii*tanPi8 L; Rii*tanPi8 Rii L; 0 Rii L;
    Ri 0 L;Ri Ri*tanPi8 L; Ri*tanPi8 Ri L; 0 Ri L];

p = 2;
q = 2;
r = 2;

% only one quadratic element along the thickness (q)
% uKnot will be refined using knot insertion

uKnot = [0 0 0 0.5 1 1 1];
vKnot = [0 0 0 1 1 1];
wKnot = [0 0 0 1 1 1];

noPtsX = length(uKnot)-p-1;
noPtsY = length(vKnot)-q-1;
noPtsZ = length(wKnot)-r-1;

weights = ones(1,noPtsX*noPtsY*noPtsZ)';

fac = 0.5*(1+1/sqrt(2));

weights([2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35]) = fac;

% refinement along p direction

refineCount = 5;

for i=1:refineCount
    oldControlPts = controlPts;
    
    uKnotVectorU = unique(uKnot);
    uKnotVectorV = unique(vKnot);
    uKnotVectorW = unique(wKnot);
    
    % new knots along three directions
    
    newKnotsX = uKnotVectorU(1:end-1) + 0.5*diff(uKnotVectorU);
    newKnotsY = uKnotVectorV(1:end-1) + 0.5*diff(uKnotVectorV);
    newKnotsZ = uKnotVectorW(1:end-1) + 0.5*diff(uKnotVectorW);
    
    %% h-refinement (NURBS) in x-direction
    dim = size(controlPts,2);
    
    nonewkX      = size(newKnotsX,2);
    newprojcoord = zeros(noPtsX*noPtsY+nonewkX*noPtsY,dim+1);
    
    rstart = 1;
    wstart = 1;
    
    for j=1:noPtsY
        rstop = rstart + noPtsX-1;
        wstop = wstart + noPtsX-1 + nonewkX;
        
        locCP        = controlPts(rstart:rstop,:);
        locweights   = weights   (rstart:rstop);
        locprojcoord = nurb2proj(noPtsX, locCP, locweights);
        
        % refinement of x
        [tempknotVectorX,tempControlPts] = ...
            RefineKnotVectCurve(noPtsX-1,p,uKnot,locprojcoord,newKnotsX,nonewkX-1);
        
        newprojcoord(wstart:wstop,:)=tempControlPts;
        wstart = wstop+1;
        rstart = rstop+1;
    end
    
    uKnot                 = tempknotVectorX;
    [controlPts, weights] = proj2nurbs(newprojcoord);
    noPtsX                = noPtsX+nonewkX;
end

% refinement along r direction

wKnot = uKnot;

noPtsZ = length(wKnot)-r-1;

for i=1:noPtsZ-1
    a = i*L/(noPtsZ-1);
    controlPts=[controlPts;
        [controlPts(1:noPtsX,1:2)          a*ones(noPtsX,1)];
        [controlPts(noPtsX+1:2*noPtsX,1:2) a*ones(noPtsX,1)];
        [controlPts(2*noPtsX+1:3*noPtsX,1:2) a*ones(noPtsX,1)];
        ];
    
    weights=[weights;
        weights(1:noPtsX);weights(1:noPtsX);weights(1:noPtsX);        
        ];
end
%%%

% buildVisualization3dMesh
% 
% figure
% hold on
% plot_mesh(node,elementV,'B8','b.-');
% view(3)
% plot3(controlPts(:,1),controlPts(:,2),controlPts(:,3),'r*');

% numNode = size(node,1);
% 
% % normal stresses
% sigmaXX = zeros(numNode,1);
% sigmaYY = zeros(numNode,1);
% sigmaZZ = zeros(numNode,1);
% 
% % shear stresses
% sigmaXY = zeros(numNode,1);
% sigmaYZ = zeros(numNode,1);
% sigmaZX = zeros(numNode,1);
% 
% % displacements
% dispX   = zeros(numNode,1);
% dispY   = zeros(numNode,1);
% dispZ   = zeros(numNode,1);
% 
% 
% VTKPostProcess3d(node,elementV,'B8',...
%     [sigmaXX sigmaYY sigmaZZ sigmaXY sigmaYZ sigmaZX],[dispX dispY dispZ]);
