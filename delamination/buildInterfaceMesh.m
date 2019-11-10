function iMesh = buildInterfaceMesh(mesh,no,t)
%
% Build mesh of interface elements.
% mesh: mesh of a trivariate NURBS solid
% no:   number of plies
%
% VP Nguyen, Cardiff University, 2013

noElemsU = mesh.noElemsU;
noElemsV = mesh.noElemsV;
p        = mesh.p;
q        = mesh.q;
uKnot    = mesh.uKnot;
vKnot    = mesh.vKnot;
noPtsX   = mesh.noPtsX;
noPtsY   = mesh.noPtsY;
controlPts=mesh.controlPts;


iElements   = zeros(noElemsU*noElemsV,2*(p+1)*(q+1));
iElementS   = buildIGA2DMesh (uKnot,vKnot,noPtsX,noPtsY,p,q);

e = 1;

% loop over plies
for ip=1:no-1
    y0 = t/no*ip;
    delaminationNodes  =  find(abs(controlPts(:,3) - y0 ) <1e-10);
    mm                 = 0.5*length(delaminationNodes);
    lowerNodes         = delaminationNodes(1:mm);
    upperNodes         = delaminationNodes(mm+1:end);
    
    for i=1:noElemsU*noElemsV
        sctr = iElementS.globElems(i,:);
        iElements(e,1:(p+1)*(q+1))     = lowerNodes(sctr);
        iElements(e,(p+1)*(q+1)+1:end) = upperNodes(sctr);
        e = e + 1;
    end
end

iMesh.locElems  = iElements;
iMesh.globElems = iElements;


function mesh = buildIGA2DMesh (uKnot,vKnot,noPtsX,noPtsY,p,q)

uniqueUKnots   = unique(uKnot);
uniqueVKnots   = unique(vKnot);

noElemsU       = length(uniqueUKnots)-1; % # of elements xi dir.
noElemsV       = length(uniqueVKnots)-1; % # of elements eta dir.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% chan = 
%        1 2 3 4
%        5 6 7 8
%        9 10 11 12
% for a 4x3 control points


chan           = zeros(noPtsY,noPtsX);

count = 1;

for i=1:noPtsY
    for j=1:noPtsX
        chan(i,j) = count;
        count = count + 1;
    end
end

% determine our element ranges and the corresponding 
% knot indices along each direction

[elRangeU,elConnU] = buildConnectivity(p,uKnot,noElemsU);
[elRangeV,elConnV] = buildConnectivity(q,vKnot,noElemsV);

% combine info from two directions to build the elements
% element is numbered as follows
%  5 | 6 | 7 | 8
% ---------------
%  1 | 2 | 3 | 4 
% for a 4x2 mesh

noElems = noElemsU * noElemsV;

element = zeros(noElems,(p+1)*(q+1));

e = 1;
for v=1:noElemsV
    vConn = elConnV(v,:);
    for u=1:noElemsU
        c = 1;
        uConn = elConnU(u,:);
        for i=1:length(vConn)
            for j=1:length(uConn)
              element(e,c) = chan(vConn(i),uConn(j));
              c = c + 1;
            end
        end
        e = e + 1;
    end      
end

index = zeros(noElems,2);
count = 1;

for j=1:size(elRangeV,1)
    for i=1:size(elRangeU,1)
        index(count,1) = i;
        index(count,2) = j;
        
        count = count + 1;
    end
end


[C,Cxi,Cet] = bezierExtraction2D(uKnot,vKnot,p,q);

mesh.p          = p;
mesh.q          = q;
mesh.uKnot      = uKnot;
mesh.vKnot      = vKnot;
mesh.noPtsX     = noPtsX;
mesh.noPtsY     = noPtsY;
mesh.C          = C;
mesh.locElems   = element;
mesh.globElems  = element;
mesh.rangeU     = elRangeU;
mesh.rangeV     = elRangeV;
mesh.index      = index;
mesh.noElemsU   = noElemsU;
mesh.noElemsV   = noElemsV;
mesh.elConnU    = elConnU;
mesh.elConnV    = elConnV;
mesh.elemCount  = noElems;



