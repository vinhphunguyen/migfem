function vMesh=buildVMesh3DInterface(solid)
%
% build a H8 mesh (eight node brick elements) from image of knots.
% this mesh is used for visualization purpose.
% Standard FE visualization function like plot_field can then be reused.
% Target: NURBS surfaces with C^0 lines for delamination problems.
% Assumption: C^0 lines along zeta direction.
%
% Vinh Phu Nguyen
% Cardiff University

[p,q,r,uKnot,vKnot,wKnot,noPtsX,noPtsY,noPtsZ,weights,controlPts] = ...
    convert3DNURBS (solid);

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);
wKnotVec = unique(wKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);
noKnotsW = length(wKnotVec);

multiplicity  = histc(wKnot,wKnotVec); % find duplicated knots
c0LineCnt     = sum(multiplicity(2:end-1)==r+1); % number of discontinuity lines

%  using homogeneous coordinates

projcoord = [controlPts(:,1).* weights ...
    controlPts(:,2).* weights ...
    controlPts(:,3).* weights ...
    weights];

dim=size(projcoord,2);

node  = zeros(noKnotsU*noKnotsV*noKnotsW,3);
count = 1;

for wk=1:noKnotsW
    zeta = wKnotVec(wk);
    for vk=1:noKnotsV
        eta = vKnotVec(vk);
        for uk=1:noKnotsU
            xi = uKnotVec(uk);
            
            tem = SolidPoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                noPtsZ-1,r,wKnot,projcoord,dim,xi,eta,zeta);
            
            node(count,1) = tem(1)/tem(4);
            node(count,2) = tem(2)/tem(4);
            node(count,3) = tem(3)/tem(4);
            
            count         = count + 1;
        end
    end
    % if multiplicity = order+1 => C^0 line
    if (multiplicity(wk) == r+1) && (wk ~= 1) && (wk ~= noKnotsW)
        for vk=1:noKnotsV
            eta = vKnotVec(vk);
            for uk=1:noKnotsU
                xi = uKnotVec(uk);
                
                tem = SolidPoint(noPtsX-1,p,uKnot,noPtsY-1,q,vKnot,...
                    noPtsZ-1,r,wKnot,projcoord,dim,xi,eta,zeta);
                
                node(count,1) = tem(1)/tem(4);
                node(count,2) = tem(2)/tem(4);
                node(count,3) = tem(3)/tem(4);
                
                count         = count + 1;
            end
        end
    end
end

% build H8 elements

chan  = zeros(noKnotsU,noKnotsV,noKnotsW+c0LineCnt);

count = 1;
m     = 0;
for i=1:noKnotsW
    if (multiplicity(i) == r+1) && (i ~= 1) && (i ~= noKnotsW)
        for j=1:noKnotsV
            for k=1:noKnotsU
                chan(m+1,j,k) = count;
                count       = count + 1;
            end
        end
        for j=1:noKnotsV
            for k=1:noKnotsU
                chan(m+2,j,k) = count;
                count       = count + 1;
            end
        end
        m = m + 2;
    else
        m = m + 1;
        for j=1:noKnotsV
            for k=1:noKnotsU
                chan(m,j,k) = count;
                count       = count + 1;
            end
        end
    end
end

connecU = zeros(noKnotsU-1,2);
connecV = zeros(noKnotsV-1,2);
connecW = zeros(noKnotsW-1,2);

for i=1:size(connecU,1)
    connecU(i,:) = [i i+1];
end

for i=1:size(connecV,1)
    connecV(i,:) = [i i+1];
end

e=1;
for i=1:noKnotsW-1
    if (multiplicity(i) == r+1) && (i ~= 1)
        e = e + 1;
    end
    connecW(i,:) = e:e+1;
    e = e + 1;
end

noElems  = (noKnotsU-1) * (noKnotsV-1) * (noKnotsW-1);
elementV = zeros(noElems,8);

e = 1;
for w=1:noKnotsW-1
    wConn = connecW(w,:);
    for v=1:noKnotsV-1
        vConn = connecV(v,:);
        for u=1:noKnotsU-1
            c = 1;
            uConn = connecU(u,:);
            
            for i=1:length(wConn)
                for j=1:length(vConn)
                    for k=1:length(uConn)
                        elementV(e,c) = chan(wConn(i),vConn(j),uConn(k));
                        c = c + 1;
                    end
                end
            end
            e = e + 1;
        end
    end
end

% renumbering nodes according to Jack Chessa's code

% col3 = elementV(:,3);
% col4 = elementV(:,4);
% col7 = elementV(:,7);
% col8 = elementV(:,8);
% 
% elementV(:,3) = col4;
% elementV(:,4) = col3;
% elementV(:,7) = col8;
% elementV(:,8) = col7;

elementV(:,[3 4 7 8]) = elementV(:,[4 3 8 7]);


vMesh.node    =  node;
vMesh.element = elementV;



