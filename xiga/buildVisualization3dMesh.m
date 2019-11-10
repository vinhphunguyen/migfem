% function [node,elementV]=buildVisualization3dMesh(controlPts,weights,...
%                         uKnot,vKnot,wKnot,p,q,r)
%
% build a H8 mesh (eight node brick elements) from image of knots.
% this mesh is used for visualization purpose.
% Standard FE visualization function like plot_field
% can then be reused.
% Vinh Phu Nguyen
% Johns Hopkins University

% get rid of zero measure knot spans

uKnotVec = unique(uKnot);
vKnotVec = unique(vKnot);
wKnotVec = unique(wKnot);

% number of distinct knot values

noKnotsU = length(uKnotVec);
noKnotsV = length(vKnotVec);
noKnotsW = length(wKnotVec);

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
end

% build H8 elements

chan  = zeros(noKnotsU,noKnotsV,noKnotsW);

count = 1;

for i=1:noKnotsW
    for j=1:noKnotsV
        for k=1:noKnotsU
            chan(i,j,k) = count;
            count       = count + 1;
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

for i=1:size(connecW,1)
   connecW(i,:) = [i i+1];	
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

col3 = elementV(:,3);
col4 = elementV(:,4);
col7 = elementV(:,7);
col8 = elementV(:,8);

elementV(:,3) = col4;
elementV(:,4) = col3;
elementV(:,7) = col8;
elementV(:,8) = col7;





