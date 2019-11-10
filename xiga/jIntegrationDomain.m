function [Jdomain,qnode,radius] = jIntegrationDomain(tip_elem,xTip,node,elementV)
% determine elements for J integral computation
% (node,element): node and connectivity of visualization Q4 mesh
%

global controlPts element jDomainFac

numnode = size(node,1);
% -------------------------------------
% calculation of the area of the tip element
x = node(elementV(tip_elem(1),:),:);
% Area = sum of areas of each sub-triangle
x0 = x(1,1);
y0 = x(1,2);

x1 = x(2,1);
y1 = x(2,2);

x2 = x(3,1);
y2 = x(3,2);

x3 = x(4,1);
y3 = x(4,2);

A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;
A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
area = A1 + A2; 

% J radius = fac * sqrt(area);

radius = jDomainFac * sqrt(area)
center = xTip;

r=[];
% Distance from the center of tip element
for i = 1 : numnode
    sctr = node(i,:);
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
    r    = [r,rho];
end
test = r-radius;
test = test(elementV)'; % put nodal test into element test
                       % test(4,numelem) for Q4 elements
test = max(test).*min(test); % test(1,numelem): 
Jdomain = find(test<=0);

% determine nodal weight values
% use control points not vertices of Q4 mesh

% r=[];
% % Distance from the center of tip element
% for i = 1 : size(controlPts,1)
%     sctr = controlPts(i,:);
%     rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
%     r    = [r,rho];
% end

test1 = r-radius;
%test1 = test1(element(Jdomain,:))';
test1 = test1(elementV(Jdomain,:))'; % Q4 mesh!!!
test1 = (test1<=0); % nodes inside circle (<0) have weights=1
                    % nodes outside have weights=0
qnode = test1'; % qnode(no,4): no=no of elements in J domain

