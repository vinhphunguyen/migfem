function [aa] = hierarchicalGaussQuad(order,phi,pts,level)

global W Q levelMax

level = level + 1;

% first decompose the parent element [-1,1]x[-1,1]
% into 4 sub-cells

pt1 = pts(1,:);
pt2 = pts(2,:);
pt3 = pts(3,:);
pt4 = pts(4,:);
node    = square_node_array(pt1,pt2,pt3,pt4,3,3);
numnode = size(node,1);
inc_u = 1;
inc_v = 3;
node_pattern = [ 1 2 3+2 3+1 ];
element      = make_elem(node_pattern,3-1,3-1,inc_u,inc_v);


% hold on
% plot_mesh(node,element,'Q4','g.-');
% %plot(Q(:,1),Q(:,2),'*');
% axis equal

levelset    = zeros(1,9);
levelset(1) = phi(1);
levelset(3) = phi(2);
levelset(7) = phi(4);
levelset(9) = phi(3);

N = lagrange_basis('Q4',0.5*[pt1(1)+pt2(1) pt1(2)+pt2(2)]);
levelset(2) = N' * phi';

N = lagrange_basis('Q4',0.5*[pt1(1)+pt4(1) pt1(2)+pt4(2)]);
levelset(4) = N' * phi';

N = lagrange_basis('Q4',0.25*[pt1(1)+pt2(1)+pt3(1)+pt4(1)...
    pt1(2)+pt2(2)+pt3(2)+pt4(2)]);
levelset(5) = N' * phi';

N = lagrange_basis('Q4',0.5*[pt2(1)+pt3(1) pt2(2)+pt3(2)]);
levelset(6) = N' * phi';

N = lagrange_basis('Q4',0.5*[pt3(1)+pt4(1) pt3(2)+pt4(2)]);
levelset(8) = N' * phi';

% loop over sub-cells

for e = 1:4
    sctr  = element(e,:);
    ls    = levelset(sctr);
    coord = node(sctr,:);
    area  = (coord(3,1) - coord(1,1))*(coord(3,2) - coord(1,2));
    
    if ( max(ls)*min(ls) < 0 )
        if level < levelMax
            [aa] = hierarchicalGaussQuad(order,ls,coord,level);
        else
            [w,q] = quadrature(order+2,'GAUSS',2);
            
            % transform quadrature points into the parent element
            
            for n=1:length(w)
                N = lagrange_basis('Q4',q(n,:));
                Q = [Q;N'*coord];
                W = [W;w(n)*area];
            end
        end
    else
        [w,q] = quadrature(order+2,'GAUSS',2);
        for n=1:length(w)
            N = lagrange_basis('Q4',q(n,:));
            Q = [Q;N'*coord];
            W = [W;w(n)*area];
        end
    end
    aa=1;
end

% debug only

%plot(Q(:,1),Q(:,2),'*')




