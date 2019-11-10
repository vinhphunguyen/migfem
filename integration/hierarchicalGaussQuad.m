function [aa] = hierarchicalGaussQuad(order,phi,pts,gCoord,iVoid,level)

global W Q levelMax VOID

level = level + 1;

% first decompose the parent element [-1,1]x[-1,1]
% into 4 sub-cells

pt1 = pts(1,:);
pt2 = pts(2,:);
pt3 = pts(3,:);
pt4 = pts(4,:);
node  = square_node_array(pt1,pt2,pt3,pt4,3,3);
inc_u = 1;
inc_v = 3;
node_pattern = [ 1 2 3+2 3+1 ];
element      = make_elem(node_pattern,3-1,3-1,inc_u,inc_v);


% hold on
% plot_mesh(node,element,'Q4','g.-');
% %plot(Q(:,1),Q(:,2),'*');
% axis equal

levelset    = zeros(1,9);

for i=1:9
    pt  = node(i,:);
    [N] = lagrange_basis('Q4', pt);
    gpt = N' * gCoord;
    xc  = VOID(iVoid,1);
    yc  = VOID(iVoid,2);
    rc  = VOID(iVoid,3);
    levelset(i) = sqrt((gpt(1)-xc)^2+(gpt(2)-yc)^2)-rc;
end

% loop over sub-cells

for e = 1:4
    sctr  = element(e,:);
    ls    = levelset(sctr);
    coord = node(sctr,:);
    if ( max(ls)*min(ls) < 0 )   
        if level < levelMax
            [aa]  = hierarchicalGaussQuad(order,ls,coord,gCoord,iVoid,level);
        else
            [w,q] = quadrature(order,'GAUSS',2);
            
            % transform quadrature points into the parent element
            
            for n=1:length(w)
                [N,dNdxi] = lagrange_basis('Q4',q(n,:));
                J0=coord'*dNdxi;
                Q = [Q;N'*coord];
                W = [W;w(n)*det(J0)];
            end
        end
    else
        [w,q] = quadrature(order,'GAUSS',2);
        for n=1:length(w)
            [N,dNdxi] = lagrange_basis('Q4',q(n,:));
            J0=coord'*dNdxi;
            Q = [Q;N'*coord];
            W = [W;w(n)*det(J0)];
        end
    end
    aa=1;
end

% debug only

%plot(Q(:,1),Q(:,2),'*')




