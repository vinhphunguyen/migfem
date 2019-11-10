function [nodes] = assembly_GMPM(e,nnx,nny)

nodes = [];
elements = [];
global element elementmatrix elementpos
        
        pos = elementpos(e,:);
        pos(1) = pos(1)-1;
        pos(2) = pos(2)-1;
        for i=1:3
            for j=1:3
                ii = (i-1)+pos(1);
                jj = (j-1)+pos(2);

                if ii>0 && ii<=(nny-1)
                   if jj>0 && jj<=(nnx-1)
                       elements = [elements , elementmatrix(ii,jj) ];
                       nodes = [nodes , element( elementmatrix(ii,jj),:) ];
                       %{
                       for q=1:4
                           element(i,q)
                           for k=1:size(nodes,2)
                              if nodes(k)== element(i,q)
                              else
                                 nodes(k) = [nodes , element(i,q) ];
                              end
                           end
                       end
                       %}
                   end
                end
            end
        end
        nodes = unique(nodes);
        sctrB = zeros(1,2*size(nodes,2));
        for i=1:size(nodes,2)
                sctrB(2*i-1) = 2*nodes(i)-1 ;
                sctrB(2*i)   = 2*nodes(i)   ;
        end
end

    