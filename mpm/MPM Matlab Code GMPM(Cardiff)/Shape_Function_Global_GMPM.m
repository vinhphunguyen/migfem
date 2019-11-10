function [dNdX,N] = Shape_Function_Global_GMPM(nodes,xp,dim,dx,dy,shapefun)

global node

switch shapefun
    
    case 'GMPM'

    lpx = dim(1)/2;
    lpy = dim(2)/2;
    num_nodes = size(nodes,2);
    N = zeros(num_nodes,1);
    dNdX = zeros(num_nodes,2);

    for i=1: num_nodes
        distx = xp(1) - node(nodes(i),1);
        disty = xp(2) - node(nodes(i),2);

              if ( (-dx-lpx)<distx ) && ( distx<=(-dx+lpx) )

                   Nx = ((dx+lpx+distx)^2)/(4*dx*lpx);
                   dNx = (dx+lpx+distx)/(2*dx*lpx);

               elseif ( (-dx+lpx)<distx ) && ( distx<=(-lpx) )

                   Nx = 1+distx/dx;
                   dNx = 1/dx;

               elseif ( (-lpx)<distx ) && ( distx<=(lpx) )

                   Nx = 1-( (distx^2) + (lpx^2) )/(2*dx*lpx);
                   dNx = -distx/(dx*lpx);

               elseif ( (lpx)<distx ) && ( distx<=(dx-lpx) )

                   Nx = 1-distx/dx ; 
                   dNx = -1/dx;

               elseif ( (dx-lpx)<distx ) && ( distx<=(dx+lpx) )

                   Nx = ((dx+lpx-distx)^2)/(4*dx*lpx);
                   dNx = -(dx+lpx-distx)/(2*dx*lpx);
               else
                   Nx = 0;
                   dNx = 0;

              end

              if ( (-dy-lpy)<disty ) && ( disty<=(-dy+lpy) )

                   Ny = ((dy+lpy+disty)^2)/(4*dy*lpy);
                   dNy = (dy+lpy+disty)/(2*dy*lpy);

               elseif ( (-dy+lpy)<disty ) && ( disty<=(-lpy) )

                   Ny = 1+disty/dy;
                   dNy = 1/dy;

               elseif ( (-lpy)<disty ) && ( disty<=(lpy) )

                   Ny = 1-( (disty^2) + (lpy^2) )/(2*dy*lpy);
                   dNy = -disty/(dy*lpy);

               elseif ( (lpy)<disty ) && ( disty<=(dy-lpy) )

                   Ny = 1-disty/dy ; 
                   dNy = -1/dy;

               elseif ( (dy-lpy)<disty ) && ( disty<=(dy+lpy) )

                   Ny = ((dy+lpy-disty)^2)/(4*dy*lpy);
                   dNy = -(dy+lpy-disty)/(2*dy*lpy);
               else
                   Ny = 0;
                   dNy = 0;

              end

              Ni = Nx * Ny ;
              dNdxi = Ny * dNx ;
              dNdyi = Nx * dNy ;
              N(i,1) = Ni;
              dNdX(i,:) = [dNdxi , dNdyi];


    end

         B(1,1:2:(2*num_nodes)) = dNdX(:,1)';
         B(2,2:2:(2*num_nodes)) = dNdX(:,2)';
         B(3,1:2:(2*num_nodes)) = dNdX(:,2)';
         B(3,2:2:(2*num_nodes)) = dNdX(:,1)';
         Nfem(1,1:2:(2*num_nodes)) = N(:,1)';
         Nfem(2,2:2:(2*num_nodes)) = N(:,1)';
         
    case 'MPM'
     
          B = zeros(3,8);
          Nfem = zeros(2,8);
          node_elm = zeros(4,2);
          node_elm(:,1) = node(nodes(:),1); 
          node_elm(:,2) = node(nodes(:),2); 
          
          area = (abs( node_elm(1,1)-(node_elm(2,1)) ))...
          *(abs( node_elm(1,2)-(node_elm(4,2)) ));

          x = xp(1); y = xp(2);
              N=(1/area)*[ (x-node_elm(2,1))*(y-node_elm(4,2));
                          -(x-node_elm(1,1))*(y-node_elm(3,2));
                           (x-node_elm(4,1))*(y-node_elm(2,2));
                          -(x-node_elm(3,1))*(y-node_elm(1,2))];

              dNdX=(1/area)*[(y-node_elm(4,2)), (x-node_elm(2,1));
                            -(y-node_elm(3,2)),-(x-node_elm(1,1));
                             (y-node_elm(2,2)), (x-node_elm(4,1));
                            -(y-node_elm(1,2)),-(x-node_elm(3,1))];

           B(1,1:2:8) = dNdX(:,1)';
           B(2,2:2:8) = dNdX(:,2)';
           B(3,1:2:8) = dNdX(:,2)';
           B(3,2:2:8) = dNdX(:,1)';
           Nfem(1,1:2:8) = N(:,1)';
           Nfem(2,2:2:8) = N(:,1)';

end

end



