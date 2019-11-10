function [N,dNdX] = Shape_Function_Global_GMPM(nodes,xp,dim,dx,dy)

global node



lpx = dim(1)/2;
lpy = dim(2)/2;
num_nodes = length(nodes);
N    = zeros(num_nodes,1);
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



