function [Nv,dNdxi]=lagrange_basis(type,coord,dim)
  
% returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=lagrange_basis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%  
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
%   H4, H10, B8 and B27  
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University    
    
  if ( nargin == 2 )
    dim=1;
  end
  
  switch type
   case 'L2'  
    %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2
    %
    if size(coord,2) < 1
      disp('Error coordinate needed for the L2 element')
    else
     xi=coord(1);
     N=([1-xi,1+xi]/2)';
     dNdxi=[-1;1]/2;
    end
  
   case 'L3' 
    %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2----------3
    %
    if size(coord,2) < 1
      disp('Error two coordinates needed for the L3 element')
    else
     xi=coord(1);
     N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2];
     dNdxi=[xi-.5;xi+.5;-2*xi];
    end
    
   case 'T3'
    %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T3 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end
    
   case 'T3fs'
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T3fs element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end
        
   case 'T4'
    %%%%%%%%%% T4 FOUR NODE TRIANGULAR CUBIC BUBBLE ELEMENT %%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /      4     \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T4 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta-3*xi*eta;xi*(1-3*eta);eta*(1-3*xi);9*xi*eta];
      dNdxi=[-1-3*eta,-1-3*xi;
	     1-3*eta, -3*xi;
	     -3*eta,   1-3*xi;
	     9*eta,   9*xi ];
    end
    
   case 'T6'
    %%%%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         6          5
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1---------4----------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the T6 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                                  xi*(2*xi-1);
                                eta*(2*eta-1);
                              4*xi*(1-xi-eta);
                                     4*xi*eta;
                              4*eta*(1-xi-eta)];
        
      dNdxi=[4*(xi+eta)-3   4*(xi+eta)-3;
                   4*xi-1              0; 
                        0        4*eta-1;
           4*(1-eta-2*xi)          -4*xi;
                    4*eta           4*xi;
                   -4*eta  4*(1-xi-2*eta)];
    end
    
    
   case 'Q4'
    %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4--------------------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q4 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[ (1-xi)*(1-eta);
              (1+xi)*(1-eta);
              (1+xi)*(1+eta);
              (1-xi)*(1+eta)];
      dNdxi=1/4*[-(1-eta), -(1-xi);
		         1-eta,    -(1+xi);
		         1+eta,      1+xi;
                -(1+eta),   1-xi];
    end
    case 'Q8'
    %%%%%%%%%%%%%%% Q8 EIGHT NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8                    6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q8 element')
    else
        xi=coord(1); eta=coord(2);
        N=  [-0.25*(1-xi)*(1-eta)*(1+xi+eta);
             -0.25*(1+xi)*(1-eta)*(1-xi+eta);
             -0.25*(1+xi)*(1+eta)*(1-xi-eta);
             -0.25*(1-xi)*(1+eta)*(1+xi-eta);
              0.5*(1-xi^2)*(1-eta);
              0.5*(1+xi)*(1-eta^2);
              0.5*(1-xi^2)*(1+eta);
              0.5*(1-xi)*(1-eta^2)];
        
        % derivatives dN/dr
 dNdxi=[0.25*(1-eta)*(2*xi+eta), 0.25*(1-xi)*(xi+2*eta);
        0.25*(1-eta)*(2*xi-eta), -0.25*(1+xi)*(xi-2*eta);
        0.25*(1+eta)*(2*xi+eta), 0.25*(1+xi)*(xi+2*eta);
        0.25*(1+eta)*(2*xi-eta), -0.25*(1-xi)*(xi-2*eta);
        -(1-eta)*xi, -0.5*(1-xi^2);
        0.5*(1-eta^2), -(1+xi)*eta;
        -(1+eta)*xi, 0.5*(1-xi^2);
        -0.5*(1-eta^2), -(1-xi)*eta]; 
    end
    
    
   case 'Q9'
    %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8          9         6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q9 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[xi*eta*(xi-1)*(eta-1);
             xi*eta*(xi+1)*(eta-1);
             xi*eta*(xi+1)*(eta+1);
             xi*eta*(xi-1)*(eta+1);
            -2*eta*(xi+1)*(xi-1)*(eta-1);
            -2*xi*(xi+1)*(eta+1)*(eta-1);
            -2*eta*(xi+1)*(xi-1)*(eta+1);
            -2*xi*(xi-1)*(eta+1)*(eta-1);
             4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
      dNdxi=1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
                 eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
                 eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
                 eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
                -4*xi*eta*(eta-1),   -2*(xi+1)*(xi-1)*(2*eta-1);
         -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
                -4*xi*eta*(eta+1),   -2*(xi+1)*(xi-1)*(2*eta+1);
         -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
                 8*xi*(eta^2-1),      8*eta*(xi^2-1)];
    end
    
    case 'Q16'
    %%%%%%%%%%%%%%% Q16 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    10-----9------8------7
    %    |                    |
    %    |                    |
    %   11     16        15   6    
    %    |                    |
    %    |                    |
    %    |                    |
    %   12     13        14   5 
    %    |                    |
    %    |                    |
    %    1-----2--------3-----4
    %
    if size(coord,2) < 2
      disp('Error two coordinates needed for the Q16 element')
    else
      r=coord(1); s=coord(2);     

shape(1) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(1) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
         (r-1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(1)= -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));

shape(2) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
           (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(2) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
          (r-1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(2) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));  
     
shape(3) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(3) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
          (r+1.0/3.0)*(r-1.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(3) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*(-9.0/16.0)*...                    
        ((s+1.0/3.0)*(s-1.0/3.0)+ (s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 
shape(4)= 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            (-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhdr(4) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
     (r+1.0/3.0)*(r-1.0/3.0))*(-9.0/16.0)*(s+1.0/3.0)*(s-1.0/3.0)*(s-1.0);
dhds(4) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(-9.0/16.0)*...
         ((s+1.0/3.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 
shape(5) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);           
dhdr(5) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(5) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(27.0/16.0)*((s+1.0)*...
         (s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));

shape(6) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
           (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0); 
dhdr(6) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(6) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*(-27.0/16.0)*...
             ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));

shape(7) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(7) = 9.0/16.0*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0/3.0)+...
         (r+1.0/3.0)*(r-1.0/3.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(7) = 9.0/16.0*(r+1.0)*(r+1.0/3.0)*(r-1.0/3.0)*9.0/16.0*...
         ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));

shape(8) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(8) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r+1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(8) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*9.0/16.0* ((s+1.0)*...
             (s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
         
shape(9) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
            9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(9) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(9) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*9.0/16.0*((s+1.0)*...
                  (s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));

shape(10) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
           9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhdr(10) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*9.0/16.0*(s+1.0)*(s+1.0/3.0)*(s-1.0/3.0);
dhds(10)= -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*9.0/16.0*...
         ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0/3.0)+(s+1.0/3.0)*(s-1.0/3.0));
 
shape(11) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
            (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0); 
dhdr(11) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(11) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));

shape(12) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*...
           (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhdr(12) = -9.0/16.0*((r+1.0/3.0)*(r-1.0/3.0)+(r+1.0/3.0)*(r-1.0)+...
              (r-1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(12) = -9.0/16.0*(r+1.0/3.0)*(r-1.0/3.0)*(r-1.0)*(27.0/16.0)*...
                 ((s+1.0)*(s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));


shape(13) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
            (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0); 
dhdr(13) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
              (r-1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(13) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(27.0/16.0)*...
                 ((s+1.0)*(s-1.0/3.0)+(s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));

shape(14) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
           (27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhdr(14) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
              (r+1.0/3.0)*(r-1.0))*(27.0/16.0)*(s+1.0)*(s-1.0/3.0)*(s-1.0);
dhds(14) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)* (27.0/16.0)*...
                ((s+1.0)*(s-1.0/3.0)+ (s+1.0)*(s-1.0)+(s-1.0/3.0)*(s-1.0));
 

shape(15) = (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*...
            (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhdr(15) = (-27.0/16.0)*((r+1.0)*(r+1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r+1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(15) =  (-27.0/16.0)*(r+1.0)*(r+1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));

shape(16) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*...
          (-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhdr(16) = 27.0/16.0*((r+1.0)*(r-1.0/3.0)+(r+1.0)*(r-1.0)+...
             (r-1.0/3.0)*(r-1.0))*(-27.0/16.0)*(s+1.0)*(s+1.0/3.0)*(s-1.0);
dhds(16) = 27.0/16.0*(r+1.0)*(r-1.0/3.0)*(r-1.0)*(-27.0/16.0)*...
                 ((s+1.0)*(s+1.0/3.0)+(s+1.0)*(s-1.0)+(s+1.0/3.0)*(s-1.0));
N=shape(:);
dNdxi=[ dhdr(:) dhds(:)];
  end
   case 'H4'
    %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        /    |    \ 
    %       /     |     \
    %      1 -----|------3
    %         -   2  -
    if size(coord,2) < 3
      disp('Error three coordinates needed for the H4 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=[1-xi-eta-zeta;
                    xi;
                   eta;
                  zeta];
      dNdxi=[-1  -1  -1;
              1   0   0;
              0   1   0;
              0   0   1];
    end
    
   case 'H10'
    %%%%%%%%%%%%%%%% H10 TEN NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
    disp(['Element ',type,' not yet supported'])
    if size(coord,2) < 3
      disp('Error three coordinates needed for the H10 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=zeros(10,1);
      dNdxi=zeros(10,3);
    end
    
   case 'B8'
    %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
    % 
    %                  8 
    %               /    \    
    %            /          \
    %         /                \
    %      5                     \
    %      |\                     7
    %      |   \                / |
    %      |     \     4    /     |
    %      |        \    /        |
    %      |           6          |
    %      1           |          |
    %       \          |          3
    %          \       |        /
    %            \     |     /
    %               \  |  /
    %                  2
    %                
    if size(coord,2) < 3
      disp('Error three coordinates needed for the B8 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      I1=1/2-coord/2;
      I2=1/2+coord/2;
      N=[   I1(1)*I1(2)*I1(3);
            I2(1)*I1(2)*I1(3);
            I2(1)*I2(2)*I1(3);
            I1(1)*I2(2)*I1(3);
            I1(1)*I1(2)*I2(3);
            I2(1)*I1(2)*I2(3);
            I2(1)*I2(2)*I2(3);
            I1(1)*I2(2)*I2(3)   ];
      dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                 1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                 1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;      
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                 1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                 1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;
    end
    
   case 'B27'
    %%%%%%%%%%%%%% B27 TWENTY SEVEN NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%
    disp(['Element ',type,' not yet supported'])
    if size(coord,2) < 3
      disp('Error three coordinates needed for the B27 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=zeros(27,1);
      dNdxi=zeros(27,3);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   otherwise
    disp(['Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
  end
 
  I=eye(dim);
  Nv=[];
  for i=1:size(N,1)
    Nv=[Nv;I*N(i)];
  end
  
  if ( dim == 1 )
    B=dNdxi;
  elseif ( dim == 2 )
    B=zeros(dim*size(N,1),3);
    
    B(1:dim:dim*size(N,1)-1,1) = dNdxi(:,1);
    B(2:dim:dim*size(N,1),2)   = dNdxi(:,2);
    
    B(1:dim:dim*size(N,1)-1,3) = dNdxi(:,2);
    B(2:dim:dim*size(N,1),3)   = dNdxi(:,1);
  elseif ( dim == 3 )
    B=zeros(dim*size(N,1),6);
    
    disp('Error: need to add 3D N and dNdxi')
    
    B(1:dim:dim*size(N,1)-2,1) = dNdxi(:,1);
    B(2:dim:dim*size(N,1)-1,2) = dNdxi(:,2);
    B(3:dim:dim*size(N,1),3)   = dNdxi(:,3);
    
    B(2:dim:dim*size(N,1)-1,4) = dNdxi(:,3);
    B(3:dim:dim*size(N,1),4)   = dNdxi(:,2);
    
    B(3:dim:dim*size(N,1),5)   = dNdxi(:,1);
    B(1:dim:dim*size(N,1)-2,5) = dNdxi(:,3);
    
    B(1:dim:dim*size(N,1)-2,6) = dNdxi(:,2);
    B(2:dim:dim*size(N,1)-1,6) = dNdxi(:,1);
    
  end  
  
  % end of function