function B = getBmatrix3D(nn,dRdx)
% Compute the 3D strain-displacement matrix
% using the order [epsilon_xx epsilon_yy epsilon_zz epsilon_xy epsilon_yz
% epsilon_zx] with U =[u_x^1
%                      u_y^1 
%                      u_z^1 
%                      u_x^n 
%                      u_y^1 u_z^n]
% Note that the way displacements are stored was changed.
% In former versions, U = [u_x1 u_x2 ... u_y1 uy_2 .. uz_1 uz_2...]
%
% Vinh Phu Nguyen
% Delft University of Technology

nn3 = 3 * nn;

B(1,1:3:nn3)    = dRdx(:,1)';
B(2,2:3:nn3)    = dRdx(:,2)';
B(3,3:3:nn3)    = dRdx(:,3)';

B(4,1:3:nn3)    = dRdx(:,2)';
B(4,2:3:nn3)    = dRdx(:,1)';

B(5,3:3:nn3)    = dRdx(:,2)';
B(5,2:3:nn3)    = dRdx(:,3)';

B(6,1:3:nn3)    = dRdx(:,3)';
B(6,3:3:nn3)    = dRdx(:,1)';


