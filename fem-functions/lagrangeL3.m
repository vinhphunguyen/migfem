function [N dNdxi]=lagrangeL3(type,coord)

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
        xi   =coord(1);
        N    =[(1-xi)*xi/(-2);1-xi^2;(1+xi)*xi/2];
        dNdxi=[xi-.5;-2*xi;xi+.5];
end