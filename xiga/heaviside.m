% *************************************************************************
%                 TWO DIMENSIONAL ELEMENT FREE GALERKIN CODE
%                            Nguyen Vinh Phu
%                        LTDS, ENISE, Juillet 2006
% *************************************************************************
function [f] = heaviside(distance)
% Compute the Heaviside function
% Input : distance is the signed distance to the crack line

f = sign(distance) ;

