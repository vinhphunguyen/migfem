function [f] = biMatCrackBranchNode(r,theta,e,alpha)
% Compute the branch functions spanning the near tip field for bimaterial
% cracks.
% Inputs:
%   (r,theta) : polar coordinates of points where the branch
%               functions are to be evaluated
%   alpha     : inclination of the crack tip segment w.r.t x axis
% Vinh Phu Nguyen, nvinhphu@gmail.com
% PhD at Delft University of Technology, The Netherlands
% April, 2012, Ton Duc Thang University, Saigon, Vietnam
% Adapted from MXFEM of Mathiew Pais, Florida Univ.

sr       = sqrt(r); 
st       = sin(theta);
st2      = sin(theta/2);
ct2      = cos(theta/2);
lgr      = log(r);
coselogr = cos(e*lgr);
sinelogr = sin(e*lgr);
expeth   = exp(e*theta);
expmeth  = exp(-e*theta);


% 12 branch functions (checked)

f  = [sr*coselogr*expmeth*st2;         % Alpha 1 crack tip enrichment value
      sr*coselogr*expmeth*ct2;         % Alpha 2 crack tip enrichment value
      sr*coselogr*expeth*st2;          % Alpha 3 crack tip enrichment value
      sr*coselogr*expeth*ct2;          % Alpha 4 crack tip enrichment value
      sr*coselogr*expeth*st2*st;       % Alpha 5 crack tip enrichment value
      sr*coselogr*expeth*ct2*st;       % Alpha 6 crack tip enrichment value
      sr*sinelogr*expmeth*st2;         % Alpha 7 crack tip enrichment value
      sr*sinelogr*expmeth*ct2;         % Alpha 8 crack tip enrichment value
      sr*sinelogr*expeth*st2;          % Alpha 9 crack tip enrichment value
      sr*sinelogr*expeth*ct2;          % Alpha 10 crack tip enrichment value
      sr*sinelogr*expeth*st2*st;       % Alpha 11 crack tip enrichment value
      sr*sinelogr*expeth*ct2*st];

