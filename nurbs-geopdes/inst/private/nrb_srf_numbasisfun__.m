function idx = nrb_srf_numbasisfun__ (points, nrb)

%  __NRB_SRF_NUMBASISFUN__: Undocumented internal function
%
%   Copyright (C) 2009 Carlo de Falco
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.
  
  m   = nrb.number(1)-1;
  n   = nrb.number(2)-1;
  
  npt = size(points,2);
  u   = points(1,:);
  v   = points(2,:);

  U   = nrb.knots{1};
  V   = nrb.knots{2};

  p   = nrb.order(1)-1;
  q   = nrb.order(2)-1;

  spu = findspan (m, p, u, U); 
  Ik  = numbasisfun (spu, u, p, U);

  spv = findspan (n, q, v, V);
  Jk  = numbasisfun (spv, v, q, V);
  
  for k=1:npt
    [Jkb, Ika] = meshgrid(Jk(k, :), Ik(k, :)); 
    idx(k, :)  = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
  end

end
