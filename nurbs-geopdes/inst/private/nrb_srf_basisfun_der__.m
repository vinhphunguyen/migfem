function [Bu, Bv, N] = nrb_srf_basisfun_der__ (points, nrb);

%  __NRB_SRF_BASISFUN_DER__: Undocumented internal function
%
%   Copyright (C) 2009 Carlo de Falco
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  m    = size (nrb.coefs, 2) -1;
  n    = size (nrb.coefs, 3) -1;
  
  p    = nrb.order(1) -1;
  q    = nrb.order(2) -1;
  
  u = points(1,:);
  v = points(2,:);
  npt = length(u);
  
  U    = nrb.knots{1};
  V    = nrb.knots{2};
  
  w    = squeeze(nrb.coefs(4,:,:));
  
  spu  =  findspan (m, p, u, U); 
  spv  =  findspan (n, q, v, V);
  N    =  nrbnumbasisfun (points, nrb);

  NuIkuk = basisfun (spu, u, p, U); 
  NvJkvk = basisfun (spv, v, q, V);
  
  NuIkukprime = basisfunder (spu, p, u, U, 1);
  NuIkukprime = reshape (NuIkukprime(:,2,:), npt, []);
  
  NvJkvkprime = basisfunder (spv, q, v, V, 1);
  NvJkvkprime = reshape (NvJkvkprime(:,2,:), npt, []);
  
  for k=1:npt
    wIkaJkb(1:p+1, 1:q+1) = reshape (w(N(k, :)), p+1, q+1);
    
    Num    = (NuIkuk(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
    Num_du = (NuIkukprime(k, :).' * NvJkvk(k, :)) .* wIkaJkb;
    Num_dv = (NuIkuk(k, :).' * NvJkvkprime(k, :)) .* wIkaJkb;
    Denom  = sum(sum(Num));
    Denom_du = sum(sum(Num_du));
    Denom_dv = sum(sum(Num_dv));
    
    Bu(k, :) = reshape((Num_du/Denom - Denom_du.*Num/Denom.^2),1,[]);
    Bv(k, :) = reshape((Num_dv/Denom - Denom_dv.*Num/Denom.^2),1,[]);
  end
  
end