function [B, N] = nrb_srf_basisfun__ (points, nrb);

%  __NRB_SRF_BASISFUN__: Undocumented internal function
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

    spu    = findspan (m, p, u, U); 
    spv    = findspan (n, q, v, V);
    NuIkuk = basisfun (spu, u, p, U);
    NvJkvk = basisfun (spv, v, q, V);

    indIkJk = nrbnumbasisfun (points, nrb);

    for k=1:npt
      wIkaJkb(1:p+1, 1:q+1) = reshape (w(indIkJk(k, :)), p+1, q+1); 
      NuIkukaNvJkvk(1:p+1, 1:q+1) = (NuIkuk(k, :).' * NvJkvk(k, :));
      RIkJk(k, :) = reshape((NuIkukaNvJkvk .* wIkaJkb ./ sum(sum(NuIkukaNvJkvk .* wIkaJkb))),1,[]);
    end
    
    B = RIkJk;
    N = indIkJk;
    
  end