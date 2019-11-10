  function [B, nbfu] = nrb_crv_basisfun__ (points, nrb);
%  __NRB_CRV_BASISFUN__: Undocumented internal function
%
%   Copyright (C) 2009 Carlo de Falco
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.
    n    = size (nrb.coefs, 2) -1;
    p    = nrb.order -1;
    u    = points;
    U    = nrb.knots;
    w    = nrb.coefs(4,:);
    
    spu  =  findspan (n, p, u, U); 
    nbfu =  numbasisfun (spu, u, p, U);
    
    N     = w(nbfu+1) .* basisfun (spu, u, p, U);
    B     = bsxfun (@(x,y) x./y, N, sum (N,2));

  end
