  function [Bu, nbfu] = nrb_crv_basisfun_der__ (points, nrb);
%  __NRB_CRV_BASISFUN_DER__: Undocumented internal function
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
    
    N     = basisfun (spu, u, p, U);

    Nprime = basisfunder (spu, p, u, U, 1);
    Nprime = squeeze(Nprime(:,2,:));

    
    [Dpc, Dpk]  = bspderiv (p, w, U);
    D           = bspeval  (p, w, U, u);
    Dprime      = bspeval  (p-1, Dpc, Dpk, u);
    

    Bu1   = bsxfun (@(np, d) np/d , Nprime.', D);
    Bu2   = bsxfun (@(n, dp)  n*dp, N.', Dprime./D.^2);
    Bu    = w(nbfu+1) .* (Bu1 - Bu2).';

  end
