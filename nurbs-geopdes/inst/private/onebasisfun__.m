function N = onebasisfun__ (u, p, U)

%  __ONEBASISFUN__: Undocumented internal function
%
%   Copyright (C) 2009 Carlo de Falco
%   This software comes with ABSOLUTELY NO WARRANTY; see the file
%   COPYING for details.  This is free software, and you are welcome
%   to distribute it under the conditions laid out in COPYING.

  N = 0;
  if (~ any (U <= u)) || (~ any (U > u))
    return;
  elseif (p == 0)
      N = 1;
      return;
  end

 ln = u - U(1);
 ld = U(end-1) - U(1);
 if (ld ~= 0)
   N = N + ln * onebasisfun__ (u, p-1, U(1:end-1))/ ld; 
 end

 dn = U(end) - u;
 dd = U(end) - U(2);
 if (dd ~= 0)
   N = N + dn * onebasisfun__ (u, p-1, U(2:end))/ dd;
 end
  
end
