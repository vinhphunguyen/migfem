function j = jacobianPaPaMapping3d(rangeU,rangeV,rangeW)
  J2xi    = 0.5 * ( rangeU(2) - rangeU(1) ); 
  J2eta   = 0.5 * ( rangeV(2) - rangeV(1) ); 
  J2zeta  = 0.5 * ( rangeW(2) - rangeW(1) ); 
  j       = J2xi * J2eta * J2zeta;
end