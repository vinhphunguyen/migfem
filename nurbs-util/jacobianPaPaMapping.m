function j = jacobianPaPaMapping(rangeU,rangeV)
  J2xi    = 0.5 * ( rangeU(2) - rangeU(1) ); 
  J2eta   = 0.5 * ( rangeV(2) - rangeV(1) ); 
  j       = J2xi * J2eta;
end