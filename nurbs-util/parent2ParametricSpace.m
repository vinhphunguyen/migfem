function xi = parent2ParametricSpace(range,xibar)
  xi = 0.5 * ( ( range(2) - range(1) ) * xibar + range(2) + range(1)); 
end