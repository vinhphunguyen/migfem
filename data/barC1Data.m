knotVec    = [0 0 0 1 1 1];
controlPts = [0 0;
              0.5 0;
	          1 0];
p       = 2;
noGPs   = 3;

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));