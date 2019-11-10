knotVec    = [0 0 0 0 1 1 1 1];
controlPts = [0 0;
              0.3 0;
              0.5 0;
	          1 0];
          
p       = 3;  
noGPs   = 4; 

knotVec = knotVec/max(knotVec);
weights = ones(1,size(controlPts,1));