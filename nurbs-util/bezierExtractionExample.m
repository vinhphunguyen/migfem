% Example showing how to the Bezier extraction 
%

uknot = [0 0 0 1/3 2/3 1 1 1];
vknot = [0 0 0 1/3 2/3 1 1 1];
p     = 2;
q     = 2;

[C1,nb1] = bezierExtraction(uknot,p);
[C2,nb2] = bezierExtraction(vknot,q);

C        = bezierExtraction2D(uknot,vknot,p,q);