h=0.02;

r1 = 0.3;
r2 = 0.5;

Point(1)={0,0,0,h};
Point(2)={r1,0,0,h};
Point(3)={0,r1,0,h};
Point(4)={r2,0,0,h};
Point(5)={0,r2,0,h};

Circle(1)={2,1,3};
Circle(2)={4,1,5};
Line(3)={2,4};
Line(4)={3,5};

Line Loop(1)={3,2,-4,-1};

Plane Surface(1)={1};

Physical Surface(888)={1};
Physical Line(111)={3};
Physical Line(222)={4};
Physical Line(333)={1};

