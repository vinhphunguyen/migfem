h=3.0;
a=50;

Point(1)={-a,a,0,h};
Point(2)={-a,-a,0,h};
Point(3)={a,-a,0,h};
Point(4)={a,0,0,h};
Point(5)={0,0,0,0.02*h};
Point(6)={0,a,0,h};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,6};
Line(6)={6,1};

Line Loop(1)={1,2,3,4,5,6};

Plane Surface(1)={1};

Physical Surface(100)={1};

Physical Line(111)={3};
Physical Line(222)={6};
Physical Line(333)={1};
