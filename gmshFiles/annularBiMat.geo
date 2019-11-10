h=0.05;

cos45=0.707106781186548;

r1 = 1;
r2 = 2;

Point(1)={0,0,0,h};
Point(2)={0,r1,0,h};
Point(3)={-r1,0,0,h};
Point(4)={-r2,0,0,h};
Point(5)={0,r2,0,h};
Point(6)={-cos45,cos45,0,h};
Point(7)={-cos45*2,cos45*2,0,h};

Circle(1)={2,1,6};
Circle(2)={6,1,3};
Circle(3)={5,1,7};
Circle(4)={7,1,4};
Line(5)={3,4};
Line(6)={5,2};
Line(7)={6,7};

Line Loop(1)={1,7,-3,6};
Line Loop(2)={2,5,-4,-7};

Plane Surface(1)={1};
Plane Surface(2)={2};

Physical Surface(888)={2};
Physical Surface(999)={1};
Physical Line(111)={5};
Physical Line(222)={6};

