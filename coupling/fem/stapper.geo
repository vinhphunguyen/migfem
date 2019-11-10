l1=24;
l2=24;
c1=3;
c2=2;

h=0.5;

Point(1)={0,-c1,0,h};
Point(2)={l1,-c1,0,h};
Point(3)={l1, c1,0,h};
Point(4)={0, c1,0,h};

Point(5)={l1,-c2,0,0.5*h};
Point(6)={l1+l2,-c2,0,h};
Point(7)={l1+l2, c2,0,h};
Point(8)={l1, c2,0,0.5*h};

Line(1)={1,2};
Line(2)={2,5};
Line(3)={5,6};
Line(4)={6,7};
Line(5)={7,8};
Line(6)={8,3};
Line(7)={3,4};
Line(8)={4,1};
Line(9)={8,5};

Line Loop(1) = {1,2,-9,6,7,8};
Line Loop(2) = {3,4,5,9};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line(333)  = {8};
Physical Line(444)  = {4};

Physical Surface(888) = {1,2};


