h = 0.3;
h1 = 0.1;

r = 1;
l = 4; 

Point(1)={0,0,0,h1};
Point(2)={0,r,0,h1};
Point(3)={-r,0,0,h1};
Point(4)={-l,0,0,h};
Point(5)={-l,l,0,h};
Point(6)={0,l,0,h};

Circle(1)={2,1,3};
Line(2)={3,4};
Line(3)={4,5};
Line(4)={5,6};
Line(5)={6,2};

Line Loop(1)={1,2,3,4,5};
Plane Surface(1)={1};

Physical Surface(1)={1};

Physical Line(111)={2};
Physical Line(222)={5};
Physical Line(333)={3};
Physical Line(444)={4};
