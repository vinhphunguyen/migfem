l=1;
w=10;

h=1;

n1 = 30;
n2 = 30;
n3 = 60;
n4 = 100;

Point(1)={0,0,0,h};
Point(2)={l,0,0,h};
Point(3)={l,w-0.5*l,0,h};
Point(4)={0.5*(w+l),w-0.5*l,0,h};
Point(5)={0.5*(w+l),w+0.5*l,0,h};
Point(6)={0,w+0.5*l,0,h};
Point(7)={0,w-0.5*l,0,h};
Point(8)={l,w+0.5*l,0,h};
Point(9)={0,8.5,0,h};
Point(10)={l,8.5,0,h};
Point(11)={l*2,w-0.5*l,0,h};
Point(12)={l*2,w+0.5*l,0,h};

Line(1)={1,2};
Line(2)={2,10};
Line(3)={10,3};
Line(4)={3,11};
Line(5)={11,4};
Line(6)={4,5};
Line(7)={5,12};
Line(8)={12,8};
Line(9)={8,6};
Line(10)={6,7};
Line(11)={7,9};
Line(12)={9,1};
Line(13)={9,10};
Line(14)={7,3};
Line(15)={3,8};
Line(16)={11,12};

Line Loop(1) = {1,2,-13,12};
Line Loop(2) = {13,3,-14,11};
Line Loop(3) = {14,15,9,10};
Line Loop(4) = {4,16,8,-15};
Line Loop(5) = {5,6,7,-16};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Transfinite Line{1,13,14,9,4,8} = n1;
Transfinite Line{10,15,16,6,11,3} = n2;
Transfinite Line{7,5}          = n3;
Transfinite Line{2,12}          = n4;

Transfinite Surface{1} = {1,2,10,9};
Transfinite Surface{2} = {9,10,3,7};
Transfinite Surface{3} = {7,3,8,6};
Transfinite Surface{4} = {3,11,12,8};
Transfinite Surface{5} = {11,4,5,12};

Recombine Surface{1,2,3,4,5};

Physical Point(111) = {4}; // positive force
Physical Point(222) = {5}; // negative force

Physical Line(333)  = {1};
Physical Line(444)  = {6};

Physical Surface(888) = {1,2,3,4,5};


