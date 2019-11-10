l=1;
w=10;

h=1;

n1 = 30;
n2 = 30;
n3 = 30;

Point(1)={0,8.5,0,h};
Point(2)={l,8.5,0,h};
Point(3)={l,w-0.5*l,0,h};
Point(4)={2,w-0.5*l,0,h};
Point(5)={2,w+0.5*l,0,h};
Point(6)={0,w+0.5*l,0,h};
Point(7)={0,w-0.5*l,0,h};
Point(8)={l,w+0.5*l,0,h};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,5};
Line(5)={5,8};
Line(6)={8,6};
Line(7)={6,7};
Line(8)={7,1};
Line(9)={7,3};
Line(10)={3,8};

Line Loop(1) = {1,2,-9,8};
Line Loop(2) = {3,4,5,-10};
Line Loop(3) = {9,10,6,7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Transfinite Line{1,9,4,6,7,10} = n1;
Transfinite Line{2,8}          = n2;
Transfinite Line{3,5}          = n3;

Transfinite Surface{1} = {1,2,3,7};
Transfinite Surface{3} = {7,3,8,6};
Transfinite Surface{2} = {3,4,5,8};

Recombine Surface{1,2,3};

Physical Point(111) = {4}; // positive force
Physical Point(222) = {5}; // negative force

Physical Line(333)  = {1};
Physical Line(444)  = {4};

Physical Surface(888) = {1,2,3};


