a = 8;
r = 0.5;

cos45 = 0.70710678118654757;

h = 1;
n1 = 17;
n2 = 9;

// first one eight part

Point(1) = {0,0,0,h};
Point(2) = {r,0,0,h};
Point(3) = {0.5*a,0,0,h};
Point(4) = {0.5*a,0.5*a,0,h};
Point(5) = {cos45*r,cos45*r,0,h};

Line(1) = {2,3};
Line(2) = {3,4};
Line(3) = {4,5};

Circle(4) = {2,1,5};

Line Loop(1) = {1,2,3,-4};

Plane Surface(10) = {1};

Transfinite Line{1,3} = n1;
Transfinite Line{2,4} = n2;

Transfinite Surface{10} = {2,3,4,5};
Recombine Surface{10};

// second one-eight part

Point(6) = {0,0.5*a,0,h};
Point(7) = {0,r,0,h};

Line(5) = {6,7};
Circle(6) = {5,1,7};
Line(7) = {6,4};

Line Loop(2) = {6,-5,7,3};

Plane Surface(20) = {2};

Transfinite Line{6,7} = n2;
Transfinite Line{3,5} = n1;

Transfinite Surface{20} = {4,5,6,7};
Recombine Surface{20};

// third part = rotate second part
// 180 degree

Rotate {{0, 1, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{20}; }
}

Transfinite Line{22,24} = n2;
Transfinite Line{25,5} = n1;

Transfinite Surface{21} = {7,6,18,8};
Recombine Surface{21};

// fourth part = rotate third part

Rotate {{-0.70710678118654757, 0.70710678118654757, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{21}; }
}

Transfinite Line{29,27} = n2;
Transfinite Line{28,25} = n1;

Transfinite Surface{26} = {25,21,8,18};
Recombine Surface{26};

// fifth part = rotate fourth part

Rotate {{1, 0, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{26}; }
}

Transfinite Line{31,33} = n2;
Transfinite Line{34,28} = n1;

Transfinite Surface{30} = {36,26,21,25};
Recombine Surface{30};

// sixth part = rotate fifth part

Rotate {{0.70710678118654757, 0.70710678118654757, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{30}; }
}

Transfinite Line{36,38} = n2;
Transfinite Line{34,37} = n1;

Transfinite Surface{35} = {36,43,39,26};
Recombine Surface{35};

// seventh part = rotate sixth part

Rotate {{0, 1, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{35}; }
}

Transfinite Line{40,42} = n2;
Transfinite Line{37,43} = n1;

Transfinite Surface{39} = {43,54,44,39};
Recombine Surface{39};

// final part = rotate seventh part

Rotate {{-0.70710678118654757, 0.70710678118654757, 0}, {0, 0, 0}, Pi} {
  Duplicata { Surface{39}; }
}

Transfinite Line{47,45} = n2;
Transfinite Line{43,1} = n1;

Transfinite Surface{44} = {44,54,3,2};
Recombine Surface{44};

Physical Surface(100) = {10,20,21,26,30,39,44,35};

// left, lower, right, upper edges
// opposite edges = same direction

Physical Line(1) = {-33,29};
Physical Line(2) = {-38,-42};
Physical Line(3) = {-47,2};
Physical Line(4) = {-24,7};

// four corner nodes
// left bottom, anti-clockwise

Physical Point(10) = {36};
Physical Point(20) = {54};
Physical Point(30) = {4};
Physical Point(40) = {18};


