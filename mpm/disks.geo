h=0.1;
r=0.2;
l=1;
d=0.05;

Point(1) = {r+d, r+d, 0, h};
Point(2) = {2*r+d, r+d, 0, h};
Point(3) = {r+d, 2*r+d, 0, h};
Point(4) = {0+d, r+d, 0, h};
Point(5) = {r+d, 0+d, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Point(10) = {l-r-d, l-r-d, 0, h};
Point(20) = {l-d, l-r-d, 0, h};
Point(30) = {l-r-d, l-d, 0, h};
Point(40) = {l-2*r-d, l-r-d, 0, h};
Point(50) = {l-r-d, l-2*r-d, 0, h};
Circle(10) = {20, 10, 30};
Circle(20) = {30, 10, 40};
Circle(30) = {40, 10, 50};
Circle(40) = {50, 10, 20};
Line Loop(50) = {10, 20, 30, 40};
Plane Surface(60) = {50};
