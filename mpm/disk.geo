h=5;
r1 = 30;
r2 = 40;

Point(1) = {0,0, 0, h};
Point(2) = {r1,0, 0, h};
Point(3) = {0, r1, 0, h};
Point(4) = {-r1,0, 0, h};
Point(5) = {0, -r1, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};

Point(20) = {r2,0, 0, h};
Point(30) = {0, r2, 0, h};
Point(40) = {-r2,0, 0, h};
Point(50) = {0, -r2, 0, h};
Circle(10) = {20, 1, 30};
Circle(20) = {30, 1, 40};
Circle(30) = {40, 1, 50};
Circle(40) = {50, 1, 20};
Line Loop(50) = {10, 20, 30, 40};
Plane Surface(60) = {5,50};
