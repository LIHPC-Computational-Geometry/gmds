// Gmsh project created on Sat Jul 27 11:56:00 2019
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 5, 0, 1.0};
//+
Point(4) = {0, 5, 0, 1.0};
//+
Point(5) = {2.5, 2.5, 0, 1.0};
//+
Point(6) = {2, 2.5, 0, 1.0};
//+
Point(7) = {3, 2.5, 0, 1.0};
//+
Point(8) = {3.5, 2.5, 0, 1.0};
//+
Point(9) = {1.5, 2.5, 0, 1.0};
//+
Delete {
  Point{6}; Point{7}; 
}
//+
Circle(1) = {9, 5, 8};
//+
Circle(2) = {8, 5, 9};
//+
Point(10) = {10, 1, 0, 1.0};
//+
Point(11) = {10, 4, 0, 1.0};
//+
Point(12) = {9, 5, 0, 1.0};
//+
Point(13) = {9, 4, 0, 1.0};
//+
Point(14) = {9, 0, 0, 1.0};
//+
Point(15) = {9, 1, 0, 1.0};
//+
Delete {
  Point{2}; Point{3}; 
}
//+
Circle(3) = {12, 13, 11};
//+
Circle(4) = {10, 15, 14};
//+
Line(5) = {4, 12};
//+
Line(6) = {11, 10};
//+
Line(7) = {14, 1};
//+
Line(8) = {4, 1};
//+
Line Loop(1) = {8, -7, -4, -6, -3, -5};
//+
Line Loop(2) = {2, 1};
//+
Plane Surface(1) = {1, 2};
