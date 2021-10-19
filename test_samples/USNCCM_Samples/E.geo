// Gmsh project created on Sat Jul 27 10:09:49 2019
//+
Point(1) = {4, 0, 0, 1.0};
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Point(5) = {4, 2, 0, 1.0};
//+
Point(6) = {4, 1, 0, 1.0};
//+
Line(1) = {4, 5};
//+
Line(2) = {1, 2};
//+
Circle(4) = {1, 6, 5};//+
Circle(5) = {4, 3, 2};//+
Point(7) = {2, 1, 0, 1.0};
//+
Point(8) = {1.5, 1, 0, 1.0};
//+
Point(9) = {2.5, 1, 0, 1.0};
//+
Point(10) = {4.5, 1, 0, 1.0};
//+
Point(11) = {3.5, 1, 0, 1.0};
//+
Point(12) = {0.5, 1, 0, 1.0};
//+
Point(13) = {-0.5, 1, 0, 1.0};
//+
Circle(6) = {13, 3, 12};
//+
Circle(7) = {12, 3, 13};
//+
Circle(8) = {8, 7, 9};
//+
Circle(9) = {9, 7, 8};
//+
Circle(10) = {11, 6, 10};
//+
Circle(11) = {10, 6, 11};
//+
Line Loop(1) = {5, -2, 4, -1};
//+
Line Loop(2) = {7, 6};
//+
Line Loop(3) = {8, 9};
//+
Line Loop(4) = {10, 11};
//+
Plane Surface(1) = {1, 2, 3, 4};
//+
Delete {
  Surface{1}; 
}
//+
Delete {
  Line{11}; 
}
//+
Delete {
  Line{10}; 
}
//+
Delete {
  Point{11}; 
}
//+
Delete {
  Point{10}; 
}
//+
Plane Surface(1) = {1, 2, 3};
