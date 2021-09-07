// Gmsh project created on Mon Jul 29 08:56:47 2019
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 5, 0, 1.0};
//+
Point(3) = {1, 5, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Point(5) = {2, 5, 0, 1.0};
//+
Point(6) = {5, 2, 0, 1.0};
//+
Point(7) = {5, 3, 0, 1.0};
//+
Point(8) = {5, 0, 0, 1.0};
//+
Point(9) = {7, 0, 0, 1.0};
//+
Point(10) = {8, 0, 0, 1.0};
//+
Point(11) = {8, -2, 0, 1.0};
//+
Point(12) = {7, -2, 0, 1.0};
//+
Point(13) = {7, -3, 0, 1.0};
//+
Point(14) = {6, -3, 0, 1.0};
//+
Point(15) = {-2, -3, 0, 1.0};
//+
Point(16) = {6, -2, 0, 1.0};
//+
Recursive Delete {
  Point{15}; Point{5}; 
}
//+
Line(1) = {1, 2};
//+
Line(2) = {12, 9};
//+
Line(3) = {10, 11};
//+
Line(4) = {12, 16};
//+
Line(5) = {14, 16};
//+
Line(6) = {13, 14};
//+
Circle(7) = {6, 8, 9};
//+
Circle(8) = {7, 8, 10};
//+
Circle(9) = {11, 12, 13};
//+
BSpline(10) = {2, 3, 7};
//+
Delete {
  Line{10}; 
}
//+
Point(17) = {4, 3, 0, 1.0};
//+
Point(18) = {4, 2, 0, 1.0};
//+
BSpline(10) = {2, 3, 17, 7};
//+
BSpline(11) = {1, 4, 18, 18, 6};
//+
Delete {
  Line{11}; 
}
//+
BSpline(11) = {1, 4, 18, 6};
//+
Line Loop(1) = {1, 10, 8, 3, 9, 6, 5, -4, 2, -7, -11};
//+
Plane Surface(1) = {1};
//+
Line(12) = {7, 6};
//+
Delete {
  Surface{1}; 
}
//+
Delete {
  Line{8}; Line{7}; Point{8}; 
}
//+
Delete {
  Line{2}; 
}
//+
Delete {
  Point{9}; 
}
//+
Recursive Delete {
  Line{4}; 
}
//+
Recursive Delete {
  Line{5}; 
}
//+
Recursive Delete {
  Point{10}; 
}
//+
Recursive Delete {
  Point{12}; 
}
//+
Recursive Delete {
  Point{10}; 
}
//+
Recursive Delete {
  Point{10}; Line{3}; 
}
//+
Recursive Delete {
  Line{9}; 
}
//+
Recursive Delete {
  Line{6}; 
}
//+
Line Loop(2) = {1, 10, 12, -11};
//+
Plane Surface(1) = {2};
//+
Delete {
  Surface{1}; 
}
//+
Point(19) = {1.5, 2.5, 0, 1.0};
//+
Point(20) = {1.5, 1.5, 0, 1.0};
//+
Point(21) = {1.5, 3, 0, 1.0};
//+
Point(22) = {1.5, 3.5, 0, 1.0};
//+
Delete {
  Point{21}; 
}
//+
Ellipse(13) = {22, 19, 20, 20};
//+
Delete {
  Line{13}; 
}
//+
Point(23) = {1, 2.5, 0, 1.0};
//+
Point(24) = {2, 2.5, 0, 1.0};
//+
Ellipse(13) = {24, 19, 22, 22};
//+
Ellipse(14) = {22, 19, 23, 23};
//+
Ellipse(15) = {23, 19, 20, 20};
//+
Ellipse(16) = {20, 19, 24, 24};
//+
Line Loop(3) = {14, 15, 16, 13};
//+
Plane Surface(1) = {2, 3};
//+
Delete {
  Surface{1}; 
}
//+
Delete {
  Line{14}; 
}
//+
Delete {
  Point{23}; 
}
//+
Delete {
  Line{15}; 
}
//+
Delete {
  Point{23}; 
}
//+
Delete {
  Line{1}; 
}
//+
Point(25) = {0, 3, 0, 1.0};
//+
Point(26) = {0, 2, 0, 1.0};
//+
Point(27) = {0.5, 2, 0, 1.0};
//+
Point(28) = {0.5, 3, 0, 1.0};
//+
Point(29) = {1, 3.5, 0, 1.0};
//+
Point(30) = {1, 1.5, 0, 1.0};
//+
BSpline(17) = {25, 28, 29, 22};
//+
BSpline(18) = {20, 30, 27, 26};
//+
Line(19) = {2, 25};
//+
Line(20) = {26, 1};
//+
Line Loop(4) = {17, -13, -16, 18, 20, 11, -12, -10, 19};
//+
Plane Surface(1) = {4};
