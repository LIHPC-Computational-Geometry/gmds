// Gmsh project created on Mon Jul 29 07:11:13 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 5, 0, 1.0};
//+
Point(4) = {0, 5, 0, 1.0};
//+
Circle(1) = {2.2, 1.3, 0, 0.5, 0, 2*Pi};
//+
Circle(2) = {2.8, 3.4, 0, 0.5, 0, 2*Pi};
//+
Circle(3) = {5, 1.2, 0, 0.5, 0, 2*Pi};
//+
Circle(4) = {5.7, 4, 0, 0.5, 0, 2*Pi};
//+
Circle(5) = {7.1, 2.7, 0, 0.5, 0, 2*Pi};
//+
Circle(6) = {8.4, 1, 0, 0.5, 0, 2*Pi};
//+
Circle(7) = {8.7, 4.1, 0, 0.5, 0, 2*Pi};
//+
Circle(8) = {1, 4, 0, 0.5, 0, 2*Pi};
//+
Line(9) = {4, 3};
//+
Line(10) = {3, 2};
//+
Line(11) = {2, 1};
//+
Line(12) = {1, 4};
//+
Line Loop(1) = {9, 10, 11, 12};
//+
Line Loop(2) = {8};
//+
Line Loop(3) = {2};
//+
Line Loop(4) = {1};
//+
Line Loop(5) = {3};
//+
Line Loop(6) = {5};
//+
Line Loop(7) = {4};
//+
Line Loop(8) = {7};
//+
Line Loop(9) = {6};
//+
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
