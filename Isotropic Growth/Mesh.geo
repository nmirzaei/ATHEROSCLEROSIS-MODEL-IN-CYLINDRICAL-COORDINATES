// Gmsh project created on Thu Dec 13 16:33:06 2018
SetFactory("OpenCASCADE");
//+
Point(1) = {0.126,0, 0, 1.0};
//+
Point(2) = {0.126, 3.36, 0, 1.0};
//+
Point(3) = {0.147, 0, 0, 1.0};
//+
Point(4) = {0.147, 3.36, 0, 1.0};
//+
Point(5) = {0.189, 0, 0, 1.0};
//+
Point(6) = {0.189, 3.36, 0, 1.0};
//+
Point(7) = {0.231, 0, 0, 1.0};
//+
Point(8) = {0.231, 3.36, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {2, 1};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 2};
//+
Line(5) = {4, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 3};
//+
Line(8) = {5, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {8, 6};
//+
Curve Loop(1) = {2, 1, 3, 4};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {3, 5, 6, 7};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {10, 6, 8, 9};
//+
Surface(3) = {5};
//+
Physical Curve(1) = {2};
//+
Physical Curve(2) = {3};
//+
Physical Curve(3) = {6};
//+
Physical Curve(4) = {9};
//+
Physical Curve(5) = {4, 5, 10};
//+
Physical Curve(6) = {1, 7, 8};
//+
Physical Surface(7) = {1};
//+
Physical Surface(8) = {2};
//+
Physical Surface(9) = {3};
//+

//+
Field[1] = MathEval;
//+
//we went to 247 with the following
//Field[1].F = "0.015*exp(0.4*(abs(x-0.3)+abs(y-4)))";
//we went to 299 with the following
//Field[1].F = "0.013*exp(0.4*(abs(x-0.3)+abs(y-4)))";
//The following works for the parameters (1,0.8,0.06) up to 3000 with order 2 elements but for (3.2,0.8,0.192) it only goes to 160 max 
Field[1].F = "0.00462*exp(0.9523809524*(abs(x-0.126)+abs(y-1.68)))";
//only to 158 for (3.2,0.8,0.192) via the following
//Field[1].F = "0.009*exp(0.4*(abs(x-0.3)+abs(y-4)))";
//Field[1].F = "0.025*exp(5.4*(abs(x-0.3)))";
Background Field = 1;
