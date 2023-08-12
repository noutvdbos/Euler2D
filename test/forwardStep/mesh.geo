D = 0.2;

inlet_h = 49;
inlet_v = 81;
outlet_h = 193;
outlet_v = 65;
step_v = 17;

Point(1) = {3*D, 0, 0};
Point(2) = {0, 0, 0};
Point(3) = {0, 5*D, 0};
Point(4) = {3*D, 5*D,0};
Point(5) =  {15*D,5*D,0};
Point(6) =  {15*D,D,0};
Point(7) = {3*D ,D, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 4};
Line(9) = {7, 1};
//Line Loop(11) = {1,2,5,6,7,9,4};
Line Loop(11) = {4,1,2,3};
Line Loop(12) = {5,6,7,8};
Plane Surface(11) = {11};
Plane Surface(12) = {12};

//+
Transfinite Curve {1,3} = inlet_v Using Progression 1;
//+
Transfinite Curve {2, 4} = inlet_h Using Progression 1;
//+
Transfinite Curve {5, 7} = outlet_h Using Progression 1;
//+
Transfinite Curve {6, 8} = outlet_v Using Progression 1;
//+
Transfinite Curve {9} = step_v Using Progression 1;
//+

Transfinite Surface "*";
//Recombine Surface "*";

Transfinite Volume "*";
//Recombine Volume "*";



//+
Physical Curve("Inlet") = {1};
//+
Physical Curve("top") = {2,5};
//+
Physical Curve("outlet") = {6};
//+
Physical Curve("bottom") = {4,9};
//+
Physical Curve("longbottom") = {7};
//+
Physical Surface("domain") = {11,12};
//+

Mesh 2;  // Generate 2D mesh
Coherence Mesh;  // Remove duplicate entities


