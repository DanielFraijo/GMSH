// Gmsh project created on Sun Aug 18 16:02:40 2024
SetFactory("OpenCASCADE");


//+ All units in m
radius = 0.001;
flow_thick = 3 * radius;
radial_cells = 50;
axial_cells = 75;
ds = 1E-4; //+ First cell spacing. Currently not used but should be in future vrsn
gr = ((radius)/(ds))^(1/(radial_cells));

Point(1) = {0, 0, 0, 1.0};
Point(2) = {radius, radius, 0, 1.0};
Point(3) = {radius, 0, 0, 1.0};
Point(4) = {2*radius, 0, 0, 1.0};

Point(5) = {(-flow_thick+ radius), 0, 0, 1.0};
Point(6) = {radius, flow_thick, 0, 1.0};
Point(7) = {radius, 0, 0, 1.0};
Point(8) = {(radius + flow_thick), 0, 0, 1.0};


Circle(1) = {1, 3, 2}; 
Circle(2) = {2, 3, 4};
Circle(4) = {6, 3, 5};
Circle(5) = {6, 3, 8};


Line(6) = {1, 5};
Line(7) = {4, 8};
Line(8) = {2, 6};
//+
Curve Loop(1) = {4, -6, 1, -8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -7, -2, -8};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {6, 8, 7} = radial_cells Using Progression gr;
//+
Transfinite Curve {4, 1} = axial_cells Using Progression 1;
//+
Transfinite Curve {2, 5} = axial_cells Using Progression 1;
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Recombine Surface {1, 2};

Physical Curve("symmetry") = {6,7};
Physical Curve("farfield") = {4};
Physical Curve("outlet") = {5};
Physical Curve("wall") = {1,2};
Physical Surface("domain") = {1,4};
