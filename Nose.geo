// Original Script
// Parameters for geometry
radius = 0.01;
theta = 15 * Pi / 180;
length = 0.04;
flow_radius = 3.5 * radius;
flow_length = flow_radius + length;

radial_cells = 150;
axial_cells = 200;
ds = 1E-6;
gr=((radius)/(ds))^(1/(radial_cells));

// Points for geometry
Point(1) = {0, 0, 0};
Point(2) = {radius, 0, 0};
Point(3) = {radius - radius * Sin(theta), radius * Cos(theta), 0};
Point(4) = {length, 0, 0};
Point(5) = {length, (length - (radius - radius * Sin(theta))) * Tan(theta) + radius * Cos(theta), 0};

// Geometrical shapes
Circle(2) = {1, 2, 3};
Line(3) = {3, 5};

// Additional points for flow geometry
Point(6) = {-flow_radius + radius, 0, 0};
Point(8) = {length, (flow_length - (flow_radius - flow_radius * Sin(theta))) * Tan(theta) + flow_radius * Cos(theta), 0};

// Flow lines
Line(5) = {1, 6};
Line(6) = {5, 8};

// Additional points for flow
Point(10) = {-(flow_radius * Sin(theta)) + radius, flow_radius * Cos(theta), 0};

// More lines and circular arcs
Line(7) = {3, 10};
Circle(8) = {6, 2, 10};
Line(9) = {8, 10};
//Line(10) ={2, 10};

// Surface definitions for geometry

// Physical entities for boundary conditions
Physical Curve("farfield", 10) = {8, 9};
Physical Curve("symmetry", 11) = {5};
Physical Curve("outlet", 12) = {6};
Physical Curve("wall", 13) = {2, 3};


//+
Curve Loop(1) = {9, -7, 3, 6};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {7, 6} = radial_cells Using Progression gr;
//+
Transfinite Curve {9, 3} = axial_cells Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Transfinite Curve {8, 2} = axial_cells Using Progression 1;
//+
Transfinite Curve {5} = radial_cells Using Progression gr;
//+
Curve Loop(2) = {8, -7, -2, 5};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {2};
//+
Recombine Surface {2};
