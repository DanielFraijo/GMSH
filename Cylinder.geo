// Gmsh project updated for cylinder with splines
SetFactory("OpenCASCADE");

// All units in meters
radius = 0.001;        // Radius of the geometry
flow_thick = 3 * radius; // Thickness of the flow domain around the geometry
rear_length = radius * 2;
front_length = radius * 0.3;
ogive_length = front_length + radius; // Total length of the ogive
ogive_radius = 0.002;
points = 100;
C = 0.5; // Haack constant for the ogive shape (adjust as needed)

// Mesh Parameters
axial = 200; // Axial cells
ds = 1E-6;
gr = 1.01;
radial = Ceil(1 + Log(1 - 0.005 * (1 - gr)/ds)/Log(gr));

// Define Points for Geometry
Point(1) = {0, 0, 0, 1.0};
Point(2) = {radius, radius, 0, 1.0};
Point(3) = {radius, 0, 0, 1.0};
Point(4) = {2*radius, 0, 0, 1.0};
Point(5) = {-front_length, 0, 0, 1.0};
Point(6) = {2 * radius + rear_length, 0, 0, 1.0};
Point(7) = {front_length + radius * 2 + rear_length, ogive_radius, 0, 1.0};


// Generate Points for Farfield Ogive with Haack Series
For j In {0:points-1}
    x = j * ogive_length / (points - 1);  // Distribute points evenly from 0 to ogive_length
    Y = (ogive_radius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / ogive_length) - (Sin(2 * Acos(1 - (2 * x) / ogive_length)) / 2) + C * Sin(Acos(1 - (2 * x) / ogive_length))^3);

// Define a point with calculated coordinates, aligning with the lengths
Point(8 + j) = {x - front_length, Y, 0, 1.0};
EndFor

// Define Geometrical Shapes
Circle(1) = {1, 3, 2}; // Circular arc from Point 1 to 2 via 3
Circle(2) = {2, 3, 4}; // Circular arc from Point 2 to 4 via 3
Line(3) = {2, points + 7}; // Vertical line from cylinder to midpoint of ogive
Line(4) = {1, 5};  // Line from center to front of domain
Line(5) = {4, 6};  // Line from end of circle to rear_length
Spline(6) = {5, 9:(points + 7)};
//BSpline(7) = {points + 7, 8, 7, 9 , 6};
BSpline(7) = {points + 7, 7, 6};

//+
Transfinite Curve {6, 1} = 15 Using Progression 1;
//+
Transfinite Curve {7, 2} = 15 Using Progression 1;
//+
Transfinite Curve {4, 3, 5} = 10 Using Progression 1;
//+
Curve Loop(1) = {6, -3, -1, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -5, -2, 3};
//+
Plane Surface(2) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Recombine Surface {1, 2};
