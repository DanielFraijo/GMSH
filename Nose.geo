// Original Script

// Geometry Parameters
radius = 0.01; // Radius for the wall
theta = 15 * Pi / 180; // Angle in radians
length = 0.04; // Wall length
Far_R = 5 * radius; // Farfield radius
Far_L = length * 1.3; // Farfield length
points = 1000; // Number of points for farfield approximation
midpoint = (points * 0.2) + 5;

// Mesh Parameters
axial = 200; // Axial cells
ds = 1E-6;
gr = 1.01;
radial = Ceil(1 + Log(1 - 0.005 * (1 - gr)/ds)/Log(gr));

// Define Points for Geometry
Point(1) = {0, 0, 0};
Point(2) = {radius, 0, 0};
Point(3) = {radius - radius * Sin(theta), radius * Cos(theta), 0};
Point(4) = {length, 0, 0};
Point(5) = {length, (length - (radius - radius * Sin(theta))) * Tan(theta) + radius * Cos(theta), 0};

// Define Geometrical Shapes
Circle(2) = {1, 2, 3}; // Creates a circular arc
Line(3) = {3, 5};       // Connects the end of the arc to another point

// Generate Points for Farfield Ogive
For j In {0:points-1}
    x = j * Far_L / (points - 1);  // Distribute points evenly from 0 to Far_L
    Y = (Far_R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / Far_L) - (Sin(2 * Acos(1 - (2 * x) / Far_L)) / 2));
    
    // Define a point with calculated coordinates, adjusting for alignment
    Point(6 + j) = {x - (Far_L - length), Y, 0, 1.0};
EndFor

// Define Flow Lines
Line(5) = {1, 6}; // Start of flow domain
Line(6) = {5, points + 5}; // End of flow domain
Line(7) = {3, midpoint}; // Line from wall to midpoint of farfield

// Create Splines for Farfield Ogive
Spline(8) = {6:midpoint};  // First part of the farfield boundary
Spline(9) = {midpoint:(points+5)};  // Second part of the farfield boundary

// Define Surface Boundaries
Curve Loop(1) = {9, -6, -3, 7}; // Loop for Plane Surface 1
Plane Surface(1) = {1};         // Define Plane Surface 1

Curve Loop(2) = {8, -7, -2, 5}; // Loop for Plane Surface 2
Plane Surface(2) = {2};         // Define Plane Surface 2

// Mesh Refinement Settings
Transfinite Curve {9, 3} = axial * 2.77 Using Progression 1; // Refine mesh along curves 9 and 3
Transfinite Curve {8, 2} = axial Using Progression 1;        // Refine mesh along curves 8 and 2
Transfinite Curve {5, 7, 6} = radial Using Progression gr;   // Refine mesh along curves 5, 7, and 6

// Mesh Surface Control
Transfinite Surface {1};  // Apply transfinite meshing to Surface 1
Transfinite Surface {2};  // Apply transfinite meshing to Surface 2
Recombine Surface {2, 1}; // Recombine surfaces for quadrilateral elements
