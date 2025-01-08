// Gmsh project updated for cylinder with splines
SetFactory("OpenCASCADE");

// Geometry Parameters (all units in meters)
radius = 0.001;           // Radius of the core geometry
flow_thick = 3 * radius; // Thickness of the flow domain around the geometry
rear_length = radius * 2; 
front_length = radius * 0.3;
ogive_length = front_length + radius; // Total length of the ogive
ogive_radius = 0.0025;
points = 100;
C = 0.2; // Haack constant for the ogive shape

// Mesh Parameters
axial = 200; // Number of axial cells
ds = 1E-6;   // Minimum characteristic length
gr = 1.01;   // Growth rate
radial = Ceil(1 + Log(1 - 0.005 * (1 - gr)/ds)/Log(gr)); // Number of radial cells

// Define Points for Geometry
Point(1) = {0, 0, 0, 1.0};            // Center of the circle
Point(2) = {radius, radius, 0, 1.0}; // Top of the cylinder
Point(3) = {radius, 0, 0, 1.0};      // Right of the circle's center
Point(4) = {2 * radius, 0, 0, 1.0};  // End of the circle
Point(5) = {-front_length, 0, 0, 1.0}; // Front of the domain
Point(6) = {2 * radius + rear_length, 0, 0, 1.0}; // End of the rear
Point(7) = {radius * 2 + rear_length, ogive_radius, 0, 1.0}; // Top of ogive at rear
Point(8) = {radius + rear_length, ogive_radius, 0, 1.0};     // Midpoint of ogive at rear
Point(9) = {radius * 2 + rear_length, ogive_radius - radius, 0, 1.0}; // Lower point at rear

// Generate Points for Farfield Ogive with Haack Series
For j In {0:points-1}
    x = j * ogive_length / (points - 1);  // Distribute points evenly from 0 to ogive_length
    Y = (ogive_radius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / ogive_length) - (Sin(2 * Acos(1 - (2 * x) / ogive_length)) / 2) + C * Sin(Acos(1 - (2 * x) / ogive_length))^3);

    // Define a point with calculated coordinates, aligning with the lengths
    Point(11 + j) = {x - front_length, Y, 0, 1.0};
EndFor

// Define Geometrical Shapes
Circle(1) = {1, 3, 2}; // Circular arc from Point 1 to 2 via 3
Circle(2) = {2, 3, 4}; // Circular arc from Point 2 to 4 via 3
Line(3) = {2, points + 10}; // Vertical line from cylinder to midpoint of ogive
Line(4) = {1, 5};  // Line from center to front of domain
Line(5) = {4, 6};  // Line from end of circle to rear_length
Spline(6) = {5, 12:(points + 10)}; // First part of the ogive boundary
BSpline(7) = {points + 10, 8, 7, 9, 6}; // BSpline for the rear ogive shape

// Mesh Refinement Settings
Transfinite Curve {6, 1} = axial Using Progression 1;
Transfinite Curve {7, 2} = axial Using Progression 1;
Transfinite Curve {4, 3, 5} = radial Using Progression gr;

// Define Surface Boundaries
Curve Loop(1) = {6, -3, -1, 4}; // Loop for the first surface
Plane Surface(1) = {1};         // Define first Surface

Curve Loop(2) = {7, -5, -2, 3}; // Loop for the second surface
Plane Surface(2) = {2};         // Define second Surface

// Mesh Surface Control
Transfinite Surface {1};
Transfinite Surface {2};
Recombine Surface {1, 2}; // Recombine surfaces for quadrilateral elements

Physical Curve("Farfield", 5) = {6};
Physical Curve("Symmetry", 6) = {4, 5};
Physical Curve("Wall", 7) = {1, 2};
Physical Curve("Outlet", 8) = {7};
