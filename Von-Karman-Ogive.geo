// Parameters for the geometry
R = 0.5;  // Radius of the base
H = R * 5;  // Height/base radius of the second ogive
C = 0;    // Haack constant (set to 0 for von Kármán ogive)
L = 1;    // Length of the first ogive
U = L * 2; // Extended length for the second ogive
N = 10000; // Number of points for approximation of the curve

// Start point for the first ogive
Point(1) = {L, 0, 0, 1.0};

// Create points for the first ogive (wall)
For i In {0:N-1}
    x = i * L / (N - 1);  // Distribute points evenly from 0 to L
    Y = (R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / L) - (Sin(2 * Acos(1 - (2 * x) / L)) / 2) + C * Sin(Acos(1 - (2 * x) / L))^3);
    
    // Define a point with calculated coordinates
    Point(i + 2) = {x, Y, 0, 1.0};
EndFor

// Connect points for the first ogive with a spline
Spline(1) = {2:N+1};  // Corrected to include all points from 2 to N+1

// Create points for the second ogive (farfield)
For j In {0:N-1}
    x = j * U / (N - 1);  // Distribute points evenly from 0 to U
    Y = (H / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / U) - (Sin(2 * Acos(1 - (2 * x) / U)) / 2) + C * Sin(Acos(1 - (2 * x) / U))^3);
    
    // Define a point with calculated coordinates, adjusting for the shift in x
    Point(j + N + 2) = {x - (U - L), Y, 0, 1.0};
EndFor

// Connect points for the second ogive with a spline
Spline(2) = {N+2:2*N+1};  // Creates a spline from points N+2 to 2*N+1

// Define lines to connect splines
Line(3) = {2, N + 2}; // Connect start of first ogive to start of second ogive
Line(4) = {N + 1, 2*N + 1}; // Connect end of first ogive to end of second ogive

// Define the surface
Curve Loop(1) = {2, -4, -1, 3}; // Loop defining the outer boundary of the surface
Plane Surface(1) = {1}; // Create the plane surface from the curve loop

// Mesh refinement settings
Transfinite Curve {2, 1} = 100 Using Progression 1; // Refine mesh along splines
Transfinite Curve {3, 4} = 150 Using Progression 1; // Refine mesh along connecting lines
Transfinite Surface {1}; // Apply transfinite meshing to the surface
Recombine Surface {1};   // Recombine the surface for quadrilateral elements

// Define physical entities for boundary conditions
Physical Curve("Farfield", 5) = {2};
Physical Curve("Symmetry", 6) = {3};
Physical Curve("Wall", 7) = {1};
Physical Curve("Outlet", 8) = {4};
