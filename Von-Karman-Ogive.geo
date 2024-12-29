// Parameters for the geometry
R = 0.5;  // Radius of the base
H = R * 3;
C = 0;  // Haack constant (set to 0 for von Kármán ogive)
L = 1; // Length of the ogive
U = L * 1.5;
N = 10000; // Number of points for approximation of the curve

// Start point
Point(1) = {L, 0, 0, 1.0};

// Loop to create points at intervals along x for the first ogive
For i In {0:N-1}
    x = i * L / (N - 1);  // Distribute points evenly from 0 to L
    Y = (R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / L) - (Sin(2 * Acos(1 - (2 * x) / L)) / 2) + C * Sin(Acos(1 - (2 * x) / L))^3);
    
    // Define a point with calculated coordinates
    Point(i + 2) = {x, Y, 0, 1.0};
EndFor

// Create spline for the first ogive
Spline(1) = {2:N+1};  // Corrected to include all points from 2 to N+1

// Loop to create points at intervals along x for the second ogive
For j In {0:N-1}
    x = j * U / (N - 1);  // Distribute points evenly from 0 to U
    Y = (H / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / U) - (Sin(2 * Acos(1 - (2 * x) / U)) / 2) + C * Sin(Acos(1 - (2 * x) / U))^3);
    
    // Define a point with calculated coordinates, adjusting for the shift in x
    Point(j + N + 2) = {x - (U - L), Y, 0, 1.0};
EndFor

// Optionally, connect points to form a line for the second ogive
Spline(2) = {N+2:2*N+1};  // Creates a spline from points N+2 to 2*N+1
Line(3) = {2, N + 2};
Line(4) = {N + 1, 2*N + 1};//+
Curve Loop(1) = {2, -4, -1, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {2, 1} = 100 Using Progression 1;
//+
Transfinite Curve {3, 4} = 150 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Physical Curve("Farfield", 5) = {2};
//+
Physical Curve("Symmetry", 6) = {3};
//+
Physical Curve("Wall", 7) = {1};
//+
Physical Curve("Outlet", 8) = {4};
