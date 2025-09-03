// Parameters for the geometry
R = 0.05;  // Radius of the base
H = R * 15;  // Height/base radius of the second ogive
C = 0;    // Haack constant (set to 0 for von Kármán ogive)
L = 1;    // Length of the first ogive
U = L * 1.2; // Extended length for the second ogive
N = 10000; // Number of points for approximation of the curve

// Mesh parameters
axial_nodes = 6000;
radial_nodes = 250;
first_cell_height = 1e-6;

// Growth rates calculated for first cell height
growth_left = 1.0363552776644864;  // For inlet/symmetry (Curve 3)
growth_right = 1.0428266347260264;  // For outlet (Curve 4)

// Create points for the first ogive (wall)
For i In {0:N-1}
    x = i * L / (N - 1);
    Y = (R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / L) - (Sin(2 * Acos(1 - (2 * x) / L)) / 2) + C * Sin(Acos(1 - (2 * x) / L))^3);
    Point(i + 2) = {x, Y, 0, 1.0};
EndFor

// Connect points for the first ogive with a spline
Spline(1) = {2:N+1};

// Create points for the second ogive (farfield)
For j In {0:N-1}
    x = j * U / (N - 1);
    Y = (H / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / U) - (Sin(2 * Acos(1 - (2 * x) / U)) / 2) + C * Sin(Acos(1 - (2 * x) / U))^3);
    Point(j + N + 2) = {x - (U - L), Y, 0, 1.0};
EndFor

// Connect points for the second ogive with a spline
Spline(2) = {N+2:2*N+1};

// Define lines to connect splines
Line(3) = {2, N + 2}; // From nose to inlet farfield
Line(4) = {N + 1, 2*N + 1}; // From wall end to outlet farfield

// Define the curve loop in counter-clockwise order
Curve Loop(1) = {-3, 1, 4, -2};

// Create the plane surface
Plane Surface(1) = {1};

// Mesh refinement settings
Transfinite Curve {1, 2} = axial_nodes Using Progression 1; // Uniform along axial directions
Transfinite Curve {3} = radial_nodes Using Progression growth_left; // Graded along inlet/symmetry
Transfinite Curve {4} = radial_nodes Using Progression growth_right; // Graded along outlet
Transfinite Surface {1} = {N+2, 2, N+1, 2*N+1} Left; // Explicit corners in order

// Recombine for quadrilateral elements
Recombine Surface {1};

// Additional settings for better mesh quality
Mesh.Algorithm = 8;  // Frontal-Delaunay for quads (corrected from Algorithm2D)
Mesh.RecombinationAlgorithm = 1;  // Blossom for quad quality
Mesh.Smoothing = 100;  // Laplacian smoothing

// Generate the mesh
Mesh 2;

// Apply elliptic smoothing to the transfinite surface (fsolve-like elliptic solver)
Smoother Surface {1} = 500;  // Number of elliptic smoothing iterations

// Further optimize the mesh
Mesh.Optimize = 1;

// Define physical entities for boundary conditions
Physical Curve("Farfield", 5) = {2};
Physical Curve("Symmetry", 6) = {3};
Physical Curve("Wall", 7) = {1};
Physical Curve("Outlet", 8) = {4};

Mesh.Format = 42;
