//=============================================================================
// WEDGE GEOMETRY GENERATOR
// This script generates a 2D wedge geometry with farfield boundary for CFD analysis
//=============================================================================

//-----------------------------------------------------------------------------
// GEOMETRY PARAMETERS
//-----------------------------------------------------------------------------
wallRadius = 0.01;                    // Radius of the wall curvature (meters)
wedgeAngle = 15 * Pi / 180;          // Wedge angle (converted from degrees to radians)
wedgeLength = 0.04;                  // Total length of the wedge (meters)
farfieldRadius = 3.85 * wallRadius;  // Radius of the farfield boundary
farfieldLength = wedgeLength * 1.5;  // Length of the farfield domain
ogivePoints = 1000;                  // Number of points for farfield ogive curve
ogiveMidpoint = (ogivePoints * 0.2) + 5;  // Midpoint for ogive curve splitting
haackConstant = 0.7;                 // Haack series constant for ogive shape

//-----------------------------------------------------------------------------
// MESH CONTROL PARAMETERS
//-----------------------------------------------------------------------------
axialCells = 200;                    // Number of cells along axial direction
firstLayerHeight = 1E-6;             // Height of first cell layer (for boundary layer)
growthRate = 1.01;                   // Cell growth rate
// Calculate required radial cells based on growth rate and first layer height
radialCells = Ceil(1 + Log(1 - 0.005 * (1 - growthRate)/firstLayerHeight)/Log(growthRate));

// Use fixed growth rates for each line
gr5 = 1.01;  // Growth rate for line 5
gr6 = 1.01;  // Growth rate for line 6
gr7 = 1.01;  // Growth rate for line 7

// Define Points for Geometry
Point(1) = {0, 0, 0};
Point(2) = {wallRadius, 0, 0};
Point(3) = {wallRadius - wallRadius * Sin(wedgeAngle), wallRadius * Cos(wedgeAngle), 0};
Point(4) = {wedgeLength, 0, 0};
Point(5) = {wedgeLength, (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle), 0};

// Define Geometrical Shapes
Circle(2) = {1, 2, 3}; // Creates a circular arc
Line(3) = {3, 5};       // Connects the end of the arc to another point

// Generate Points for Farfield Ogive with Haack Series
For j In {0:ogivePoints-1}
    x = j * farfieldLength / (ogivePoints - 1);  // Distribute points evenly from 0 to farfieldLength
    Y = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / farfieldLength) - (Sin(2 * Acos(1 - (2 * x) / farfieldLength)) / 2) + haackConstant * Sin(Acos(1 - (2 * x) / farfieldLength))^3);
    
    // Define a point with calculated coordinates, adjusting for alignment
    Point(6 + j) = {x - (farfieldLength - wedgeLength), Y, 0, 1.0};
EndFor

// Define Flow Lines
Line(5) = {1, 6}; // Start of flow domain
Line(6) = {5, ogivePoints + 5}; // End of flow domain
Line(7) = {3, ogiveMidpoint}; // Line from wall to midpoint of farfield

// Create Splines for Farfield Ogive
Spline(8) = {6:ogiveMidpoint};  // First part of the farfield boundary
Spline(9) = {ogiveMidpoint:(ogivePoints+5)};  // Second part of the farfield boundary

// Define Surface Boundaries
Curve Loop(1) = {9, -6, -3, 7}; // Loop for Plane Surface 1
Plane Surface(1) = {1};         // Define Plane Surface 1

Curve Loop(2) = {8, -7, -2, 5}; // Loop for Plane Surface 2
Plane Surface(2) = {2};         // Define Plane Surface 2

// Mesh Refinement Settings
Transfinite Curve {9, 3} = axialCells * 2.77 Using Progression 1;
Transfinite Curve {8, 2} = axialCells Using Progression 1;
Transfinite Curve {5} = radialCells Using Progression gr5;
Transfinite Curve {7} = radialCells Using Progression gr7;
Transfinite Curve {6} = radialCells Using Progression gr6;

// Mesh Surface Control
Transfinite Surface {1};  // Apply transfinite meshing to Surface 1
Transfinite Surface {2};  // Apply transfinite meshing to Surface 2
Recombine Surface {2, 1}; // Recombine surfaces for quadrilateral elements

// Mesh Boundaries 
Physical Curve("farfield", 5) = {8, 9};
Physical Curve("symmetry", 6) = {5};
Physical Curve("wall", 7) = {2, 3};
Physical Curve("outlet", 8) = {6};

// Calculate line lengths
// Line 5: From origin to first ogive point
line5_x = -(farfieldLength - wedgeLength); // x-coordinate of first ogive point
line5_y = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - 0) - (Sin(2 * Acos(1 - 0)) / 2) + haackConstant * Sin(Acos(1 - 0))^3);
line5_length = Sqrt(line5_x^2 + line5_y^2);
Printf("Line 5 length: %g", line5_length);

// Line 6: From wedge end to last ogive point
line6_x = 0; // x-distance is 0 since both points are at wedgeLength
line6_y_start = (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle);
line6_y_end = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - 2) - (Sin(2 * Acos(1 - 2)) / 2) + haackConstant * Sin(Acos(1 - 2))^3);
line6_length = Abs(line6_y_end - line6_y_start);
Printf("Line 6 length: %g", line6_length);

// Line 7: From wall curve end to midpoint of ogive
line7_x_start = wallRadius - wallRadius * Sin(wedgeAngle);
line7_y_start = wallRadius * Cos(wedgeAngle);
line7_x_end = (ogiveMidpoint - 6) * farfieldLength / (ogivePoints - 1) - (farfieldLength - wedgeLength);
line7_y_end = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * line7_x_end) / farfieldLength) - (Sin(2 * Acos(1 - (2 * line7_x_end) / farfieldLength)) / 2) + haackConstant * Sin(Acos(1 - (2 * line7_x_end) / farfieldLength))^3);
line7_length = Sqrt((line7_x_end - line7_x_start)^2 + (line7_y_end - line7_y_start)^2);
Printf("Line 7 length: %g", line7_length);


