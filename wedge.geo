//=============================================================================
// WEDGE GEOMETRY GENERATOR
// 
// Purpose: Generates a 2D wedge geometry with farfield boundary for CFD analysis
// Author: [Your Name]
// Date: [Current Date]
//
// This script creates:
// 1. A wedge with curved leading edge
// 2. Farfield boundary using Haack series ogive shape
// 3. Structured quadrilateral mesh with boundary layer refinement
//=============================================================================

//-----------------------------------------------------------------------------
// WEDGE GEOMETRY PARAMETERS
//-----------------------------------------------------------------------------
// Basic wedge dimensions
wallRadius = 0.01;                    // Leading edge radius [m]
wedgeAngle = 15 * Pi / 180;          // Wedge angle [rad]
wedgeLength = 0.04;                   // Total wedge length [m]

// Farfield domain parameters
farfieldRadius = 3.85 * wallRadius;   // Farfield boundary radius [m]
farfieldLength = wedgeLength * 1.3;   // Farfield domain length [m]

// Ogive curve discretization
ogivePoints = 1000;                   // Number of points for farfield ogive
ogiveMidpoint = (ogivePoints * 0.2) + 5;  // Splitting point for mesh control
haackConstant = 0.7;                  // Haack series shape parameter [0-1]

//-----------------------------------------------------------------------------
// MESH CONTROL PARAMETERS
//-----------------------------------------------------------------------------
// Mesh density controls
radialCells = 200;
axialCells = 200;                     // Cells along axial direction
firstLayerHeight = 1E-6;             // First cell height for boundary layer [m]

// Growth rates for specific mesh lines
gr5 = 1.01;  // Inlet boundary growth rate
gr6 = 1.01;  // Outlet boundary growth rate
gr7 = 1.01;  // Middle connector growth rate

//-----------------------------------------------------------------------------
// GEOMETRY CONSTRUCTION
//-----------------------------------------------------------------------------
// Define key geometry points
Point(1) = {0, 0, 0};                // Origin/leading edge center
Point(2) = {wallRadius, 0, 0};        // Leading edge radius point
Point(3) = {wallRadius - wallRadius * Sin(wedgeAngle),    // End of curved section
            wallRadius * Cos(wedgeAngle), 0};
Point(4) = {wedgeLength, 0, 0};       // Wedge end point (bottom)
Point(5) = {wedgeLength,              // Wedge end point (top)
            (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) 
            * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle), 0};

// Define Geometrical Shapes
Circle(2) = {1, 2, 3}; // Creates a circular arc
Line(3) = {3, 5};       // Connects the end of the arc to another point

//-----------------------------------------------------------------------------
// FARFIELD BOUNDARY GENERATION
//-----------------------------------------------------------------------------
// Generate Haack series ogive points
For j In {0:ogivePoints-1}
    x = j * farfieldLength / (ogivePoints - 1);
    // Haack series equation for ogive shape
    Y = (farfieldRadius / Sqrt(Pi)) * 
        Sqrt(Acos(1 - (2 * x) / farfieldLength) - 
             (Sin(2 * Acos(1 - (2 * x) / farfieldLength)) / 2) + 
             haackConstant * Sin(Acos(1 - (2 * x) / farfieldLength))^3);
    
    Point(6 + j) = {x - (farfieldLength - wedgeLength), Y, 0, 1.0};
EndFor

// Define Flow Lines
Line(5) = {1, 6}; // Start of flow domain
Line(6) = {5, ogivePoints + 5}; // End of flow domain
Line(7) = {3, ogiveMidpoint}; // Line from wall to midpoint of farfield

// Create Splines for Farfield Ogive
Spline(8) = {6:ogiveMidpoint};  // First part of the farfield boundary
Spline(9) = {ogiveMidpoint:(ogivePoints+5)};  // Second part of the farfield boundary

//-----------------------------------------------------------------------------
// DIAGNOSTIC CALCULATIONS
//-----------------------------------------------------------------------------
// Calculate and output important geometric measurements
// Length for Line 5
length_5 = farfieldLength - wedgeLength;
Printf("Line 5 length: %g", length_5);

line6_x = 0; // x-distance is 0 since both points are at wedgeLength
line6_y_start = (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle);
line6_y_end = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - 2) - (Sin(2 * Acos(1 - 2)) / 2) + haackConstant * Sin(Acos(1 - 2))^3);
line6_length = Abs(line6_y_end - line6_y_start);
Printf("Line 6 length: %g", line6_length);

// length of line 7
// Calculate x and y for midpoint
x_calculation = (ogiveMidpoint - 6) * farfieldLength / (ogivePoints - 1);  
x_midpoint = x_calculation - (farfieldLength - wedgeLength);
y_midpoint = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * (x_calculation) / farfieldLength)) - (Sin(2 * Acos(1 - (2 * (x_calculation) / farfieldLength))) / 2) + haackConstant * Sin(Acos(1 - (2 * (x_calculation) / farfieldLength)))^3);
x_point_3 = wallRadius - wallRadius * Sin(wedgeAngle);
y_point_3 = wallRadius * Cos(wedgeAngle);
// calculations for length
length_7 = Sqrt((Abs(x_midpoint) + x_point_3)^2 + (y_midpoint - y_point_3)^2);
Printf("Line 7 length: %g", length_7);
Printf("Point 3: x = %g, y = %g", x_point_3, y_point_3);
Printf("Midpoint: x = %g, y = %g", x_midpoint, y_midpoint);

//-----------------------------------------------------------------------------
// MESH GENERATION AND CONTROL
//-----------------------------------------------------------------------------
// Define mesh density on curves
Transfinite Curve {9, 3} = axialCells * 2.77 Using Progression 1;
Transfinite Curve {8, 2} = axialCells Using Progression 1;
Transfinite Curve {5} = radialCells Using Progression gr5;
Transfinite Curve {7} = radialCells Using Progression gr7;
Transfinite Curve {6} = radialCells Using Progression gr6;

// Define Surface Boundaries
Curve Loop(1) = {9, -6, -3, 7}; // Loop for Plane Surface 1
Plane Surface(1) = {1};         // Define Plane Surface 1

Curve Loop(2) = {8, -7, -2, 5}; // Loop for Plane Surface 2
Plane Surface(2) = {2};         // Define Plane Surface 2

// Mesh Surface Control
Transfinite Surface {1};  // Apply transfinite meshing to Surface 1
Transfinite Surface {2};  // Apply transfinite meshing to Surface 2
Recombine Surface {2, 1}; // Recombine surfaces for quadrilateral elements

//-----------------------------------------------------------------------------
// BOUNDARY CONDITIONS
//-----------------------------------------------------------------------------
Physical Curve("farfield", 5) = {8, 9};   // Farfield boundary
Physical Curve("symmetry", 6) = {5};       // Symmetry plane
Physical Curve("wall", 7) = {2, 3};        // Wedge surface
Physical Curve("outlet", 8) = {6};         // Outlet boundary
