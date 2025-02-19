import os
import subprocess
import numpy as np
from scipy.optimize import root_scalar

# Define file names
geo_filename = "python_wedge.geo"
su2_filename = "iteration_2.su2"

# Define geometry parameters
r_wall = 0.01
angle = np.radians(15)
length = 0.04
r_far = 3.85 * r_wall
l_far = length * 1.3
haack = 0.7
n_ogive = 1000
h_first = 1e-6
n_radial = 200

# Calculate line lengths
l5 = l_far - length
print(f"Length 5: {l5}")

y_start = (length - (r_wall - r_wall * np.sin(angle))) * np.tan(angle) + r_wall * np.cos(angle)
y_end = (r_far / np.sqrt(np.pi)) * np.sqrt(np.arccos(1 - 2) - (np.sin(2 * np.arccos(1 - 2)) / 2) + haack * np.sin(np.arccos(1 - 2))**3)
l6 = abs(y_end - y_start)
print(f"Length 6: {l6}")

mid_point = int(n_ogive * 0.2) + 5
x_calc = (mid_point - 6) * l_far / (n_ogive - 1)
x_mid = x_calc - (l_far - length)
y_mid = (r_far / np.sqrt(np.pi)) * np.sqrt(np.arccos(1 - (2 * x_calc / l_far)) - (np.sin(2 * np.arccos(1 - (2 * x_calc / l_far))) / 2) + haack * np.sin(np.arccos(1 - (2 * x_calc / l_far)))**3)
x3 = r_wall - r_wall * np.sin(angle)
y3 = r_wall * np.cos(angle)
l7 = np.sqrt((abs(x_mid) + x3)**2 + (y_mid - y3)**2)
print(f"Length 7: {l7}")

# Calculate growth rates
def calculate_growth_rate(total_length, ds, N):
    def f(gr):
        if np.isclose(gr, 1.0):
            return ds * (N - 1) - total_length
        else:
            return ds * (1 - gr**(N-1)) / (1 - gr) - total_length
    result = root_scalar(f, bracket=[1.0 + 1e-6, 1.2], method='brentq')
    return result.root

GR5 = calculate_growth_rate(l5, h_first, n_radial)
GR6 = calculate_growth_rate(l6, h_first, n_radial)
GR7 = calculate_growth_rate(l7, h_first, n_radial)

# Print the growth rates
print(f"Growth rate 5: {GR5:.6f}")
print(f"Growth rate 6: {GR6:.6f}")
print(f"Growth rate 7: {GR7:.6f}")

# Create the .geo file with Wedge geometry
geo_content = f"""//=============================================================================
// WEDGE GEOMETRY GENERATOR
// 
// Purpose: Generates a 2D wedge geometry with farfield boundary for CFD analysis
// Author: Daniel Fraijo
// Date: Febuuary 14, 2025
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
gr5 = {GR5};  // Inlet boundary growth rate
gr6 = {GR6};  // Outlet boundary growth rate
gr7 = {GR7};  // Middle connector growth rate

//-----------------------------------------------------------------------------
// GEOMETRY CONSTRUCTION
//-----------------------------------------------------------------------------
// Define key geometry points
Point(1) = {{0, 0, 0}};                // Origin/leading edge center
Point(2) = {{wallRadius, 0, 0}};        // Leading edge radius point
Point(3) = {{wallRadius - wallRadius * Sin(wedgeAngle), wallRadius * Cos(wedgeAngle), 0}};
Point(4) = {{wedgeLength, 0, 0}};       // Wedge end point (bottom)
Point(5) = {{wedgeLength, (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle), 0}};

// Define Geometrical Shapes
Circle(2) = {{1, 2, 3}}; // Creates a circular arc
Line(3) = {{3, 5}};       // Connects the end of the arc to another point

//-----------------------------------------------------------------------------
// FARFIELD BOUNDARY GENERATION
//-----------------------------------------------------------------------------
// Generate Haack series ogive points
For j In {{0:ogivePoints-1}}
    x = j * farfieldLength / (ogivePoints - 1);
    // Haack series equation for ogive shape
    Y = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / farfieldLength) - (Sin(2 * Acos(1 - (2 * x) / farfieldLength)) / 2) + haackConstant * Sin(Acos(1 - (2 * x) / farfieldLength))^3);
    Point(6 + j) = {{x - (farfieldLength - wedgeLength), Y, 0, 1.0}};
EndFor

// Define Flow Lines
Line(5) = {{1, 6}}; // Start of flow domain
Line(6) = {{5, ogivePoints + 5}}; // End of flow domain
Line(7) = {{3, ogiveMidpoint}}; // Line from wall to midpoint of farfield

// Create Splines for Farfield Ogive
Spline(8) = {{6:ogiveMidpoint}};  // First part of the farfield boundary
Spline(9) = {{ogiveMidpoint:(ogivePoints+5)}};  // Second part of the farfield boundary

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
Transfinite Curve {{9, 3}} = axialCells * 2.77 Using Progression 1;
Transfinite Curve {{8, 2}} = axialCells Using Progression 1;
Transfinite Curve {{5}} = radialCells Using Progression gr5;
Transfinite Curve {{7}} = radialCells Using Progression gr7;
Transfinite Curve {{6}} = radialCells Using Progression gr6;

// Define Surface Boundaries
Curve Loop(1) = {{9, -6, -3, 7}}; // Loop for Plane Surface 1
Plane Surface(1) = {{1}};         // Define Plane Surface 1
Curve Loop(2) = {{8, -7, -2, 5}}; // Loop for Plane Surface 2
Plane Surface(2) = {{2}};         // Define Plane Surface 2

// Mesh Surface Control
Transfinite Surface {{1}};  // Apply transfinite meshing to Surface 1
Transfinite Surface {{2}};  // Apply transfinite meshing to Surface 2
Recombine Surface {{2, 1}}; // Recombine surfaces for quadrilateral elements

//-----------------------------------------------------------------------------
// BOUNDARY CONDITIONS
//-----------------------------------------------------------------------------
Physical Curve("farfield", 5) = {{8, 9}};   // Farfield boundary
Physical Curve("symmetry", 6) = {{5}};       // Symmetry plane
Physical Curve("wall", 7) = {{2, 3}};        // Wedge surface
Physical Curve("outlet", 8) = {{6}};         // Outlet boundary
"""

# Write .geo file and run Gmsh
try:
    # Write the .geo file
    with open(geo_filename, "w") as f:
        f.write(geo_content)
    print(f"Successfully created {geo_filename}")

    # Check if gmsh is available
    gmsh_check = subprocess.run(["which", "gmsh"], capture_output=True, text=True)
    if gmsh_check.returncode != 0:
        raise FileNotFoundError("Gmsh not found in system PATH")

    # Run Gmsh with detailed output
    result = subprocess.run(
        ["gmsh", "-2", geo_filename, "-format", "su2", "-o", su2_filename],
        capture_output=True,
        text=True,
        check=True
    )
    
    # Check if output file was created
    if os.path.exists(su2_filename):
        print(f"Successfully created SU2 mesh file: {su2_filename}")
    else:
        print("Warning: SU2 file not created despite successful Gmsh execution")
        
    # Print Gmsh output if there's any
    if result.stdout:
        print("Gmsh output:")
        print(result.stdout)

except FileNotFoundError as e:
    print(f"Error: {e}")
    print("Please ensure Gmsh is installed and available in your system PATH")
except subprocess.CalledProcessError as e:
    print(f"Error running Gmsh (exit code {e.returncode}):")
    print(e.stderr)
except Exception as e:
    print(f"Unexpected error: {e}")
