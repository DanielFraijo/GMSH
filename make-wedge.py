import os
import subprocess
import shutil
import numpy as np
from scipy.optimize import root_scalar

# Define file names
geo_filename = "python_wedge.geo"
su2_filename = "iteration_2.su2"

# Define geometry parameters
leading_edge_radius = 0.01               # Radius of the leading edge [m]
wedge_angle_deg = 15                     # Wedge angle in degrees
wedge_length = 0.04                      # Total length of the wedge [m]
farfield_radius_factor = 3.85            # Factor to calculate farfield radius
farfield_length_factor = 1.3             # Factor to calculate farfield length
haack_series_coefficient = 0.7           # Haack series shape parameter [0-1]
num_ogive_points = 1000                  # Number of points for farfield ogive
boundary_layer_height = 1e-6             # Height of the first boundary layer cell [m]
num_radial_cells = 200                    # Number of radial cells
num_axial_cells = 200                    # 

# Calculate Derived Lengths
farfield_length = wedge_length * farfield_length_factor
length_difference = farfield_length - wedge_length
print(f"Length Difference: {length_difference}")

# Calculate Y-Coordinates for Line 6
y_start = (wedge_length - (leading_edge_radius - leading_edge_radius * np.sin(np.radians(wedge_angle_deg)))) \
          * np.tan(np.radians(wedge_angle_deg)) + leading_edge_radius * np.cos(np.radians(wedge_angle_deg))
y_end_expression = 1 - 2  # This seems to be a constant; verify if intended
y_end = (farfield_radius_factor * leading_edge_radius / np.sqrt(np.pi)) \
        * np.sqrt(np.arccos(y_end_expression) - (np.sin(2 * np.arccos(y_end_expression)) / 2) \
        + haack_series_coefficient * np.sin(np.arccos(y_end_expression))**3)
length_line6 = abs(y_end - y_start)
print(f"Length Line 6: {length_line6}")

# Calculate Midpoint Coordinates for Line 7
mid_ogive_index = int(num_ogive_points * 0.2) + 5
x_calculation = (mid_ogive_index - 6) * farfield_length / (num_ogive_points - 1)
x_midpoint = x_calculation - (farfield_length - wedge_length)
y_midpoint = (farfield_radius_factor * leading_edge_radius / np.sqrt(np.pi)) \
             * np.sqrt(np.arccos(1 - (2 * x_calculation / farfield_length)) \
             - (np.sin(2 * np.arccos(1 - (2 * x_calculation / farfield_length))) / 2) \
             + haack_series_coefficient * np.sin(np.arccos(1 - (2 * x_calculation / farfield_length)))**3)
x_point_wall = leading_edge_radius - leading_edge_radius * np.sin(np.radians(wedge_angle_deg))
y_point_wall = leading_edge_radius * np.cos(np.radians(wedge_angle_deg))
length_line7 = np.sqrt((abs(x_midpoint) + x_point_wall)**2 + (y_midpoint - y_point_wall)**2)
print(f"Length Line 7: {length_line7}")

# Growth Rate Calculation Function
def calculate_growth_rate(total_length, initial_cell_size, num_cells):
    def growth_rate_equation(growth_rate):
        if np.isclose(growth_rate, 1.0):
            return initial_cell_size * (num_cells - 1) - total_length
        else:
            return initial_cell_size * (1 - growth_rate**(num_cells - 1)) / (1 - growth_rate) - total_length
    result = root_scalar(growth_rate_equation, bracket=[1.0 + 1e-6, 1.2], method='brentq')
    return result.root

GR5 = calculate_growth_rate(length_difference, boundary_layer_height, num_radial_cells)
GR6 = calculate_growth_rate(length_line6, boundary_layer_height, num_radial_cells)
GR7 = calculate_growth_rate(length_line7, boundary_layer_height, num_radial_cells)

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
wallRadius = {leading_edge_radius};                    // Leading edge radius [m]
wedgeAngle = {wedge_angle_deg} * Pi / 180;          // Wedge angle [rad]
wedgeLength = {wedge_length};                   // Total wedge length [m]

Printf("wall radius: %g", wallRadius);
Printf("wedge angle: %g", wedgeAngle);
Printf("wedge length: %g", wedgeLength);

// Farfield domain parameters
farfieldRadius = {farfield_radius_factor} * wallRadius;   // Farfield boundary radius [m]
farfieldLength = wedgeLength * {farfield_length_factor};   // Farfield domain length [m]

Printf("farfield radius: %g", farfieldRadius);
Printf("farfield length: %g", farfieldLength);

// Ogive curve discretization
ogivePoints = {num_ogive_points};                   // Number of points for farfield ogive
ogiveMidpoint = (ogivePoints * 0.2) + 5;  // Splitting point for mesh control
haackConstant = {haack_series_coefficient};                  // Haack series shape parameter [0-1]

Printf("ogive points: %g", ogivePoints);
Printf("ogive midpoint: %g", ogiveMidpoint);
Printf("haack constant: %g", haackConstant);

//-----------------------------------------------------------------------------
// MESH CONTROL PARAMETERS
//-----------------------------------------------------------------------------
// Mesh density controls
radialCells = {num_radial_cells};
axialCells = {num_axial_cells};                     // Cells along axial direction
firstLayerHeight = {boundary_layer_height};             // First cell height for boundary layer [m]

Printf("radial cells: %g", radialCells);
Printf("axial cells: %g", axialCells);
Printf("first layer height: %g", firstLayerHeight);

// Growth rates for specific mesh lines
gr5 = {GR5};  // Inlet boundary growth rate
gr6 = {GR6};  // Outlet boundary growth rate
gr7 = {GR7};  // Middle connector growth rate

Printf("growth rate 5: %g", gr5);
Printf("growth rate 6: %g", gr6);
Printf("growth rate 7: %g", gr7);

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

    # Check if Gmsh is available
    gmsh_path = shutil.which("gmsh")
    if gmsh_path is None:
        raise FileNotFoundError("Gmsh not found in system PATH")

    # Run Gmsh to generate the SU2 mesh
    result = subprocess.run(
        ["gmsh", "-2", geo_filename, "-format", "su2", "-o", su2_filename, "-v", "4"],
        capture_output=True,
        text=True,
        check=True
    )
    
    # Verify output file
    if os.path.exists(su2_filename):
        print(f"Successfully created SU2 mesh file: {su2_filename}")
    else:
        print("Warning: SU2 file not created despite successful Gmsh execution")
    
    # Print Gmsh output for debugging
    if result.stdout:
        print("Gmsh Output:")
        print(result.stdout)
    if result.stderr:
        print("Gmsh Errors:")
        print(result.stderr)

except FileNotFoundError as e:
    print(f"Error: {e}")
    print("Please ensure Gmsh is installed and available in your system PATH")
except subprocess.CalledProcessError as e:
    print(f"Error running Gmsh (exit code {e.returncode}):")
    print("Gmsh Output:")
    print(e.stdout)
    print("Gmsh Errors:")
    print(e.stderr)
except Exception as e:
    print(f"Unexpected error: {e}")
