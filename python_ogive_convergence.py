import os
import subprocess
import numpy as np
from scipy.optimize import root_scalar

# Define file names
geo_filename = "python_ogive.geo"

# Define geometry parameters
base_ogive_radius = 0.5               # Radius of the base [m]
outer_ogive_radius = base_ogive_radius * 5  # Farfield radius factor
base_ogive_length = 1                 # Length of the first ogive [m]
outer_ogive_length = base_ogive_length * 2  # Farfield length factor
haack_coefficient = 0                 # Haack constant (0 for von Kármán ogive)
num_ogive_points = 1000               # Number of points for ogive curves
boundary_layer_height = 1e-6          # Height of the first boundary layer cell [m]

# Grid convergence parameters
radial_cells_start = 50 
axial_cells_start = 50
delta_cells = 25
max_cells = 75

# Calculate derived lengths
line_3_length = outer_ogive_length - base_ogive_length  # Symmetry boundary length
line_4_length = outer_ogive_radius - base_ogive_radius  # Outlet boundary length

print(f"Line 3 Length: {line_3_length:.6f}")
print(f"Line 4 Length: {line_4_length:.6f}")

# Growth Rate Calculation Function
def calculate_growth_rate(total_length, initial_cell_size, num_cells):
    def growth_rate_equation(growth_rate):
        if np.isclose(growth_rate, 1.0):
            return initial_cell_size * (num_cells - 1) - total_length
        else:
            return initial_cell_size * (1 - growth_rate**(num_cells - 1)) / (1 - growth_rate) - total_length

    # Initial bracket
    r_min, r_max = 1.0 + 1e-6, 1.2
    f_min = growth_rate_equation(r_min)
    f_max = growth_rate_equation(r_max)

    # Adjust bracket if necessary
    if f_min * f_max > 0:
        if f_min > 0:
            r_min_new = 1.0 + 1e-10
            f_min_new = growth_rate_equation(r_min_new)
            if f_min_new < 0:
                r_min = r_min_new
            else:
                raise ValueError(f"No root found: f({r_min_new}) = {f_min_new}, f({r_max}) = {f_max}")
        elif f_max < 0:
            r_max_new = 1.5
            f_max_new = growth_rate_equation(r_max_new)
            if f_max_new > 0:
                r_max = r_max_new
            else:
                raise ValueError(f"No root found: f({r_min}) = {f_min}, f({r_max_new}) = {f_max_new}")

    result = root_scalar(growth_rate_equation, bracket=[r_min, r_max], method='brentq')
    return result.root

# Create meshes for different grid resolutions
for i in range(0, (max_cells - radial_cells_start) // delta_cells + 1):
    num_radial_cells = radial_cells_start + (i * delta_cells)
    num_axial_cells = axial_cells_start + (i * delta_cells)
    
    # Update SU2 filename
    su2_filename = f"ogive_r{num_radial_cells}_a{num_axial_cells}.su2"
    
    # Calculate growth rates
    GR3 = calculate_growth_rate(line_3_length, boundary_layer_height, num_radial_cells)
    GR4 = calculate_growth_rate(line_4_length, boundary_layer_height, num_radial_cells)
    
    print(f"\nGenerating mesh with {num_radial_cells} radial cells and {num_axial_cells} axial cells")
    print(f"Growth rates - GR3: {GR3:.6f}, GR4: {GR4:.6f}")
    
    # Generate .geo content with wedge-style comments
    geo_content = f"""//=============================================================================
// OGIVE GEOMETRY GENERATOR
// 
// Purpose: Generates a 2D ogive geometry with farfield boundary for CFD analysis
// Author: [Your Name]
// Date: [Current Date]
//
// This script creates:
// 1. An ogive-shaped body with a specified base radius and length
// 2. A farfield boundary using a larger ogive shape
// 3. Structured quadrilateral mesh with boundary layer refinement
//=============================================================================

//-----------------------------------------------------------------------------
// OGIVE GEOMETRY PARAMETERS
//-----------------------------------------------------------------------------
R = {base_ogive_radius};  // Base radius of the ogive [m]
H = {outer_ogive_radius};  // Farfield radius [m]
C = {haack_coefficient};   // Haack series shape parameter [0-1]
L = {base_ogive_length};   // Length of the ogive body [m]
U = {outer_ogive_length};  // Farfield length [m]
N = {num_ogive_points};    // Number of points for ogive curves

//-----------------------------------------------------------------------------
// MESH CONTROL PARAMETERS
//-----------------------------------------------------------------------------
radial_cells = {num_radial_cells};  // Number of cells in radial direction
axial_cells = {num_axial_cells};    // Number of cells in axial direction
first_layer_height = {boundary_layer_height};  // First cell height [m]
gr3 = {GR3};  // Growth rate for symmetry boundary (Line 3)
gr4 = {GR4};  // Growth rate for outlet boundary (Line 4)

// Pre-calculated lengths for reference:
// Line 3 (symmetry): {line_3_length:.6f} m
// Line 4 (outlet): {line_4_length:.6f} m

//-----------------------------------------------------------------------------
// GEOMETRY CONSTRUCTION
//-----------------------------------------------------------------------------
// Define points for the ogive body (wall)
For i In {{0:N-1}}
    x = i * L / (N - 1);
    Y = (R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / L) - (Sin(2 * Acos(1 - (2 * x) / L)) / 2) + C * Sin(Acos(1 - (2 * x) / L))^3);
    Point(i + 2) = {{x, Y, 0, 1.0}};
EndFor

// Create spline for the ogive body (wall)
Spline(1) = {{2:N+1}};

// Define points for the farfield ogive
For j In {{0:N-1}}
    x = j * U / (N - 1);
    Y = (H / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / U) - (Sin(2 * Acos(1 - (2 * x) / U)) / 2) + C * Sin(Acos(1 - (2 * x) / U))^3);
    Point(j + N + 2) = {{x - (U - L), Y, 0, 1.0}};
EndFor

// Create spline for the farfield boundary
Spline(2) = {{N+2:2*N+1}};

// Connect the splines with lines
Line(3) = {{2, N + 2}};  // Symmetry boundary: connects wall start to farfield start
Line(4) = {{N + 1, 2*N + 1}};  // Outlet boundary: connects wall end to farfield end

// Define the surface using the curve loop
Curve Loop(1) = {{2, -4, -1, 3}};
Plane Surface(1) = {{1}};
Physical Surface("domain") = {{1}};

//-----------------------------------------------------------------------------
// MESH GENERATION AND CONTROL
//-----------------------------------------------------------------------------
// Set transfinite curves for mesh density control
Transfinite Curve {{2, 1}} = axial_cells Using Progression 1;  // Along axial direction (farfield and wall)
Transfinite Curve {{3}} = radial_cells Using Progression gr3;  // Radial direction (symmetry boundary)
Transfinite Curve {{4}} = radial_cells Using Progression gr4;  // Radial direction (outlet boundary)

// Apply transfinite meshing to the surface
Transfinite Surface {{1}};

// Recombine for quadrilateral elements
Recombine Surface {{1}};

//-----------------------------------------------------------------------------
// BOUNDARY CONDITIONS
//-----------------------------------------------------------------------------
Physical Curve("farfield", 5) = {{2}};   // Farfield boundary
Physical Curve("symmetry", 6) = {{3}};   // Symmetry plane
Physical Curve("wall", 7) = {{1}};       // Ogive surface
Physical Curve("outlet", 8) = {{4}};     // Outlet boundary

// Generate the mesh
Mesh 2;

// Set the mesh format to SU2 (format code 42)
Mesh.Format = 42;
"""
    
    try:
        # Write the .geo file
        with open(geo_filename, "w") as f:
            f.write(geo_content)
        print(f"Successfully created {geo_filename}")

        # Run Gmsh to generate the SU2 mesh
        result = subprocess.run(
            ["gmsh", "-2", geo_filename, "-format", "su2", "-o", su2_filename, "-v", "4"],
            capture_output=True,
            text=True,
            check=True
        )
        print("Gmsh stdout:", result.stdout)
        print("Gmsh stderr:", result.stderr)
        
        if os.path.exists(su2_filename):
            print(f"Successfully created SU2 mesh file: {su2_filename}")
        else:
            print(f"Mesh file {su2_filename} was not created.")
        
    except Exception as e:
        print(f"Error generating mesh {su2_filename}: {e}")
