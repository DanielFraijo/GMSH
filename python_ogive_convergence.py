import os
import subprocess
import shutil
import numpy as np
from scipy.optimize import root_scalar

# Define file names
geo_filename = "python_ogive.geo"

# Define geometry parameters
base_ogive_radius = 0.5               # Radius of the leading edge [m]
outer_ogive_radius = base_ogive_radius * 5                     # Wedge angle in degrees
base_ogive_length = 1            # Factor to calculate farfield radius
outer_ogive_length = base_ogive_length * 2             # Factor to calculate farfield length
haack_coefficient = 0           # Haack series shape parameter [0-1]
num_ogive_points = 1000                  # Number of points for farfield ogive
boundary_layer_height = 1e-6             # Height of the first boundary layer cell [m]
# Grid convergence parameters
radial_cells_start = 50 
axial_cells_start = 50
delta_cells = 25
max_cells = 75

# Calculate Derived Lengths
line_3_length = outer_ogive_length - base_ogive_length
line_4_length = outer_ogive_radius - base_ogive_radius

# Print the calculated lengths
print(f"Line 3 Length: {line_3_length:.6f}")
print(f"Line 4 Length: {line_4_length:.6f}")

# Growth Rate Calculation Function
def calculate_growth_rate(total_length, initial_cell_size, num_cells):
    def growth_rate_equation(growth_rate):
        if np.isclose(growth_rate, 1.0):
            return initial_cell_size * (num_cells - 1) - total_length
        else:
            return initial_cell_size * (1 - growth_rate**(num_cells - 1)) / (1 - growth_rate) - total_length
    result = root_scalar(growth_rate_equation, bracket=[1.0 + 1e-6, 1.2], method='brentq')
    return result.root


# Create meshes for different grid resolutions
for i in range(0, (max_cells - radial_cells_start) // delta_cells + 1):
    num_radial_cells = radial_cells_start + (i * delta_cells)
    num_axial_cells = axial_cells_start + (i * delta_cells)
    
    # Update SU2 filename for this iteration
    su2_filename = f"ogive_r{num_radial_cells}_a{num_axial_cells}.su2"
    
    # Recalculate growth rates for new cell counts
    GR3 = calculate_growth_rate(line_3_length, boundary_layer_height, num_radial_cells)
    GR4 = calculate_growth_rate(line_4_length, boundary_layer_height, num_radial_cells)
    
    print(f"\nGenerating mesh with {num_radial_cells} radial cells and {num_axial_cells} axial cells")
    print(f"Growth rates - GR3: {GR3:.6f}, GR4: {GR4:.6f}")
    
    # Generate geo_content with updated parameters
    geo_content = f"""// Parameters for the geometry
R = {{base_ogive_radius}};  // Radius of the base
H = {{outer_ogive_radius}};  // Height/base radius of the second ogive
C = {{haack_coefficient}};    // Haack constant (set to 0 for von Kármán ogive)
L = {{base_ogive_length}};    // Length of the first ogive
U = {{outer_ogive_length}}; // Extended length for the second ogive
N = {{num_ogive_points}}; // Number of points for approximation of the curve

radial_cells = {{num_radial_cells}};
axial_cells = {{num_axial_cells}};

gr3 = {{GR3}};
gr4 = {{GR4}};

// Start point for the first ogive
Point(1) = {{L, 0, 0, 1.0}};

// Create points for the first ogive (wall)
For i In {{0:N-1}}
    x = i * L / (N - 1);  // Distribute points evenly from 0 to L
    Y = (R / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / L) - (Sin(2 * Acos(1 - (2 * x) / L)) / 2) + C * Sin(Acos(1 - (2 * x) / L))^3);
    
    // Define a point with calculated coordinates
    Point(i + 2) = {{x, Y, 0, 1.0}};
EndFor

// Connect points for the first ogive with a spline
Spline(1) = {{2:N+1}};  // Corrected to include all points from 2 to N+1

// Create points for the second ogive (farfield)
For j In {{0:N-1}}
    x = j * U / (N - 1);  // Distribute points evenly from 0 to U
    Y = (H / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / U) - (Sin(2 * Acos(1 - (2 * x) / U)) / 2) + C * Sin(Acos(1 - (2 * x) / U))^3);
    
    // Define a point with calculated coordinates, adjusting for the shift in x
    Point(j + N + 2) = {{x - (U - L), Y, 0, 1.0}};
EndFor

// Connect points for the second ogive with a spline
Spline(2) = {{N+2:2*N+1}};  // Creates a spline from points N+2 to 2*N+1

// Define lines to connect splines
Line(3) = {{2, N + 2}}; // Connect start of first ogive to start of second ogive
Line(4) = {{N + 1, 2*N + 1}}; // Connect end of first ogive to end of second ogive

// Define the surface
Curve Loop(1) = {{2, -4, -1, 3}}; // Loop defining the outer boundary of the surface
Plane Surface(1) = {{1}}; // Create the plane surface from the curve loop
Physical Surface("domain") = {{1}};  // Include both surfaces in the domain


// Mesh refinement settings
Transfinite Curve {{2, 1}} = axial_cells Using Progression 1; // Refine mesh along splines
Transfinite Curve {{3}} = radial_cells Using Progression gr3; // Refine mesh along connecting lines
Transfinite Curve {{4}} = radial_cells Using Progression gr4; // Refine mesh along connecting lines
Transfinite Surface {{1}}; // Apply transfinite meshing to the surface
Recombine Surface {{1}};   // Recombine the surface for quadrilateral elements

// Define physical entities for boundary conditions
Physical Curve("farfield", 5) = {{2}};
Physical Curve("symmetry", 6) = {{3}};
Physical Curve("wall", 7) = {{1}};
Physical Curve("outlet", 8) = {{4}};

// Generate the mesh
Mesh 2;

// Set the mesh format to SU2 (optional, format code 42)
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
        
        if os.path.exists(su2_filename):
            print(f"Successfully created SU2 mesh file: {su2_filename}")
        
    except Exception as e:
        print(f"Error generating mesh {su2_filename}: {e}")
