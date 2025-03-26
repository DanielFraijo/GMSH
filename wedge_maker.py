import os
import subprocess
import numpy as np
from scipy.optimize import root_scalar

# Define file names
geo_filename = "python_wedge.geo"

# Define geometry parameters
leading_edge_radius = 0.01               # Radius of the leading edge [m]
wedge_angle_deg = 15                     # Wedge angle in degrees
wedge_length = 0.04                      # Total length of the wedge [m]
farfield_radius_factor = 4.5            # Factor to calculate farfield radius
farfield_length_factor = 1.3             # Factor to calculate farfield length
haack_series_coefficient = 0.7           # Haack series shape parameter [0-1]
num_ogive_points = 1000                  # Number of points for farfield ogive

# Boundary layer height sweep parameters
min_boundary_layer_height = 1e-6         # Minimum height of first boundary layer cell [m]
max_boundary_layer_height = 1e-3         # Maximum height of first boundary layer cell [m]
num_sweep_points = 3                     # Number of points in the sweep

# Grid convergence parameters
radial_cells_start = 100 
axial_cells_start = 100
delta_cells = 100
max_cells = 300

# Generate logarithmic sweep of boundary layer heights
boundary_layer_heights = np.logspace(np.log10(min_boundary_layer_height), 
                                    np.log10(max_boundary_layer_height), 
                                    num_sweep_points)

# Calculate Derived Lengths
farfield_length = wedge_length * farfield_length_factor
length_difference = farfield_length - wedge_length
print(f"Length Difference: {length_difference}")

# Calculate Y-Coordinates for Line 6
y_start = (wedge_length - (leading_edge_radius - leading_edge_radius * np.sin(np.radians(wedge_angle_deg)))) \
          * np.tan(np.radians(wedge_angle_deg)) + leading_edge_radius * np.cos(np.radians(wedge_angle_deg))
y_end_expression = 1 - 2
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

# Updated Growth Rate Calculation Function
def calculate_growth_rate(total_length, initial_cell_size, num_cells):
    def growth_rate_equation(growth_rate):
        if np.isclose(growth_rate, 1.0):
            return initial_cell_size * (num_cells - 1) - total_length
        else:
            return initial_cell_size * (1 - growth_rate**(num_cells - 1)) / (1 - growth_rate) - total_length

    # Check if uniform sizing is feasible
    uniform_size = initial_cell_size * (num_cells - 1)
    if uniform_size > total_length:
        print(f"Warning: Initial cell size {initial_cell_size} with {num_cells} cells exceeds total length {total_length}. Adjusting...")
        return 1.0  # Return uniform growth rate if overshoot occurs

    # Initial bracket
    r_min = 1.0 + 1e-6  # Slightly above 1
    r_max = 1.2         # Reasonable upper bound
    f_min = growth_rate_equation(r_min)
    f_max = growth_rate_equation(r_max)

    # Adjust bracket dynamically
    if f_min * f_max > 0:  # Same sign
        if f_min > 0:  # Both positive, need smaller r
            r_min_new = 1.0 + 1e-10
            f_min_new = growth_rate_equation(r_min_new)
            if f_min_new < 0:
                r_min = r_min_new
            else:
                print(f"Warning: No root found with initial bracket. Returning approximate growth rate.")
                return 1.0  # Fallback to uniform sizing
        elif f_max < 0:  # Both negative, need larger r
            r_max = 2.0
            f_max = growth_rate_equation(r_max)
            if f_max < 0:
                r_max = 10.0  # Extreme case
                f_max = growth_rate_equation(r_max)
                if f_max < 0:
                    print(f"Warning: No root found with extended bracket. Returning approximate growth rate.")
                    return 1.1  # Arbitrary reasonable growth rate

    try:
        result = root_scalar(growth_rate_equation, bracket=[r_min, r_max], method='brentq')
        return result.root
    except ValueError as e:
        print(f"Root finding failed: {e}. Using fallback growth rate.")
        return 1.05  # Fallback growth rate

# Create meshes for different grid resolutions and boundary layer heights
for bl_height in boundary_layer_heights:
    for i in range(0, (max_cells - radial_cells_start) // delta_cells + 1):
        num_radial_cells = radial_cells_start + (i * delta_cells)
        num_axial_cells = axial_cells_start + (i * delta_cells)
        
        # Format boundary layer height as scientific notation without decimal (e.g., 1e-06)
        bl_str = f"{bl_height:.0e}".replace(".", "")  # Remove decimal point
        # Update SU2 filename with boundary layer height without decimal
        su2_filename = f"wedge_r{num_radial_cells}_a{num_axial_cells}_{bl_str}.su2"
        
        # Recalculate growth rates
        GR5 = calculate_growth_rate(length_difference, bl_height, num_radial_cells)
        GR6 = calculate_growth_rate(length_line6, bl_height, num_radial_cells)
        GR7 = calculate_growth_rate(length_line7, bl_height, num_radial_cells)
        
        print(f"\nGenerating mesh with {num_radial_cells} radial cells and {num_axial_cells} axial cells")
        print(f"Boundary layer height: {bl_height:.1e}")
        print(f"Growth rates - GR5: {GR5:.6f}, GR6: {GR6:.6f}, GR7: {GR7:.6f}")
        
        # Generate geo_content with updated parameters
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
        wallRadius = {leading_edge_radius};
        wedgeAngle = {wedge_angle_deg} * Pi / 180;
        wedgeLength = {wedge_length};

        Printf("wall radius: %g", wallRadius);
        Printf("wedge angle: %g", wedgeAngle);
        Printf("wedge length: %g", wedgeLength);

        farfieldRadius = {farfield_radius_factor} * wallRadius;
        farfieldLength = wedgeLength * {farfield_length_factor};

        Printf("farfield radius: %g", farfieldRadius);
        Printf("farfield length: %g", farfieldLength);

        ogivePoints = {num_ogive_points};
        ogiveMidpoint = (ogivePoints * 0.2) + 5;
        haackConstant = {haack_series_coefficient};

        Printf("ogive points: %g", ogivePoints);
        Printf("ogive midpoint: %g", ogiveMidpoint);
        Printf("haack constant: %g", haackConstant);

        //-----------------------------------------------------------------------------
        // MESH CONTROL PARAMETERS
        //-----------------------------------------------------------------------------
        radialCells = {num_radial_cells};
        axialCells = {num_axial_cells};
        firstLayerHeight = {bl_height};

        Printf("radial cells: %g", radialCells);
        Printf("axial cells: %g", axialCells);
        Printf("first layer height: %g", firstLayerHeight);

        gr5 = {GR5};
        gr6 = {GR6};
        gr7 = {GR7};

        Printf("growth rate 5: %g", gr5);
        Printf("growth rate 6: %g", gr6);
        Printf("growth rate 7: %g", gr7);

        //-----------------------------------------------------------------------------
        // GEOMETRY CONSTRUCTION
        //-----------------------------------------------------------------------------
        Point(1) = {{0, 0, 0}};
        Point(2) = {{wallRadius, 0, 0}};
        Point(3) = {{wallRadius - wallRadius * Sin(wedgeAngle), wallRadius * Cos(wedgeAngle), 0}};
        Point(4) = {{wedgeLength, 0, 0}};
        Point(5) = {{wedgeLength, (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle), 0}};

        Circle(2) = {{1, 2, 3}};
        Line(3) = {{3, 5}};

        //-----------------------------------------------------------------------------
        // FARFIELD BOUNDARY GENERATION
        //-----------------------------------------------------------------------------
        For j In {{0:ogivePoints-1}}
            x = j * farfieldLength / (ogivePoints - 1);
            Y = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * x) / farfieldLength) - (Sin(2 * Acos(1 - (2 * x) / farfieldLength)) / 2) + haackConstant * Sin(Acos(1 - (2 * x) / farfieldLength))^3);
            Point(6 + j) = {{x - (farfieldLength - wedgeLength), Y, 0, 1.0}};
        EndFor

        Line(5) = {{1, 6}};
        Line(6) = {{5, ogivePoints + 5}};
        Line(7) = {{3, ogiveMidpoint}};

        Spline(8) = {{6:ogiveMidpoint}};
        Spline(9) = {{ogiveMidpoint:(ogivePoints+5)}};

        //-----------------------------------------------------------------------------
        // DIAGNOSTIC CALCULATIONS
        //-----------------------------------------------------------------------------
        length_5 = farfieldLength - wedgeLength;
        Printf("Line 5 length: %g", length_5);

        line6_x = 0;
        line6_y_start = (wedgeLength - (wallRadius - wallRadius * Sin(wedgeAngle))) * Tan(wedgeAngle) + wallRadius * Cos(wedgeAngle);
        line6_y_end = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - 2) - (Sin(2 * Acos(1 - 2)) / 2) + haackConstant * Sin(Acos(1 - 2))^3);
        line6_length = Abs(line6_y_end - line6_y_start);
        Printf("Line 6 length: %g", line6_length);

        x_calculation = (ogiveMidpoint - 6) * farfieldLength / (ogivePoints - 1);  
        x_midpoint = x_calculation - (farfieldLength - wedgeLength);
        y_midpoint = (farfieldRadius / Sqrt(Pi)) * Sqrt(Acos(1 - (2 * (x_calculation) / farfieldLength)) - (Sin(2 * Acos(1 - (2 * (x_calculation) / farfieldLength))) / 2) + haackConstant * Sin(Acos(1 - (2 * (x_calculation) / farfieldLength)))^3);
        x_point_3 = wallRadius - wallRadius * Sin(wedgeAngle);
        y_point_3 = wallRadius * Cos(wedgeAngle);
        length_7 = Sqrt((Abs(x_midpoint) + x_point_3)^2 + (y_midpoint - y_point_3)^2);
        Printf("Line 7 length: %g", length_7);
        Printf("Point 3: x = %g, y = %g", x_point_3, y_point_3);
        Printf("Midpoint: x = %g, y = %g", x_midpoint, y_midpoint);

        //-----------------------------------------------------------------------------
        // MESH GENERATION AND CONTROL
        //-----------------------------------------------------------------------------
        Transfinite Curve {{9, 3}} = axialCells * 2.77 Using Progression 1;
        Transfinite Curve {{8, 2}} = axialCells Using Progression 1;
        Transfinite Curve {{5}} = radialCells Using Progression gr5;
        Transfinite Curve {{7}} = radialCells Using Progression gr7;
        Transfinite Curve {{6}} = radialCells Using Progression gr6;

        Curve Loop(1) = {{9, -6, -3, 7}};
        Plane Surface(1) = {{1}};
        Curve Loop(2) = {{8, -7, -2, 5}};
        Plane Surface(2) = {{2}};
        Physical Surface("domain") = {{1, 2}};

        Transfinite Surface {{1}};
        Transfinite Surface {{2}};
        Recombine Surface {{2, 1}};

        //-----------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //-----------------------------------------------------------------------------
        Physical Curve("farfield", 5) = {{8, 9}};
        Physical Curve("symmetry", 6) = {{5}};
        Physical Curve("wall", 7) = {{2, 3}};
        Physical Curve("outlet", 8) = {{6}};

        Mesh 2;
        Mesh.Format = 42;
        """
        
        try:
            with open(geo_filename, "w") as f:
                f.write(geo_content)
            print(f"Successfully created {geo_filename}")

            result = subprocess.run(
                ["gmsh", "-2", geo_filename, "-format", "su2", "-o", su2_filename, "-v", "4"],
                capture_output=True,
                text=True,
                check=True
            )
            if os.path.exists(su2_filename):
                print(f"Successfully created SU2 mesh file: {su2_filename}")
        
        except subprocess.CalledProcessError as e:
            print(f"Error generating mesh {su2_filename}: {e}")
            print(f"Gmsh output: {e.stdout}")
            print(f"Gmsh error: {e.stderr}")
        except Exception as e:
            print(f"Unexpected error generating mesh {su2_filename}: {e}")
