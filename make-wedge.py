import os
import subprocess

# Define file names
geo_filename = "mesh.geo"
su2_filename = "mesh.su2"

# Create the .geo file with Wedge geometry
geo_content = """//=============================================================================
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

// Function to calculate growth rate using Bisection method
Function GrowthRate
    // Set bounds for growth rate
    r_min = 1.001;  // Minimum allowable growth rate
    r_max = 2.0;    // Maximum allowable growth rate
    r = 1.5;        // Initial guess
    
    n = radialCells;     // Use user-specified number of points
    
    // Target length function
    Function GetLength
        sum = firstLayerHeight * (1 - r^(n-1))/(1-r);
        Return sum - wedgeLength;
    EndFunction
    
    // Check if solution exists in bounds
    If (Call GetLength~{r = r_min} * Call GetLength~{r = r_max} > 0)
        Printf("Warning: No solution exists within bounds!");
        Return 1.001;
    EndIf
    
    // Bisection iteration
    For iter In {1:max_iterations}
        r = (r_min + r_max) / 2;
        length_error = Call GetLength;
        
        // Check convergence
        If (Abs(length_error) < tolerance)
            Break;
        EndIf
        
        // Update bounds
        If (length_error > 0)
            r_max = r;
        Else
            r_min = r;
        EndIf
        
        // Warning if not converged
        If (iter == max_iterations)
            Printf("Warning: Growth rate calculation reached maximum iterations!");
        EndIf
    EndFor
    
    Return r;
EndFunction

// Calculate growth rates for each line using fixed number of points
// Line 5 length: distance from origin to farfield start
L5 = Sqrt(6^2 + Y^2);
L = L5;
gr5 = Call GrowthRate;

// Line 6 length: distance from wall end to farfield end
L6 = Sqrt((ogivePoints+5-5)^2 + Y^2);
L = L6;
gr6 = Call GrowthRate;

// Line 7 length: distance from wall to midpoint of farfield
L7 = Sqrt((ogiveMidpoint-3)^2 + Y^2);
L = L7;
gr7 = Call GrowthRate;

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
"""

# Write .geo file
with open(geo_filename, "w") as f:
    f.write(geo_content)

print(f"✅ Created {geo_filename}")

# Run Gmsh to generate the SU2 file
try:
    subprocess.run(["gmsh", "-2", geo_filename, "-format", "su2", "-o", su2_filename], check=True)
    print(f"✅ SU2 mesh file created: {su2_filename}")
except subprocess.CalledProcessError as e:
    print(f"❌ Error running Gmsh: {e}")

# Add growth rate analysis
def calculate_growth_rate(L1, total_length, N):
    """Calculate the growth rate for geometric grid spacing"""
    def objective(r):
        return L1 * (1 - r**(N-1))/(1-r) - total_length
    
    # Initial guess for growth rate
    r_initial = 1.1
    
    # Simple iterative solver (similar to bisection method in .geo file)
    r_min, r_max = 1.001, 2.0
    r = r_initial
    max_iter = 100
    tolerance = 1e-6
    
    for _ in range(max_iter):
        error = objective(r)
        if abs(error) < tolerance:
            break
        if error > 0:
            r_max = r
        else:
            r_min = r
        r = (r_min + r_max) / 2
    
    return r

def generate_points(L1, total_length, N):
    """Generate points with geometric spacing"""
    r = calculate_growth_rate(L1, total_length, N)
    points = np.zeros(N)
    for i in range(1, N):
        points[i] = points[i-1] + L1 * r**(i-1)
    return points, r

# Analyze each line's growth rate
lines = {
    'Line 5 (Start of flow domain)': {'length': L5, 'cells': radialCells},
    'Line 6 (End of flow domain)': {'length': L6, 'cells': radialCells},
    'Line 7 (Wall to farfield)': {'length': L7, 'cells': radialCells}
}

plt.figure(figsize=(15, 10))

for idx, (name, params) in enumerate(lines.items(), 1):
    points, growth_rate = generate_points(firstLayerHeight, params['length'], params['cells'])
    
    plt.subplot(3, 1, idx)
    plt.plot(points, np.zeros_like(points), '-b', label='Domain')
    plt.plot(points, np.zeros_like(points), 'or', markersize=2, label='Grid points')
    
    plt.title(f'{name} (Growth rate = {growth_rate:.3f})')
    plt.xlabel('Position (m)')
    plt.grid(True)
    plt.legend()
    plt.ylim(-0.1, 0.1)
    
    print(f"\nAnalysis for {name}:")
    print(f"Calculated growth rate: {growth_rate:.3f}")
    print(f"First spacing: {points[1] - points[0]:.2e} m")
    print(f"Last spacing: {points[-1] - points[-2]:.2e} m")
    print(f"Total length covered: {points[-1]:.3f} m")

plt.tight_layout()
plt.show()
