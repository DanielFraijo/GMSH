import numpy as np
from scipy.optimize import root_scalar

# Parameters
T = 0.005       # Total thickness
ds = 1e-6       # Initial spacing
N = 120         # Number of points

# Define the function for the geometric series sum
def f(gr):
    if np.isclose(gr, 1.0):
        # Handle gr=1 case (uniform spacing)
        return ds * (N - 1) - T
    else:
        return ds * (1 - gr**(N-1)) / (1 - gr) - T

# Solve using Brent's method (bracketed root-finding)
# Bracket gr between 1.0 and 1.1 (adjust if needed)
result = root_scalar(f, bracket=[1.0 + 1e-6, 1.1], method='brentq')
gr = result.root

print(f"The growth rate gr is: {gr:.10f}")  # 10 decimal places
