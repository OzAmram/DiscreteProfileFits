import ROOT
import ctypes  # Import the ctypes module
from array import array

# Create a sample TGraph
n = 5
x_vals_orig = array('d', [1.0, 2.0, 3.0, 4.0, 5.0])
y_vals_orig = array('d', [2.0, 4.0, 6.0, 8.0, 10.0])
graph = ROOT.TGraph(n, x_vals_orig, y_vals_orig)

# 1. Use ctypes.c_double for variables that will be passed by reference
x_val = ctypes.c_double(0.0)
y_val = ctypes.c_double(0.0)

# Iterate through the points and get their coordinates
print("Points in the TGraph:")
for i in range(n):
    # 2. Pass the ctypes objects directly to GetPoint()
    graph.GetPoint(i, x_val, y_val)
    
    # 3. Access the actual value using the .value attribute of the ctypes object
    print(f"Point {i}: x = {x_val.value:.2f}, y = {y_val.value:.2f}")

# The GetPointX/Y methods are still the easiest way to get individual values:
point_index = 2
x_single = graph.GetPointX(point_index)
y_single = graph.GetPointY(point_index)
print(f"\nCoordinates of point {point_index} using GetPointX/Y: x = {x_single:.2f}, y = {y_single:.2f}")

