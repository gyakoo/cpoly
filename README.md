# cpoly
Testing some algorithms with polygons.<br/>
OpenGL viewer provided. Depends on GLFW which is included but also in https://github.com/glfw/glfw

# Polygon partitioning

# Convex generation

# Test for convexity
```c++
// checks if all consecutive edges have same sign of the z component of their cross products
// all negative => convex CW  (return 1)
// all positive => convex CCW (return 2)
// mix sign => not convex (return 0)
int cpoly_is_convex(void* pts, int npts, int stride);
```
