# cpoly
Testing some algorithms with polygons.<br/>
Sample in OpenGL provided, which depends on GLFW (included and also in https://github.com/glfw/glfw).

# Polygon partitioning

# Convex generation

# Test for point inside a convex polygon
```c++
// Returns 0 if it's not convex. 1 for CW for CCW
  int cpoly_is_convex(void* pts, int npts, int stride);
```
# Test for convex polygon intersection
```c++
// Return 1 if two convex polygons intersects (Separating Axis Theorem)
int cpoly_cv_intersects_SAT(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);
```

# Test for convexity
```c++
// checks if all consecutive edges have same sign of the z component of their cross products
// all negative => convex CW  (return 1)
// all positive => convex CCW (return 2)
// mix sign => not convex (return 0)
int cpoly_is_convex(void* pts, int npts, int stride);
```
