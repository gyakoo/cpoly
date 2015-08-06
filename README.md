# cpoly
Testing some algorithms with polygons.<br/>

# Some operations

```c
    // Returns 0 if it's not convex. 1 for CW for CCW
  int cpoly_is_convex(void* pts, int npts, int stride);

  // Returns 1 if the point is inside a convex polygon
  int cpoly_cv_point_inside(void* pts, int npts, int stride, float x, float y);

  // Return 1 if two convex polygons intersects (Separating Axis Theorem)
  int cpoly_cv_intersects_SAT(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // Return 1 if two segment a and b intersects  (optionally returns the intersection point)
  // 's' is the intersection fraction 0..1 from x0,y0 to x1,y1 when there's an intersection
  int cpoly_seg_isec(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float* s);

  // Return 1 if point x,y is inside the triangle given by x0,y0 .. x1,y1 .. x2,y2
  int cpoly_point_in_triangle(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2);
  
  // calculates the bounding polygons for a set of circle points (x,y,radius)...using marching squares algorithm
  // sqside is the side of the square. no interpolation.
  // returns number of polygons created. See NOTES for how to get the results.
  int cpoly_marchingsq_nointerp(void* pts, int npts, int stride, float sqside);

  // triangulate by Ear Clipping method
  int cpoly_triangulate_EC(void* pts, int npts, int stride, void* reserved);

  // Hertel-Mehlhorn partition algorithm. Produces 
  int cpoly_partition_HM(void* pts, int npts, int stride);

  // returns 1 if segment intersects polygon. returns optionally in ix,iy the closest intersection point
  // edge is the edge index in the polygon assuming they're consecutive (edge = start_vertex)
  int cpoly_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // returns 1 if segment intersects polygon and optionally the intersection point in ix/iy and the edge. First intersection found.
  int cpoly_cv_seg_isec_poly_first(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // Union of convex polygons. Returns no. of vertices to be accessed with cpoly_pool_get_vertex
  // assumes: convex polygons, intersecting, no one inside another, no holes.
  // somehow based on Sutherland–Hodgman algorithm
  int cpoly_cv_clip_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);
  
  // Difference of convex polygons, returns no. of parts generated
  // You got the parts offsets starting from part 1 in cpoly_pool_i indices
  // All the vertices generated in cpoly_pool_v (See NOTES how to get results.)  
  // somehow based on Sutherland–Hodgman algorithm
  int cpoly_cv_clip_diff(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // rotate around a pivot, NULL to rotate around center.
  void cpoly_transform_rotate(void* pts, int npts, int stride, float angle, float* xpivot, float* ypivot);

  // scale polygon points along a pivot or using center if NULL.
  void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot);

  // translate center of polygon to x,y or the pivot if not null
  void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y, float* xpivot, float* ypivot);

  // geometric center computation (center of mass or centroid)
  void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy);

  // computes convex hull of a polygon. Returns no of indices to vertices in original polygon, use cpoly_pool_get_index.
  // Gift wrapping algorithm
  int cpoly_convex_hull(void* pts, int npts, int stride);

  // computes axis aligned bounding box
  void cpoly_aabb(void* pts, int npts, int stride, float* xmin, float* ymin, float* xmax, float* ymax);
``` 
