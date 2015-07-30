/*
Disclaimer:
THIS IS WORK IN PROGRESS, SOME FUNCTIONS AREN'T FINISHED AND/OR HEAVILY TESTED FOR BEHAVIOR AND PERFORMANCE.
 */

#ifndef CPOLY_H
#define CPOLY_H

#include <math.h>
#ifdef _DEBUG
#include <intrin.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

  // Returns 0 if it's not convex. 1 for CW for CCW
  int cpoly_is_convex(void* pts, int npts, int stride);

  // Returns 1 if the point is inside a convex polygon
  int cpoly_cv_point_inside(void* pts, int npts, int stride, float x, float y);

  // Return 1 if two convex polygons intersects (Separating Axis Theorem)
  int cpoly_cv_intersects_SAT(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // Return 1 if two segment a and b intersects  (optionally returns the intersection point)
  // 's' is the intersection fraction 0..1 from x0,y0 to x1,y1 when there's an intersection
  int cpoly_seg_isec(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float* s);

  // returns 1 if segment intersects polygon. returns optionally in ix,iy the closest intersection point
  // edge is the edge index in the polygon assuming they're consecutive (edge = start_vertex)
  int cpoly_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // returns 1 if segment intersects polygon and optionally the intersection point in ix/iy and the edge. First intersection found.
  int cpoly_cv_seg_isec_poly_first(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // Union of convex polygons. Returns no. of vertices to be accessed with cpoly_pool_get_vertex
  // assumes: convex polygons, intersecting, no one inside another, no holes.
  int cpoly_cv_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);
  
  // Difference of convex polygons, returns no. of parts generated
  // You got the parts offsets starting from part 1 in cpoly_pool_i indices
  // All the vertices generated in cpoly_pool_v (See NOTES how to get results.)  
  int cpoly_cv_diff(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // rotate around a pivot, NULL to rotate around center.
  void cpoly_transform_rotate(void* pts, int npts, int stride, float angle, float* xpivot, float* ypivot);

  // scale polygon points along a pivot or using center if NULL.
  void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot);

  // translate center of polygon to x,y or the pivot if not null
  void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y, float* xpivot, float* ypivot);

  // geometric center computation (center of mass or centroid)
  void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy);

  // computes convex hull of a polygon. Returns no of indices to vertices in original polygon, use cpoly_pool_get_index
  int cpoly_convex_hull(void* pts, int npts, int stride);

  // computes axis aligned bounding box
  void cpoly_aabb(void* pts, int npts, int stride, float* xmin, float* ymin, float* xmax, float* ymax);

  // calculates the bounding polygons for a set of circle points (x,y,radius)...using marching squares algorithm with no interpolation
  // sqside is the side of the square
  // returns number of polygons created. See NOTES for how to get the results.
  int cpoly_marchingsq_nointerp(void* pts, int npts, int stride, float sqside);

  // computes the convex partition of an arbitrary partition (Hertel-Mehlhorn). See NOTES how to get results.
  int cpoly_convex_partition_HM(void* pts, int npts, int stride);

  // convext partition
  int cpoly_convex_partition(void* pts, int npts, int stride);


  /* NOTES
    - For functions returning +1 polygon vertices, it uses poly_pool_v as vertex buffer and cpoly_pool_i as array with polygon counts    
    Example:

    cpoly_cv_diff(poly0, N0, sizeof(float)*2, poly1, N1, sizeof(float)*2); // performs difference of polygons poly0-poly1
    k=0;
    for ( i = 0; i < cpoly_pool_icount; ++i ) // for all resulting polygons
    {
      glBegin(GL_LINE_LOOP);
      for ( j=k; j<cpoly_pool_get_index(i); ++j )
      {
          cpoly_pool_get_vertex(j,&x,&y); // Polygon i, vertex j
          glVertex2f( x, y);
      }
      glEnd();
      k=j;
    }
  */


#ifdef __cplusplus
};
#endif

#endif // CPOLY_H

#ifdef CPOLY_IMPLEMENTATION

/*
TODO: 
  [DONE] http://jamie-wong.com/2014/08/19/metaballs-and-marching-squares/
  http://www.dma.fi.upm.es/docencia/trabajosfindecarrera/programas/geometriacomputacional/PiezasConvex/algoritmo_i.html
  http://bl.ocks.org/mbostock
*/

#ifdef _MSC_VER
# pragma warning (disable: 4996) // Switch off security warnings
# pragma warning (disable: 4100) // Switch off unreferenced formal parameter warnings
# pragma warning (disable: 4204) // Switch off nonstandard extension used : non-constant aggregate initializer
#endif

#ifndef CPOLY_MAXPOOLSIZE_BYTES 
#define CPOLY_MAXPOOLSIZE_BYTES (1<<20)
#else
#if CPOLY_MAXPOOLSIZE_BYTES <= 32
#undef CPOLY_MAXPOOLSIZE_BYTES
#define CPOLY_MAXPOOLSIZE_BYTES (16<<10)
#endif
#endif

// Used for bitsets and Marching Squares Grid
// MSG: Olny 4 bits are used to represent a grid cell value (0-15)
typedef struct cpolyBitPool
{
  unsigned char* values;
  int stride;
}cpolyBitPool;

// NEGATIVE is to the eye (CW) POSITIVE is towards the screen (CCW)
#define cpoly_zcross(x0,y0, x1, y1, x2, y2) ((x1-x0)*(y2-y1) - (y1-y0)*(x2-x1))

// returns xy from a stream of floats (pointer, stride, index)
#define cpoly_getx(p,s,n) (*(    (float*)( (char*)(p) + ((s)*(n)) ) ))
#define cpoly_gety(p,s,n) (*(    (float*)( (char*)(p) + ((s)*(n)+sizeof(float)) ) ))
#define cpoly_getz(p,s,n) (*(    (float*)( (char*)(p) + ((s)*(n)+sizeof(float)*2) ) ))
#define cpoly_getxy(p,s,n,x,y) { x=cpoly_getx(p,s,n); y=cpoly_gety(p,s,n); }
#define cpoly_getxyz(p,s,n,x,y,z) { x=cpoly_getx(p,s,n); y=cpoly_gety(p,s,n); z=cpoly_getz(p,s,n);}
#define cpoly_setx(p,s,n,f) cpoly_getx(p,s,n)=f
#define cpoly_sety(p,s,n,f) cpoly_gety(p,s,n)=f
#define cpoly_setxy(p,s,n,x,y) { cpoly_setx(p,s,n,x); cpoly_sety(p,s,n,y); }

// non thread safe shared pools (vertex and index)
unsigned char cpoly_pool_v[CPOLY_MAXPOOLSIZE_BYTES];
unsigned char cpoly_pool_i[CPOLY_MAXPOOLSIZE_BYTES];
const int CPOLY_POOL_VMAX = CPOLY_MAXPOOLSIZE_BYTES/(sizeof(float)*2); // max no of vertices (x,y) in the pool
const int CPOLY_POOL_IMAX= CPOLY_MAXPOOLSIZE_BYTES/sizeof(int); // max no of indices
int cpoly_pool_vcount=0;
int cpoly_pool_icount=0;

//////////////////////////////////////////////////////////////////////////
// forward decls
//////////////////////////////////////////////////////////////////////////
int   cpoly_max(int a, int b);
int   cpoly_min(int a, int b);
int   cpoly_clampi(int v, int m, int M);
void  cpoly_pool_add_vertex(float x, float y);
void  cpoly_pool_add_index(int ndx);
void  cpoly_pool_get_vertex(int n, float* x, float* y);
int   cpoly_pool_get_index(int n);
void  cpoly_bitset_create(cpolyBitPool* bs, int nbits);
void  cpoly_bitset_destroy(cpolyBitPool* bs);
void  cpoly_bitset_set(cpolyBitPool* bs, int n, int value);
int   cpoly_bitset_get(cpolyBitPool* bs, int n);
void  cpoly_msg_create(cpolyBitPool* rmg, int width, int height);
void  cpoly_msg_destroy(cpolyBitPool* rmg);
int   cpoly_msg_get(cpolyBitPool* rmg, int i, int j);
void  cpoly_msg_set(cpolyBitPool* rmg, int i, int j, int val);
float cpoly_msg_func(void* pts, int npts, int stride, float x, float y);
int   cpoly_msg_cellvalue(float c0, float c1, float c2, float c3); //corners: bottomleft, bottomright, topright, topleft

//////////////////////////////////////////////////////////////////////////
// IMPLEMENTATION 
//////////////////////////////////////////////////////////////////////////

// checks if all consecutive edges have same sign of the z component of their cross products
// all negative => convex CW  (return 1)
// all positive => convex CCW (return 2)
// mix sign => not convex (return 0)
int cpoly_is_convex(void* pts, int npts, int stride)
{
  int i;
  float s;
  float curz=.0f, lastz=.0f;
  float x0,y0,x1,y1,x2,y2;
  
  if ( npts <= 2 ) return 0; // not a polygon
  
  for (i=0;i<npts-1;++i)
  {
    cpoly_getxy(pts,stride,i,x0,y0);
    cpoly_getxy(pts,stride,i+1,x1,y1);
    cpoly_getxy(pts,stride,(i+2)%npts,x2,y2);

    // z comp of cross product of both edges
    curz = cpoly_zcross(x0,y0,x1,y1,x2,y2);

    // changed sign?
    s=curz*lastz;
    if ( s < .0f ) return 0;
    lastz = curz;
  }
  return (lastz < .0f) ? 1 : 2;
}

// check if point inside a CONVEX polygon
int cpoly_cv_point_inside(void* pts, int npts, int stride, float x, float y)
{
  float _x0,_y0,x0,y0,x1,y1,x2,y2;
  float zcross, curz;
  int i;

  if ( npts <= 2 ) return 0; // not a polygon

  // check order (cw or ccw)
  cpoly_getxy(pts,stride,0,x0,y0); _x0=x0; _y0=y0;
  cpoly_getxy(pts,stride,1,x1,y1);
  cpoly_getxy(pts,stride,2,x2,y2);
  zcross = cpoly_zcross(x0,y0,x1,y1,x2,y2);

  // we have already the first three points
  curz = cpoly_zcross(x0,y0,x1,y1,x,y); if ( (curz*zcross) < 0.0f ) return 0;
  curz = cpoly_zcross(x1,y1,x2,y2,x,y); if ( (curz*zcross) < 0.0f ) return 0;

  // iterate all remaining points
  x0=x2; y0=y2;
  for ( i=3; i<npts;++i)
  {
    cpoly_getxy(pts,stride,i,x1,y1);
    curz = cpoly_zcross(x0,y0,x1,y1,x,y); if ( (curz*zcross) < 0.0f ) return 0;
    x0=x1; y0=y1;
  }

  // closes loop with the 1st
  curz = cpoly_zcross(x0,y0,_x0,_y0,x,y); if ( (curz*zcross) < 0.0f ) return 0;
  return 1;
}

// Separating Axis Theorem
// The MTV is the minimum magnitude vector used to push the shapes out of the collision.
//#define CPOLY_ISECT_SAT_FINDMTV
int cpoly_cv_intersects_SAT(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1/*=-1*/)
{
  // for all edges in polygon 0 + polygon 1
    // normal to edge
    // project all points of polygon 0 and 1 onto normal, computing their min and max values
    // check if min/max for both in the normal, overlaps. 
    // if no overlap => return 0
  
  const int npts[2]={npts0,npts1};
  const void* pts[2]={pts0,pts1};
  const int stride[2]={stride0,stride1<=0?stride0:stride1};
  float minis[2];
  float maxis[2];
  float nx,ny;
  float x0,y0,x1,y1;
  float xp,yp;
  float d;
  int P,p;
  int e,v;
#ifdef CPOLY_ISECT_SAT_FINDMTV
  float minovlap=FLT_MAX,ovlap;
  float mtvx,mtvy;
#endif

  // edges in both polygons
  for ( P=0; P < 2; ++ P)
  {
    // for all edges of P
    for (e=0; e < npts[P]; ++e)
    {
      cpoly_getxy(pts[P], stride[P], e, x0, y0);
      cpoly_getxy(pts[P], stride[P], (e+1)%npts[P], x1, y1);
      // edge normal
      nx = y1-y0; ny = -(x1-x0);
      
      // projecting all points from P0 and P1 onto normal and getting min and maxs
      minis[0]=minis[1]=FLT_MAX; maxis[0]=maxis[1]=-FLT_MAX;
      for (p=0;p<2;++p)
      {
        for (v=0;v<npts[p];++v)
        {
          cpoly_getxy(pts[p], stride[p], v, xp, yp);
          // project along normal and getting min/max
          d=xp*nx + yp*ny;
          if ( d < minis[p] ) minis[p]=d; if ( d > maxis[p] ) maxis[p]=d;
        }
      }

      // if no overlapping of min/max for both polygons, we can return 0 (we found a separating axis)
      if ( maxis[0]<minis[1] || maxis[1]<minis[0] )
        return 0;
#ifdef CPOLY_ISECT_SAT_FINDMTV
      else
      {
        // find minimum traslation vector (mtv) 
        ovlap = (maxis[0]>minis[1])?(maxis[0]-minis[1]):(maxis[1]-minis[0]); // overlapping len
        if ( ovlap < minovlap )
        {
          minovlap = ovlap;
          mtvx=nx; mtvy=ny;
        }
      }
#endif
    }
  }
  // we went through all edges and no separating axis found, we have an intersection
  return 1;
}

int cpoly_seg_isec(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float* it)
{
  const float a=x2-x0; const float b=x1-x0;
  const float c=x3-x2; const float d=y1-y0;
  const float e=y3-y2; const float f=y2-y0;
  const float inv = 1.0f/(c*d-e*b);
  const float t = (f*c-a*e)*inv;
  const float s = (f*b-a*d)*inv;

  // if are indeterminate, they're collinear
  // if are infinite, they're parallel
  // in both cases no inters and will return 0
  // segments, so they'd be in the range
  if ( s>=0.0f && s<=1.0f && t>=0.0f && t<=1.0f )
  {
    if ( it ) *it= t;
    return 1;
  }  
  return 0;
}

int __internal_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge, char breakfirst)
{
  int v;
  float x2,y2,x3,y3;
  float s, mins=FLT_MAX;
  int tmp; 
  if (!edge) edge=&tmp;
  *edge=-1;
  for (v=0;v<npts;++v)
  {
    cpoly_getxy(pts,stride,v,x2,y2);
    cpoly_getxy(pts,stride,(v+1)%npts,x3,y3);
    if ( cpoly_seg_isec(x0,y0,x1,y1,x2,y2,x3,y3,&s) && s<mins )
    {
      mins=s;      
      *edge=v;
      if ( breakfirst )
        break;
    }
  }
  if ( *edge== -1 )
    return 0;

  if ( ix ) *ix = x0 + (x1-x0)*mins;
  if ( iy ) *iy = y0 + (y1-y0)*mins;
  return 1;
}

float g_rayOffset=0.0f;
int cpoly_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge)
{
  float dx=0,dy=0;
  if ( g_rayOffset!=0.0f )
  {
    dx=x1-x0; dy=y1-y0;
    x0 += dx*g_rayOffset; y0 += dy*g_rayOffset;
    x1 -= dx*g_rayOffset; y1 -= dy*g_rayOffset;
  }
  return __internal_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts,npts,stride,ix,iy,edge,0);
}

int cpoly_cv_seg_isec_poly_first(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge)
{
  float dx=0,dy=0;
  if ( g_rayOffset!=0.0f )
  {
    dx=x1-x0; dy=y1-y0;
    x0 += dx*g_rayOffset; y0 += dy*g_rayOffset;
    x1 -= dx*g_rayOffset; y1 -= dy*g_rayOffset;
  }

  return __internal_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts,npts,stride,ix,iy,edge,1);
}

void cpoly_pool_add_vertex(float x, float y)
{
  float* ptr = (float*)( cpoly_pool_v+cpoly_pool_vcount*sizeof(float)*2);

#ifdef _DEBUG
  if ( cpoly_pool_vcount >= CPOLY_POOL_VMAX )
    __debugbreak();
#endif
  *(ptr++) = x; *ptr = y;
  ++cpoly_pool_vcount;
}

void cpoly_pool_add_index(int ndx)
{
  int* ptr = (int*)( cpoly_pool_i+cpoly_pool_icount*sizeof(int));
  *ptr = ndx;
  ++cpoly_pool_icount;
}

void cpoly_pool_get_vertex(int n, float* x, float* y)
{
  float* ptr = (float*)( cpoly_pool_v+n*sizeof(float)*2);
  *x=*(ptr++); *y=*ptr;
}

int cpoly_pool_get_index(int n)
{
  return *(int*)( cpoly_pool_i+n*sizeof(int));
}

void cpoly_bitset_create(cpolyBitPool* bs, int nbits)
{
  bs->stride = (int)ceilf( (float)nbits/8.0f );
  bs->values = (unsigned char*)calloc(bs->stride, 1);
}

void cpoly_bitset_destroy(cpolyBitPool* bs)
{
  bs->stride=0;
  free(bs->values);
  bs->values=0;
}

void cpoly_bitset_set(cpolyBitPool* bs, int n, int value)
{
  bs->values[n>>3] ^= value<<(n%8);
}

int cpoly_bitset_get(cpolyBitPool* bs, int n)
{
  const int nmod8=n%8;
  return (bs->values[n>>3] & (1<<nmod8)) >> nmod8;
}

// assumes: convex polygons, no one inside another, no holes, clockwise????
int cpoly_cv_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1)
{
  const int npts[2]={npts0,npts1};
  void* pts[2]={pts0,pts1};
  const int stride[2]={stride0,stride1<=0?stride0:stride1};
  int P,v,e,oP;
  float minx=FLT_MAX;
  int minv, minp;
  float x0,y0,x1,y1,ix,iy;
  int counts[2]={0,0};
  int inside;
  cpolyBitPool bitsets[2];

  // some more checks here
  if (stride1<=0) stride1=stride0;
  cpoly_bitset_create(bitsets, npts0);
  cpoly_bitset_create(bitsets+1,npts1);

  // pick start vertex out of two polygons (can be the one with min X if we use a straight Vertical sweep line)
  for (P=0;P<2;++P)
  {
    oP = (P+1)%2;
    for (v=0;v<npts[P];++v )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);
      
      // computes how many points of P are outside of oP && mark flag 1 for those points that are in the overlapping region
      inside = cpoly_cv_point_inside(pts[oP],npts[oP],stride[oP],x0,y0);
      if ( !inside ) 
      { 
        ++counts[P]; 
        if ( x0<minx ) { minp=P; minv=v; minx=x0; }
      } 
      cpoly_bitset_set(bitsets+P,v, inside);
    }
  }

  if ( counts[0]==npts[0] && counts[1]==npts[1] ) return 0; // no overlap
  // start with that min vertex
  v = minv; P = minp; oP=(P+1)%2;
  cpoly_pool_vcount = 0;

  // while still outside points to process for both
  while ( counts[0]>0 || counts[1]>0 )
  {
    if ( !cpoly_bitset_get(bitsets+P,v) )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);
      cpoly_pool_add_vertex(x0,y0); --counts[P]; // adds the first one in the edge
      cpoly_bitset_set(bitsets+P,v,1); // marks this vertex as outside, so can't be further processed

      cpoly_getxy(pts[P],stride[P],(v+1)%npts[P],x1,y1);

      // for this edge, intersects with the any of the other polygon's edge?
      if ( cpoly_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts[oP],npts[oP],stride[oP],&ix,&iy,&e) )
      {
        cpoly_pool_add_vertex(ix,iy);        
        P=oP; oP=(oP+1)%2; // the other poly
        v=e; //continues with next poly vertex, the end of the edge
      }      
    } 
    else if ( !counts[P] ) // no more vertices to process in current polygon, jump to the other
    {
      P=oP; oP=(oP+1)%2;
    }
    v = (v+1)%npts[P];
  }

  cpoly_bitset_destroy(bitsets+1);
  cpoly_bitset_destroy(bitsets);
  return cpoly_pool_vcount;
}

int cpoly_max(int a, int b){ return a>b?a:b; }

int cpoly_min(int a, int b){ return a<b?a:b; }

int cpoly_clampi(int v, int m, int M){ return (v<m)?m:(v>M?M:v); }

// Result is Polygon0 - Polygon1
int cpoly_cv_diff(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1)
{
  const int npts[2]={npts0,npts1};
  void* pts[2]={pts0,pts1};
  const int stride[2]={stride0,stride1<=0?stride0:stride1};
  int P,v,e,oP,nextv;
  float minx=FLT_MAX;
  int minv;
  float x0,y0,x1,y1,ix,iy;
  int counts[2]={0,0}; // in P0 it counts no. of remaining outsiders, in P1 it marks no. of remaining insiders
  int inside;
  cpolyBitPool bitsets[2]; // in P0 it marks insiders, in P1 it marks outsiders

  cpoly_pool_vcount = 0;
  cpoly_pool_icount = 0;
  // some more checks here
  if (stride1<=0) stride1=stride0;
  cpoly_bitset_create(bitsets,npts0);
  cpoly_bitset_create(bitsets+1,npts1);

  // pick start vertex of first polygon, count for both
  for (P=0;P<2;++P)
  {
    oP = (P+1)%2;
    for (v=0;v<npts[P];++v )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);

      inside = cpoly_cv_point_inside(pts[oP],npts[oP],stride[oP],x0,y0);
      if ( P==0 )
      {
        if ( !inside ){ ++counts[0]; if ( x0<minx ) { minv=v; minx=x0; } }
        cpoly_bitset_set(bitsets,v,inside);
      }
      else 
      {
        if (inside) {++counts[1];} 
        cpoly_bitset_set(bitsets+1,v,inside^1);
      }      
    }
  }
 
  if ( counts[0]==npts[0] && counts[1]==0 ) 
    return 0; // no overlap
   
  // start with that min vertex from P0. In poly0, we go CW, in poly1 we go CCW
  v = minv; nextv=(v+1)%npts[0]; P=0; oP=1; 
  
 
  // while still outside points from P0 and inside points from P1
  while ( counts[0]>0 || counts[1]>0 )
  {
    if ( !cpoly_bitset_get(bitsets+P,v) )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);
      cpoly_pool_add_vertex(x0,y0); --counts[P]; // adds the first one in the edge
      cpoly_bitset_set(bitsets+P,v,1); // marks this vertex to further ignored

      cpoly_getxy(pts[P],stride[P],nextv,x1,y1);

againray:
      // for this edge, intersects with the any of the other polygon's edge?
      if ( cpoly_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts[oP],npts[oP],stride[oP],&ix,&iy,&e) )
      {
        cpoly_pool_add_vertex(ix,iy);        
        P=oP; oP=(oP+1)%2; // the other poly
        //continues with next poly vertex
        if (P==0) v=e; 
        else
        { 
          v=(e+1)%npts[1]; nextv=v-1; if (nextv<0) nextv=npts[1]-1;
          if ( cpoly_bitset_get(bitsets+1,v) && cpoly_bitset_get(bitsets+1,nextv) )
          {
            cpoly_getxy(pts[1],stride[1],nextv,x1,y1);
            x0=ix+(x1-x0)*0.00001f; y0=iy+(y1-y0)*0.00001f;
            oP=0;P=1;
            goto againray;
          }
        }
      }      
    } 
    else if ( !counts[P] ) // no more vertices to process in current polygon, jump to the other
    {
      P=oP; oP=(oP+1)%2;
    }

    if (P==0) { v = (v+1)%npts[0]; nextv=(v+1)%npts[P]; } 
    else      { --v; if (v<0) v=npts[1]-1; nextv=v-1; if (nextv<0) nextv=npts[1]-1; }

    
    // we reach first vertex again for the polygon0, let's check if there's
    if ( P==0 && v==minv && counts[0])
    {
      do
      {
        v=(v+1)%npts[0]; nextv=(v+1)%npts[0];
      }while( cpoly_bitset_get(bitsets,v) && v!=minv);
      
      if (v!=minv)
      {
        minv=v;
        cpoly_pool_add_index(cpoly_pool_vcount);
      }
    }
  }

  if ( cpoly_pool_vcount )
    cpoly_pool_add_index(cpoly_pool_vcount);

  cpoly_bitset_destroy(bitsets+1);
  cpoly_bitset_destroy(bitsets);
  return cpoly_pool_icount;
}

// some basic hom transformation
void cpoly_transform_rotate(void* pts, int npts, int stride, float anglerad, float* xpivot, float* ypivot)
{
  int i;
  float x,y;
  float xp=0,yp=0;
  float a, b;
  const float c=cosf(anglerad);
  const float s=sinf(anglerad);
  
  if ( !xpivot || !ypivot ) cpoly_poly_centroid(pts,npts,stride,&xp,&yp);
  if ( xpivot ) xp=*xpivot;
  if ( ypivot ) yp=*ypivot;

  for (i=0;i<npts;++i)
  {
    cpoly_getxy(pts,stride,i,x,y); 
    a=x-xp; b=y-yp; // translates to origin,
    x = (a*c+b*s)+xp; y = (b*c-a*s)+yp; //  rotate, back to original pos
    cpoly_setxy(pts,stride,i,x,y); 
  }
}

void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot)
{
  int i;
  float x,y;
  float xp=0,yp=0;
  float a, b;
  
  if ( !xpivot || !ypivot ) cpoly_poly_centroid(pts,npts,stride,&xp,&yp);
  if ( xpivot ) xp=*xpivot;
  if ( ypivot ) yp=*ypivot;

  for (i=0;i<npts;++i)
  {
    cpoly_getxy(pts,stride,i,x,y);
    a=x-xp; b=y-yp; // translates to origin,
    x = a*sx+xp; y = b*sy+yp; //  rotate, back to original pos
    cpoly_setxy(pts,stride,i,x,y); 
  }
}

void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y, float* xpivot, float* ypivot)
{
  int i;
  float xp,yp;
  float xc=0,yc=0;

  if ( !xpivot || !ypivot ) cpoly_poly_centroid(pts,npts,stride,&xc,&yc);
  if ( xpivot ) xc=*xpivot;
  if ( ypivot ) yc=*ypivot;

  for (i=0;i<npts;++i)
  {
    cpoly_getxy(pts,stride,i,xp,yp); 
    xp = (xp-xc)+x;
    yp = (yp-yc)+y;
    cpoly_setxy(pts,stride,i,xp,yp);
  }
}

void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy)
{
  int i;
  float ax=0, ay=0;
  float inv;

  for (i=0;i<npts;++i) 
  { 
    ax+=cpoly_getx(pts,stride,i); ay+=cpoly_gety(pts,stride,i); 
  }
  inv= 1.0f / npts;

  if ( cx ) *cx = ax*inv;
  if ( cy ) *cy = ay*inv;
}

// computes convex hull of a polygon. Returns no of indices to vertices in original polygon, use cpoly_pool_get_index
// gift wrapping algorithm
int cpoly_convex_hull(void* pts, int npts, int stride)
{
  int i, endp;
  int initv=0, v;
  float minx=FLT_MAX;
  float x0,y0,x1,y1,x2,y2;

  if ( npts <= 2 ) return 0;

  // init
  cpoly_pool_icount=0;

  // picks up the leftmost to start with
  for (i=0;i<npts;++i)
  {
    x0=cpoly_getx(pts,stride,i); 
    if ( x0<minx ){ minx=x0; initv=i; }
  }
  
  v = endp = initv;
  do
  {
    v=endp;
    cpoly_pool_add_index(v); // add this vertex to CH
    cpoly_getxy(pts,stride,v,x0,y0);

    // first guess
    endp=0;
    cpoly_getxy(pts,stride,endp,x1,y1);
    for (i=1;i<npts;++i)
    {
      cpoly_getxy(pts,stride,i,x2,y2);
      if ( endp==v || cpoly_zcross(x0,y0,x1,y1,x2,y2) > 0.0f )
      {
        // start and end points are the same, advance to next || it's on the left of edge
        endp=i;
        cpoly_getxy(pts,stride,endp,x1,y1);
      }
    }
  }while (initv!=endp);

  return cpoly_pool_icount;
}

void cpoly_msg_create(cpolyBitPool* rmg, int width, int height)
{
  rmg->stride = (int)ceilf(width/2.0f);
  rmg->values = (unsigned char*)calloc(rmg->stride*height,1);
}

void cpoly_msg_destroy(cpolyBitPool* rmg)
{
  if ( rmg->values ) free(rmg->values);
}

int cpoly_msg_get(cpolyBitPool* rmg, int i, int j)
{
  const int byteoffs = rmg->stride*j+(i>>1); // byte offset in memory
  const unsigned char b = rmg->values[byteoffs]; // byte with two values
  const int magic = (1-(i%2))<<2; // avoid branch, offset or not
  return (int)((b>>magic)&0xf); // only takes 4 lsb
}

void cpoly_msg_set(cpolyBitPool* rmg, int i, int j, int val)
{
  const int byteoffs = rmg->stride*j+(i>>1); // byte offset in memory
  const unsigned char b = rmg->values[byteoffs]; // byte with two values
  const char valc = val&0xf;
  rmg->values[byteoffs]= ((i%2)==1) ? ((b&0xf0)|valc) : ((b&0x0f)|(valc<<4));
}

float cpoly_msg_func(void* pts, int npts, int stride, float x, float y)
{
  float f=0.0f;
  int k;
  float x1,y1,r;
  float a,b;

  for (k=0;k<npts;++k) 
  { 
    cpoly_getxyz(pts,stride,k,x1,y1,r); 
    a=x1-x; a*=a; 
    b=y1-y; b*=b; 
    f+= (r*r)/(a+b); 
  }

  return f;
}

int cpoly_msg_cellvalue(float c0, float c1, float c2, float c3)
{
  int val=0;
  if (c0>=1.0f) val |= 1<<0;
  if (c1>=1.0f) val |= 1<<1;
  if (c2>=1.0f) val |= 1<<2;
  if (c3>=1.0f) val |= 1<<3;
  return val;
}

// computes marching squares of the points given by the circles (x0,y0,r0), (x1,y1,r1)...
int cpoly_marchingsq_nointerp(void* pts, int npts, int stride, float sqside)
{
  int i,j,k,mini,minj,dir;
  int stepsx, stepsy;
  float x0,y0,r;
  float x1,y1,x2,y2;
  float c0,c1,c2,c3;
  float minis[2]={FLT_MAX,FLT_MAX};
  float maxis[2]={-FLT_MAX,-FLT_MAX};
  cpolyBitPool grid;
  float *rowvalues;

  if ( npts <= 0 ) return 0;

  // computes AABB of circles
  for (i=0;i<npts;++i)
  {
    cpoly_getxyz(pts,stride,i,x0,y0,r);
    x1=x0-r; y1=y0-r; if ( x1 < minis[0] ) minis[0]=x1; if ( y1 < minis[1] ) minis[1]=y1;
    x1=x0+r; y1=y0+r; if ( x1 > maxis[0] ) maxis[0]=x1; if ( y1 > maxis[1] ) maxis[1]=y1;
  }
  // expanding one row/col on both sides
  minis[0]-=sqside; maxis[0]+=sqside;
  maxis[1]+=sqside; minis[1]-=sqside;

  // setting up
  stepsx = (int)ceilf((maxis[0]-minis[0])/sqside);
  stepsy = (int)ceilf((maxis[1]-minis[1])/sqside);
  maxis[0] = minis[0]+sqside*stepsx;
  minis[1] = maxis[1]-sqside*stepsy;

  // grid will store the values, 4 bits per cell
  cpoly_msg_create(&grid,stepsx,stepsy);  
  rowvalues = (float*)cpoly_pool_i; 
  mini=minj=-1;

  // computing cell values
  {
    // first row (avoid duplicate computation)
    x0=minis[0]; y0=maxis[1]; 
    for (i=0;i<=stepsx;++i,x0+=sqside) rowvalues[i]=cpoly_msg_func(pts,npts,stride,x0,y0);
    
    // computing cell corners values for the rest of rows
    y0 = maxis[1]-sqside;
    for (j=0;j<stepsy;++j,y0-=sqside)
    {
      x0 = minis[0]+sqside;
      for (i=0;i<stepsx;++i,x0+=sqside)
      {
        c0=cpoly_msg_func(pts,npts,stride,x0-sqside,y0);
        c1=cpoly_msg_func(pts,npts,stride,x0,y0);
        c2=rowvalues[i+1];
        c3=rowvalues[i];
        k=cpoly_msg_cellvalue(c0,c1,c2,c3);
        if ( k==15 )
        {
          if      ( i==0 ) k=6;
          else if ( i==stepsx-1 ) k=9;
          else if ( j==0 ) k=3;
          else if ( j==stepsy-1 ) k=12;
        }
        cpoly_msg_set(&grid,i,j,k);
        if ( k && mini==-1){ mini=i; minj=j; } // finding first valid cell from top,left
        rowvalues[i]=c0;
      }
      rowvalues[i]=c1;
    }
  }

  // adding poly vertices
  cpoly_pool_vcount=cpoly_pool_icount=0;
  c0=sqside*0.5f;

  // outlining a polygon (clock wise)
  do 
  {
    i=mini; j=minj;
    dir=0;
    do 
    {
      k=cpoly_msg_get(&grid,i,j);
      if ( k!=5 && k!=10 ) cpoly_msg_set(&grid,i,j,0);
      x0=minis[0]+sqside*i; x1=x0+c0; x2=x0+sqside; // left, med, right (HORIZ)
      y0=maxis[1]-sqside*j; y1=y0-c0; y2=y0-sqside; // top, med, bottom (VERTI)
      switch ( k )
      {
      case 1: case 13: case 9:    cpoly_pool_add_vertex(x1, y2); ++j; break;
      case 2: case 3: case 11:    cpoly_pool_add_vertex(x2, y1); ++i; break;
      case 6: case 7: case 4:     cpoly_pool_add_vertex(x1, y0); --j; break;
      case 8: case 12: case 14:   cpoly_pool_add_vertex(x0, y1); --i; break;      
      case 5:
        switch (dir)
        { 
        case 8: case 12: case 14: cpoly_pool_add_vertex(x1, y2); ++j; break;
        default:                  cpoly_pool_add_vertex(x1, y0); --j; break;
        }
      break;
      case 10: 
        switch (dir)
        {
        case 1: case 13: case 9:  cpoly_pool_add_vertex(x2,y1); ++i; break;
        default:                  cpoly_pool_add_vertex(x0,y1); --i; break;
        }
      break;
      default: 
        goto exitPoly;
      }

      // advance CW, caution with boundaries
      do
      {
        dir=0;
        if ( i<0 )      { i=0; --j; dir=1;}
        if ( i>=stepsx ){ i=stepsx-1; ++j; dir=1;}
        if ( j<0 )      { j=0; ++i; dir=1; }
        if ( j>=stepsy ){ j=stepsy-1; --i; dir=1; }
      }while(dir);

      dir = k; // last k in dir
    } while (i!=mini ||j!=minj); // until we reach first cell again
exitPoly:
    cpoly_pool_add_index(cpoly_pool_vcount);

    // find other polygon
    {
      j=minj;
      minj=mini=-1;
      for (;j<stepsy && minj==-1;++j)
      {
        for (i=0;i<stepsx;++i)
        {
          k=cpoly_msg_get(&grid,i,j);
          if ( k>0 && k<15 && k!=5 && k!=10)
          {
            mini=i; minj=j;
            break;
          }
        }
      }
    }
  } while ( mini!=-1 && minj!=-1 );
  cpoly_msg_destroy(&grid);
  return cpoly_pool_icount;
}

// computes axis aligned bounding box
void cpoly_aabb(void* pts, int npts, int stride, float* xmin, float* ymin, float* xmax, float* ymax)
{
  int i;
  float x,y;

  *xmin = *ymin = FLT_MAX;
  *xmax = *ymax = -FLT_MAX;

  // computes AABB 
  for (i=0;i<npts;++i)
  {
    cpoly_getxy(pts,stride,i,x,y);
    if ( x < *xmin ) *xmin=x; if ( y < *ymin ) *ymin=y;
    if ( x > *xmax ) *xmax=x; if ( y > *ymax ) *ymax=y;
  }
}

// computes convex partition with Hertel-Mehlhorn algorithm
int cpoly_convex_partition_HM(void* pts, int npts, int stride)
{

  return 0;
}

void cpoly_pvs_create(cpolyBitPool* pvs, int npts)
{
  const int nbits = (int)ceilf(((float)npts-1.0f)*0.5f*npts);
  cpoly_bitset_create(pvs,nbits);
  pvs->stride = npts;
}

void cpoly_pvs_destroy(cpolyBitPool* pvs)
{
  pvs->stride=0;
  free(pvs->values);
}

int cpoly_pvs_get(cpolyBitPool* pvs, int i, int j)
{
  int ndx;
#ifdef _DEBUG
  const int maxoffs=(int)( pvs->stride*( (pvs->stride-1)*0.5f ) );
#endif
  // swap to keep i < j
  if ( i>j ){ ndx=i; i=j; j=ndx; } 

  // adjacencies
  if ( j==((i+1)%pvs->stride) ) return 1;
  ndx = i-1; if (ndx<0) ndx=pvs->stride-1;
  if ( j==ndx ) return 1;

  // offset
  ndx = (int)( i * ( (2*pvs->stride - i - 3)*0.5f ) + j - 1 );
#ifdef _DEBUG
  if (ndx < 0 || ndx>=maxoffs) return 0;
#endif
  return cpoly_bitset_get(pvs,ndx);
}

void cpoly_pvs_set(cpolyBitPool* pvs, int i, int j, int v)
{
  const int offset = (int)( i * ( (2*pvs->stride - i - 3)*0.5f ) + j - 1 );
#ifdef _DEBUG
  const int maxoffs=(int)( pvs->stride*( (pvs->stride-1)*0.5f ) );
  if (offset < 0 || offset>=maxoffs) return;
#endif  
  cpoly_bitset_set(pvs,offset,v);
}

// check if ray from vertex i to vertex j is an inner ray (going from inside of the polygon)
int cpoly_is_inner_ray(void* pts, int npts, int stride, int i, int j)
{
  int iprev,inext; // index to previous and next vertices to i
  float xi,yi,xp,yp,xn,yn; // from i, ia, ib respectively
  float xj,yj; // from j
  float z,zprev,znext;
  int r =0;

  iprev = i-1; if (iprev<0) iprev=npts-1;
  inext = (i+1)%npts;
  cpoly_getxy(pts,stride,iprev, xp,yp);
  cpoly_getxy(pts,stride,inext, xn,yn);
  cpoly_getxy(pts,stride,i, xi,yi);
  cpoly_getxy(pts,stride,j, xj,yj);

  // check if i is a convex or a concave vertex (clock wise)
  z=cpoly_zcross(xp,yp,xi,yi,xn,yn);

  zprev=cpoly_zcross(xi,yi,xp,yp,xj,yj);
  znext=cpoly_zcross(xi,yi,xn,yn,xj,yj);
  if ( z <= 0.0f )
  {
    r = (zprev*znext)<=0.0f ? 1 : 0; // convex vertex. j should be in different sides of edges
  }
  else
  {
    // concave vertex
    if ( zprev<=0.0f && znext>=0.0f ) r=0;
    else r=1;
  }
  return r;
}

// compute potentially visible sets
void cpoly_pvs(void* pts, int npts, int stride, cpolyBitPool* pvs)
{
  int i,j,r;
  float x0,y0,x1,y1;

  g_rayOffset=0.001f;
  for (i=0;i<npts; ++i )
  {
    // for my adjacent point, vis=1
    cpoly_pvs_set(pvs,i,i+1,1);
    cpoly_getxy(pts,stride,i,x0,y0);
    for (j=i+2;j<npts;++j)
    {
      if ( i==0 && j==npts-1 ) { cpoly_pvs_set(pvs,i,j,1); continue; } // special case for first vertex adjacent with last one (loop) (don't like this branch here)
      if ( cpoly_is_inner_ray(pts,npts,stride, i,j) )
      {
        cpoly_getxy(pts,stride,j,x1,y1);
        r = cpoly_cv_seg_isec_poly_first(x0,y0,x1,y1,pts,npts,stride,NULL,NULL,NULL);
        cpoly_pvs_set(pvs,i,j,1-r); // vis=1 if no intersection happened (r=0) or vis=0 when r=1
      }
      else
      {
        cpoly_pvs_set(pvs,i,j,0); // no inner ray, vis=0
      }
    }
  }
  g_rayOffset=0.0f;
}

int cpoly_convex_partition(void* pts, int npts, int stride)
{
  int i, j, k, startp;
  int n,m,v,nextp,initi;
  float x0,y0,x1,y1,x2,y2;
  float z;
  cpolyBitPool pvs;

  if ( npts < 3 ) return 0;

  cpoly_pvs_create(&pvs, npts);
  cpoly_pvs(pts,npts,stride,&pvs);

  cpoly_pool_icount=0;
  startp=cpoly_pool_icount;

  i=0; j=1; k=2;
  initi=i;
  cpoly_pool_add_index(i);
  cpoly_pool_add_index(j);
  nextp=-1;
  while (k!=initi)
  {
    // k is visible for all vertices in our current convex poly?
    v=1;
    for (n=startp;n<cpoly_pool_icount && v;++n)
    {
      m = cpoly_pool_get_index(n);
      v=cpoly_pvs_get(&pvs, m, k);
    }
    
    // if it's visible, does it break the convexity if we add it?
    if (v)
    {
      cpoly_getxy(pts,stride,i,x0,y0);
      cpoly_getxy(pts,stride,j,x1,y1);
      cpoly_getxy(pts,stride,k,x2,y2);
      z=cpoly_zcross(x0,y0,x1,y1,x2,y2);
      if ( z<=0.0f )
      {
        cpoly_pool_add_index(k);
        i=j; j=k;
      }
      else goto jumpvertex; // k was visible, but breaks convexity
    }
    else
    {
      jumpvertex:
      // k wasn't visible
      if ( nextp==-1 ) nextp=k;
    }
    k=(k+1)%npts;
  }

  /*
  V = Compute Visibility From Each Vertex to other Vertices
  CP = {}

  while ( k!=i )
  {
    CP += {i,j}
    Check if: 
        a) i,j and j,k are convex
        b) no intersection from all inner rays CP -> k
    If ( a && b )
        CP += k 
        i=j, j=k
    k = next vertex
  }
      
      

  */

  cpoly_pvs_destroy(&pvs);
  return 0;
}



#endif
