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

  // returns 1 if segment intersects polygon
  // returns in ix,iy the closest intersection point
  // edge is the edge index in the polygon assuming they're consecutive (edge=start vertex)
  int cpoly_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // returns 1 if segment intersects polygon
  int cpoly_cv_seg_isec_poly_first(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge);

  // Union of convex polygons. Returns no. of vertices to be accessed with cpoly_pool_get_vertex
  // assumes: convex polygons, intersecting, no one inside another, no holes.
  int cpoly_cv_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);
  
  // Difference of convex polygons, returns no. of parts generated
  // You got the parts offsets starting from part 1 in cpoly_pool_i indices
  // All the vertices generated in cpoly_pool_v
  /* example:
    cpoly_cv_diff(poly0, N0, STRIDE, poly1, N1, STRIDE);    
    k=0;
    for ( i = 0; i < cpoly_pool_icount; ++i ) // for all parts
    {
      glBegin(GL_LINE_LOOP);
      for ( j=k; j<cpoly_pool_get_index(i); ++j )
      {
          cpoly_pool_get_vertex(j,&x,&y);
          glVertex2f( x, y);
      }
      glEnd();
      k=j;
  */
  int cpoly_cv_diff(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // some basic homogeneous transformation functions
  void cpoly_transform_rotate(void* pts, int npts, int stride, float angle, float* xpivot, float* ypivot);
  void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot);
  void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y, float* xpivot, float* ypivot);

  // geometric center computation (center of mass or centroid)
  void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy);

  // computes convex hull of a polygon. Returns no of indices to vertices in original polygon, use cpoly_pool_get_index
  int cpoly_convex_hull(void* pts, int npts, int stride);

  // computes axis aligned bounding box
  void cpoly_aabb(void* pts, int npts, int stride, float* xmin, float* ymin, float* xmax, float* ymax);

  // 
  int cpoly_marching_sq(void* pts, int npts, int stride, float sqside);

#ifdef __cplusplus
};
#endif

#endif // CPOLY_H

#ifdef CPOLY_IMPLEMENTATION
/*
TODO: 
  http://jamie-wong.com/2014/08/19/metaballs-and-marching-squares/
  http://www.dma.fi.upm.es/docencia/trabajosfindecarrera/programas/geometriacomputacional/PiezasConvex/algoritmo_i.html
*/

#ifndef CPOLY_MAXPOOLSIZE_BYTES 
#define CPOLY_MAXPOOLSIZE_BYTES (1<<20)
#else
#if CPOLY_MAXPOOLSIZE_BYTES <= 32
#undef CPOLY_MAXPOOLSIZE_BYTES
#define CPOLY_MAXPOOLSIZE_BYTES (16<<10)
#endif
#endif

#ifndef CPOLY_MAXBITSET_BYTES
#define CPOLY_MAXBITSET_BYTES (128)
#else
#if CPOLY_MAXBITSET_BYTES <= 8
#undef CPOLY_MAXBITSET_BYTES
#define CPOLY_MAXBITSET_BYTES (128)
#endif
#endif

// non thread safe shared pools
unsigned char cpoly_pool_v[CPOLY_MAXPOOLSIZE_BYTES];
unsigned char cpoly_pool_i[CPOLY_MAXPOOLSIZE_BYTES];
const int CPOLY_POOL_VMAX = CPOLY_MAXPOOLSIZE_BYTES/(sizeof(float)*2); // max no of vertices (x,y) in the pool
const int CPOLY_POOL_IMAX= CPOLY_MAXPOOLSIZE_BYTES/sizeof(int); // max no of indices
int cpoly_pool_vcount=0;
int cpoly_pool_icount=0;

// bitset type
typedef char cpoly_bitset[CPOLY_MAXBITSET_BYTES];
const int cpoly_bitset_maxbits = CPOLY_MAXBITSET_BYTES*8;

#ifdef _MSC_VER
	#pragma warning (disable: 4996) // Switch off security warnings
  #pragma warning (disable: 4100) // Switch off unreferenced formal parameter warnings
  #pragma warning (disable: 4204) // Switch off nonstandard extension used : non-constant aggregate initializer
#endif

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

// forward decls
int cpoly_max(int a, int b);
int cpoly_min(int a, int b);
void cpoly_pool_add_vertex(float x, float y);
void cpoly_pool_add_index(int ndx);
void cpoly_pool_get_vertex(int n, float* x, float* y);
int cpoly_pool_get_index(int n);
void cpoly_bitset_put(cpoly_bitset bs, int n, int value);
int cpoly_bitset_get(cpoly_bitset bs, int n);
void cpoly_bitset_allzero(cpoly_bitset bs);

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

int cpoly_cv_seg_isec_poly_closest(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge)
{
  return __internal_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts,npts,stride,ix,iy,edge,0);
}

int cpoly_cv_seg_isec_poly_first(float x0, float y0, float x1, float y1, void* pts, int npts, int stride, float* ix, float* iy, int* edge)
{
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

void cpoly_bitset_put(cpoly_bitset bs, int n, int value)
{
  bs[n>>3] ^= value<<(n%8);
}

int cpoly_bitset_get(cpoly_bitset bs, int n)
{
  const int nmod8=n%8;
  return (bs[n>>3] & (1<<nmod8)) >> nmod8;
}

void cpoly_bitset_allzero(cpoly_bitset bs)
{
  int i;
  for (i=0;i<CPOLY_MAXBITSET_BYTES;++i) bs[i]=0;
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
  cpoly_bitset bitsets[2];

  // some more checks here
  if (stride1<=0) stride1=stride0;
  if (npts0>=cpoly_bitset_maxbits || npts1>=cpoly_bitset_maxbits) return 0;

  // pick start vertex out of two polygons (can be the one with min X if we use a straight Vertical sweep line)
  for (P=0;P<2;++P)
  {
    cpoly_bitset_allzero(bitsets[P]);
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
      cpoly_bitset_put(bitsets[P],v, inside);
    }
  }

  if ( counts[0]==npts[0] && counts[1]==npts[1] ) return 0; // no overlap
  // start with that min vertex
  v = minv; P = minp; oP=(P+1)%2;
  cpoly_pool_vcount = 0;

  // while still outside points to process for both
  while ( counts[0]>0 || counts[1]>0 )
  {
    if ( !cpoly_bitset_get(bitsets[P],v) )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);
      cpoly_pool_add_vertex(x0,y0); --counts[P]; // adds the first one in the edge
      cpoly_bitset_put(bitsets[P],v,1); // marks this vertex as outside, so can't be further processed

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

  return cpoly_pool_vcount;
}

int cpoly_max(int a, int b){ return a>b?a:b; }
int cpoly_min(int a, int b){ return a<b?a:b; }

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
  cpoly_bitset bitsets[2]; // in P0 it marks insiders, in P1 it marks outsiders

  cpoly_pool_vcount = 0;
  cpoly_pool_icount = 0;
  // some more checks here
  if (stride1<=0) stride1=stride0;
  if (npts0>=cpoly_bitset_maxbits || npts1>=cpoly_bitset_maxbits) return 0;

  // pick start vertex of first polygon, count for both
  for (P=0;P<2;++P)
  {
    cpoly_bitset_allzero(bitsets[P]);
    oP = (P+1)%2;
    for (v=0;v<npts[P];++v )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);

      inside = cpoly_cv_point_inside(pts[oP],npts[oP],stride[oP],x0,y0);
      if ( P==0 )
      {
        if ( !inside ){ ++counts[0]; if ( x0<minx ) { minv=v; minx=x0; } }
        cpoly_bitset_put(bitsets[0],v,inside);
      }
      else 
      {
        if (inside) {++counts[1];} 
        cpoly_bitset_put(bitsets[1],v,inside^1);
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
    if ( !cpoly_bitset_get(bitsets[P],v) )
    {
      cpoly_getxy(pts[P],stride[P],v,x0,y0);
      cpoly_pool_add_vertex(x0,y0); --counts[P]; // adds the first one in the edge
      cpoly_bitset_put(bitsets[P],v,1); // marks this vertex to further ignored

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
          if ( cpoly_bitset_get(bitsets[1],v) && cpoly_bitset_get(bitsets[1],nextv) )
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
      }while( cpoly_bitset_get(bitsets[0],v) && v!=minv);
      
      if (v!=minv)
      {
        minv=v;
        cpoly_pool_add_index(cpoly_pool_vcount);
      }
    }
  }

  if ( cpoly_pool_vcount )
    cpoly_pool_add_index(cpoly_pool_vcount);

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

// uses bytes, but olny 4 bits per byte are used to represent a grid cell value (0-15)
typedef struct sRMGrid
{
  unsigned char* values;
  int stride;
}sRMGrid;

void cpoly_rmg_create(sRMGrid* rmg, int width, int height)
{
  rmg->stride = (int)ceilf(width/2.0f);
  rmg->values = (unsigned char*)calloc(rmg->stride*height,1);
}

void cpoly_rmg_destroy(sRMGrid* rmg)
{
  if ( rmg->values ) free(rmg->values);
}

int cpoly_rmg_get(sRMGrid* rmg, int i, int j)
{
  const int byteoffs = rmg->stride*j+(i>>1); // byte offset in memory
  const unsigned char b = rmg->values[byteoffs]; // byte with two values
  const int magic = (1-(i%2))<<2; // avoid branch, offset or not
  return (int)((b>>magic)&0xf); // only takes 4 lsb
}

void cpoly_rmg_set(sRMGrid* rmg, int i, int j, int val)
{
  const int byteoffs = rmg->stride*j+(i>>1); // byte offset in memory
  const unsigned char b = rmg->values[byteoffs]; // byte with two values
  const char valc = val&0xf;
  rmg->values[byteoffs]= ((i%2)==1) ? ((b&0xf0)|valc) : ((b&0x0f)|(valc<<4));
}

float cpoly_rmg_func(void* pts, int npts, int stride, float x, float y)
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

int cpoly_rmg_cellvalue(float c0, float c1, float c2, float c3)
{
  int val=0;
  if (c0>=1.0f) val |= 1<<0;
  if (c1>=1.0f) val |= 1<<1;
  if (c2>=1.0f) val |= 1<<2;
  if (c3>=1.0f) val |= 1<<3;
  return val;
}

int cpoly_clampi(int v, int m, int M){ return (v<m)?m:(v>M?M:v); }

// computes marching squares of the points given by the circles (x0,y0,r0), (x1,y1,r1)...
int cpoly_marching_sq(void* pts, int npts, int stride, float sqside)
{
  int i,j,k,mini,minj,dir;
  int stepsx, stepsy;
  float x0,y0,r;
  float x1,y1;
  float c0,c1,c2,c3;
  float minis[2]={FLT_MAX,FLT_MAX};
  float maxis[2]={-FLT_MAX,-FLT_MAX};
  sRMGrid grid;
  float *rowvalues;

  if ( npts <= 0 ) return 0;

  // computes AABB of circles
  for (i=0;i<npts;++i)
  {
    cpoly_getxyz(pts,stride,i,x0,y0,r);
    x1=x0-r; y1=y0-r; if ( x1 < minis[0] ) minis[0]=x1; if ( y1 < minis[1] ) minis[1]=y1;
    x1=x0+r; y1=y0+r; if ( x1 > maxis[0] ) maxis[0]=x1; if ( y1 > maxis[1] ) maxis[1]=y1;
  }

  // setting up
  stepsx = (int)ceilf((maxis[0]-minis[0])/sqside);
  stepsy = (int)ceilf((maxis[1]-minis[1])/sqside);
  cpoly_rmg_create(&grid,stepsx,stepsy);  
  rowvalues = (float*)cpoly_pool_i; 
  mini=minj=-1;
  cpoly_pool_vcount=0;

  // computing cell values
  {
    // first row (avoid duplicate computation)
    x0=minis[0]; y0=maxis[1]; 
    for (i=0;i<=stepsx;++i,x0+=sqside)
    { 
      rowvalues[i]=cpoly_rmg_func(pts,npts,stride,x0,y0); 
//       if ( rowvalues[i]>=1.0f)
//         cpoly_pool_add_vertex(x0,y0);
    }
    
    // computing cell corners values for the rest of rows
    y0 = maxis[1]-sqside;
    for (j=0;j<stepsy;++j,y0-=sqside)
    {
      x0 = minis[0]+sqside;
      for (i=0;i<stepsx;++i,x0+=sqside)
      {
        c0=cpoly_rmg_func(pts,npts,stride,x0-sqside,y0); //if ( c0>=1.0f) cpoly_pool_add_vertex(x0-sqside,y0);
        c1=cpoly_rmg_func(pts,npts,stride,x0,y0);        //if ( c1>=1.0f) cpoly_pool_add_vertex(x0,y0);
        c2=rowvalues[i+1];
        c3=rowvalues[i];
        k=cpoly_rmg_cellvalue(c0,c1,c2,c3);
        cpoly_rmg_set(&grid,i,j,k);
        if ( k && mini==-1)
        { 
          mini=i; minj=j; 
        }
        rowvalues[i]=c0;
      }
      rowvalues[i]=c1;
    }
  }

  // outlining a polygon (clock wise)
  {
    i=mini; j=minj;
    k=cpoly_rmg_get(&grid,i,j);
    if ( k!=2 )
    {
      k=k;
    }
  }

  // temporary
//   y0 = maxis[1]-sqside*0.5f;
//   for (j=0;j<stepsy/2;++j,y0-=sqside)
//   {
//     x0 = minis[0]+sqside*0.5f;
//     for (i=0;i<stepsx;++i,x0+=sqside)
//     {
//       k = cpoly_rmg_get(&grid,i,j);
//       if ( k && k<15 )
//         cpoly_pool_add_vertex(x0,y0);
//     }
//   }
  cpoly_rmg_destroy(&grid);
  return cpoly_pool_vcount;
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

#endif
