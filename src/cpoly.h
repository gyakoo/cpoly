/*

 */

#ifndef CPOLY_H
#define CPOLY_H

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  // Partition a polygon into a set of convex polygons (ClockWise)
  int cpoly_cv_partitioning_cw(void* pts, int npts, int stride, int** partndxs, int** poffsets);
  
  // Deallocates memory previously allocated by cpoly_decomp* functions
  void cpoly_free_parts(int** parts, int** psizes);

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

  // union of convex polygons
  // assumes: convex polygons, intersecting, no one inside another, no holes.
  int cpoly_cv_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

  // some basic hom transformation
  void cpoly_transform_rotate(void* pts, int npts, int stride, float angle, float* xpivot, float* ypivot);
  void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot);
  void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y);

  // geometric center computation (center of mass or centroid)
  void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy);

#ifdef __cplusplus
};
#endif

#endif // CPOLY_H

#ifdef CPOLY_IMPLEMENTATION

#ifndef CPOLY_MAXPOOLSIZE_BYTES 
#define CPOLY_MAXPOOLSIZE_BYTES (16<<10)
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

// non thread safe shared pool
unsigned char cpoly_pool[CPOLY_MAXPOOLSIZE_BYTES];
const int cpoly_poolflts = CPOLY_MAXPOOLSIZE_BYTES/sizeof(float); // max no of floats in the pool
const int cpoly_poolverts = CPOLY_MAXPOOLSIZE_BYTES/(sizeof(float)*2); // max no of vertices (x,y) in the pool
int cpoly_pool_count=0;

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

#define cpoly_setx(p,s,n,f) cpoly_getx(p,s,n)=f
#define cpoly_sety(p,s,n,f) cpoly_gety(p,s,n)=f


int cpoly_cv_partitioning_cw(void* pts, int npts, int stride, int** partndxs, int** poffsets)
{
  return 0;
  /*
  LOOP1
  if first partition, 
  select current point and next one (cp=current, np=next)
  else
  select prior and next one (cp=prior, np=current)
  P += {cp,np}

  LOOP2
  compute plane of last edge (cp-->np)
  p = get next point to np positive to plane
  if p no intersect (cp-->p) with any edge
  add p to partition P
  sp=np, np=p
  Goto LOOP2 while points

  add partition P
  add adjacent points to remaining points
  Goto LOOP1 while points (if 1+2 remaining points is last partition though)

   */
}

void cpoly_free_parts(int** partndxs, int** poffsets)
{
//  int i=0; 

  if ( !partndxs || !*partndxs || !poffsets || !*poffsets )
    return;
  free(*partndxs);
  free(*poffsets);
  *partndxs = *poffsets = 0;
}

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
    x0 = cpoly_getx(pts,stride,i); y0 = cpoly_gety(pts,stride,i);
    x1 = cpoly_getx(pts,stride,i+1); y1 = cpoly_gety(pts,stride,i+1);
    x2 = cpoly_getx(pts,stride,(i+2)%npts); y2 = cpoly_gety(pts,stride,(i+2)%npts);

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
  _x0= x0 = cpoly_getx(pts,stride,0); _y0= y0 = cpoly_gety(pts,stride,0);
  x1      = cpoly_getx(pts,stride,1); y1      = cpoly_gety(pts,stride,1);
  x2      = cpoly_getx(pts,stride,2); y2      = cpoly_gety(pts,stride,2);
  zcross = cpoly_zcross(x0,y0,x1,y1,x2,y2);

  // we have already the first three points
  curz = cpoly_zcross(x0,y0,x1,y1,x,y); if ( (curz*zcross) < 0.0f ) return 0;
  curz = cpoly_zcross(x1,y1,x2,y2,x,y); if ( (curz*zcross) < 0.0f ) return 0;

  // iterate all remaining points
  x0=x2; y0=y2;
  for ( i=3; i<npts;++i)
  {
    x1 = cpoly_getx(pts,stride,i); y1 = cpoly_gety(pts,stride,i);
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
      x0=cpoly_getx(pts[P],stride[P],e); y0=cpoly_gety(pts[P],stride[P],e);
      x1=cpoly_getx(pts[P],stride[P],(e+1)%npts[P]); y1=cpoly_gety(pts[P],stride[P],(e+1)%npts[P]);
      // edge normal
      nx = y1-y0; ny = -(x1-x0);
      
      // projecting all points from P0 and P1 onto normal and getting min and maxs
      minis[0]=minis[1]=FLT_MAX; maxis[0]=maxis[1]=-FLT_MAX;
      for (p=0;p<2;++p)
      {
        for (v=0;v<npts[p];++v)
        {
          xp=cpoly_getx(pts[p],stride[p],v); yp=cpoly_gety(pts[p],stride[p],v);
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


int cpoly_seg_isec(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3, float* is)
{
  const float a=x2-x0; const float b=x3-x2;
  const float c=x1-x0; const float d=y3-y2;
  const float e=y2-y0; const float f=y1-y0;
  const float t =  (c*e - a*f) / (b*f - d*c);
  const float s = (1.0f/c)*( t*b + a);

  // if are indeterminate, they're collinear
  // if are infinite, they're parallel
  // in both cases no inters and will return 0
  // segments, so they'd be in the range
  if ( s>=0.0f && s<=1.0f && t>=0.0f && t<=1.0f )
  {
    if ( is ) *is= s;
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
    x2=cpoly_getx(pts,stride,v);          y2=cpoly_gety(pts,stride,v);
    x3=cpoly_getx(pts,stride,(v+1)%npts); y3=cpoly_gety(pts,stride,(v+1)%npts);
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

void cpoly_pool_add(float x, float y)
{
  float* ptr = (float*)( cpoly_pool+cpoly_pool_count*sizeof(float)*2);
  *(ptr++) = x; *ptr = y;
  ++cpoly_pool_count;
}

void cpoly_pool_get(int n, float* x, float* y)
{
  float* ptr = (float*)( cpoly_pool+n*sizeof(float)*2);
  *x=*(ptr++); *y=*ptr;
}

void cpoly_bitset_reset(cpoly_bitset bs, int n)
{
  bs[n>>3] &= ~( 1<<(n%8) );
}

void cpoly_bitset_set(cpoly_bitset bs, int n)
{
  bs[n>>3] |= 1<<(n%8);
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

// assumes: convex polygons, no one inside another, no holes.
int cpoly_cv_union(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1)
{
  const int npts[2]={npts0,npts1};
  void* pts[2]={pts0,pts1};
  const int stride[2]={stride0,stride1<=0?stride0:stride1};
  int P,v,e,oP;
  float minx=FLT_MAX;
  int minv, minp;
  float x0,y0,x1,y1,ix,iy,s;
  int counts[2]={0,0};
  int inside;
  cpoly_bitset bitsets[2];

  // some more checks here
  if (stride1<=0) stride1=stride0;
  if (npts0>=cpoly_bitset_maxbits || npts1>=cpoly_bitset_maxbits) return 0;
  
  // pick start vertex out of two polygons (can be the one with min X if we use a straight Vertical sweep line)
  for (P=0;P<2;++P)
  {
    oP = (P+1)%2;
    for (v=0;v<npts[P];++v )
    {
      x0 = cpoly_getx(pts[P],stride[P],v); y0 = cpoly_gety(pts[P],stride[P],v);
      if ( x0<minx ) { minp=P; minv=v; minx=x0; }

      // computes how many points of P are outside of oP
      // mark flag 1 for those points that are in the overlapping region
      inside = cpoly_cv_point_inside(pts[oP],npts[oP],stride[oP],x0,y0);
      if ( !inside ) { ++counts[P]; cpoly_bitset_reset(bitsets[P],v); } else cpoly_bitset_set(bitsets[P],v);
    }
  }

  if ( counts[0]==npts[0] && counts[1]==npts[1] ) return 0; // no overlap
  // start with that min vertex
  v = minv; P = minp; oP=(P+1)%2;
  cpoly_pool_count = 0;

  // while still outside points to process for both
  while ( counts[0] || counts[1] )
  {
    if ( !cpoly_bitset_get(bitsets[P],v) )
    {
      x0 = cpoly_getx(pts[P],stride[P],v); y0 = cpoly_gety(pts[P],stride[P],v);
      cpoly_pool_add(x0,y0); // adds the first one in the edge
      --counts[P];
      x1 = cpoly_getx(pts[P],stride[P],(v+1)%npts[P]); y1 = cpoly_gety(pts[P],stride[P],(v+1)%npts[P]);

      // for this edge, intersects with the any of the other polygon's edge?
      if ( cpoly_cv_seg_isec_poly_closest(x0,y0,x1,y1,pts[oP],npts[oP],stride[oP],&ix,&iy,&e) )
      {
        cpoly_pool_add(ix,iy);
        P=oP; oP=(oP+1)%2; // the other poly
        v=e; //continues with next poly vertex, the end of the edge
      }      
    }    
    v = (v+1)%npts[P];
  }

  return cpoly_pool_count;
}

// some basic hom transformation
void cpoly_transform_rotate(void* pts, int npts, int stride, float anglerad, float* xpivot, float* ypivot)
{
  int i;
  float x,y;
  float xp,yp;
  const float c=cosf(anglerad);
  const float s=sinf(anglerad);
  float hx, hy;

  if ( !xpivot || !ypivot )
  {
    cpoly_poly_centroid(pts,npts,stride,&xp,&yp);
    if ( xpivot ) xp=*xpivot;
    if ( ypivot ) yp=*ypivot;
  }

  for (i=0;i<npts;++i)
  {
    x=cpoly_getx(pts,stride,i); y=cpoly_gety(pts,stride,i);
    x-=xp; y-=yp;
    x= x*c+y*s; y= -x*s+y*c;
    x+=xp; y+=yp;
    cpoly_setx(pts,stride,i,x); cpoly_sety(pts,stride,i,y);
  }
}

void cpoly_transform_scale(void* pts, int npts, int stride, float sx, float sy, float* xpivot, float* ypivot)
{
}

void cpoly_transform_translate(void* pts, int npts, int stride, float x, float y)
{
}

void cpoly_poly_centroid(void* pts, int npts, int stride, float* cx, float* cy)
{
  int i;
  float ax=0, ay=0;
  float inv;

  for (i=0;i<npts;++i) ax+=cpoly_getx(pts,stride,i); ay=cpoly_gety(pts,stride,i);
  inv= 1.0f / npts;

  if ( cx ) *cx = ax*inv;
  if ( cy ) *cy = ay*inv;
}


#endif
