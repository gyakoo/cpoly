/*

 */

#ifndef CPOLY_H
#define CPOLY_H

#include <math.h>
#ifdef __cplusplus
extern "C" {
#endif

  // Partition a polygon into a set of convex polygons (ClockWise)
  int cpoly_partitioning_cw(void* pts, int npts, int stride, int** partndxs, int** poffsets);
  // Deallocates memory previously allocated by cpoly_decomp* functions
  void cpoly_free_parts(int** parts, int** psizes);

  // Returns 0 if it's not convex. 1 for CW for CCW
  int cpoly_is_convex(void* pts, int npts, int stride);

  // Returns 1 if the point is inside a convex polygon
  int cpoly_cv_point_inside(void* pts, int npts, int stride, float x, float y);

  // Return 1 if two convex polygons intersects (Separating Axis Theorem)
  int cpoly_cv_intersects_SAT(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

#ifdef __cplusplus
};
#endif

#endif // CPOLY_H

#ifdef CPOLY_IMPLEMENTATION


#ifdef _MSC_VER
	#pragma warning (disable: 4996) // Switch off security warnings
  #pragma warning (disable: 4100) // Switch off unreferenced formal parameter warnings
  #pragma warning (disable: 4204) // Switch off nonstandard extension used : non-constant aggregate initializer
#endif

// NEGATIVE is to the eye (CW) POSITIVE is towards the screen (CCW)
#define cpoly_zcross(x0,y0, x1, y1, x2, y2) ((x1-x0)*(y2-y1) - (y1-y0)*(x2-x1))
#define cpoly_getx(p,s,n) * (    (float*)( (char*)(p) + ((s)*(n)) ) )
#define cpoly_gety(p,s,n) * (    (float*)( (char*)(p) + ((s)*(n)+sizeof(float)) ) )


int cpoly_partitioning_cw(void* pts, int npts, int stride, int** partndxs, int** poffsets)
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
  int i0, i1, i2;
  int count;
  float curz=.0f, lastz=.0f, s;
  float x0,y0,x1,y1,x2,y2;
  char* ptr=(char*)pts, *curptr;
  
  if ( npts <= 2 ) return 0; // not a polygon
  
  count=npts-1;
  for (i0=0, i1=1; i0<count; ++i0, ++i1, ptr+=stride)
  {
    curptr = ptr;
    i2=(i1+1)%npts;

    // z comp of cross product of both edges
    x0 = *(float*)(curptr+0); y0 = *(float*)(curptr+sizeof(float)); curptr+=stride;
    x1 = *(float*)(curptr+0); y1 = *(float*)(curptr+sizeof(float)); curptr+=stride;
    x2 = *(float*)(curptr+0); y2 = *(float*)(curptr+sizeof(float));
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
  char* curptr=(char*)pts;
  float _x0,_y0,x0,y0,x1,y1,x2,y2;
  float zcross, curz;
  int i;

  if ( npts <= 2 ) return 0; // not a polygon

  // check order (cw or ccw)
  _x0= x0 = *(float*)(curptr+0); _y0 = y0 = *(float*)(curptr+sizeof(float)); curptr+=stride;
  x1 = *(float*)(curptr+0); y1 = *(float*)(curptr+sizeof(float)); curptr+=stride;
  x2 = *(float*)(curptr+0); y2 = *(float*)(curptr+sizeof(float)); curptr+=stride;
  zcross = cpoly_zcross(x0,y0,x1,y1,x2,y2);

  // we have already the first three points
  curz = cpoly_zcross(x0,y0,x1,y1,x,y); if ( (curz*zcross) < 0.0f ) return 0;
  curz = cpoly_zcross(x1,y1,x2,y2,x,y); if ( (curz*zcross) < 0.0f ) return 0;

  // iterate all remaining points
  x0=x2; y0=y2;
  for ( i=3; i<npts;++i,curptr+=stride)
  {
    x1 = *(float*)(curptr+0); y1 = *(float*)(curptr+sizeof(float));
    curz = cpoly_zcross(x0,y0,x1,y1,x,y); if ( (curz*zcross) < 0.0f ) return 0;
    x0=x1; y0=y1;
  }

  // closes loop with the 1st
  curz = cpoly_zcross(x0,y0,_x0,_y0,x,y); if ( (curz*zcross) < 0.0f ) return 0;
  return 1;
}

// Separating Axis Theorem
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
        // The MTV is the minimum magnitude vector used to push the shapes out of the collision.
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


#endif
