/*

 */

#ifndef CPOLY_H
#define CPOLY_H

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

  // Return 1 if two convex polygons intersects
  int cpoly_cv_intersects(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1);

#ifdef __cplusplus
};
#endif

#endif // CPOLY_H

#ifdef CPOLY_IMPLEMENTATION


#ifdef _MSC_VER
	#pragma warning (disable: 4996) // Switch off security warnings
	#pragma warning (disable: 4100) // Switch off unreferenced formal parameter warnings
#endif


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

// NEGATIVE is to the eye (CW) POSITIVE is towards the screen (CCW)
#define cpoly_zcross(x0,y0, x1, y1, x2, y2) ((x1-x0)*(y2-y1) - (y1-y0)*(x2-x1))

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

// sweeping to axes and finding
int cpoly_cv_intersects(void* pts0, int npts0, int stride0, void* pts1, int npts1, int stride1/*=-1*/)
{
  char* ptr[2]={(char*)pts0, (char*)pts1};
  const int pts[2]={npts0,npts1};
  const int strides[2]={stride0, stride1<=0?stride0:stride1};
  float minis[2][2]={FLT_MAX,FLT_MAX};
  float maxis[2][2]={-FLT_MAX,-FLT_MAX};
  float x,y;
  int i;
  char isect=0;
  register float a,b,c,d;

  // sweeping to axes to early discard far away polygons
  // compute maximum and minimum of both polys in x and y
//   for (i=0;i<2;++i)
//   {
//     x=*(float*)(ptr[i]); y=*(float*)(ptr[i]+sizeof(float));
//     if ( x<=minis[i][0] ) minis[i][0]=x; if ( x>maxis[i][0] ) maxis[i][0]=x;
//     if ( y<=minis[i][1] ) minis[i][1]=y; if ( y>maxis[i][1] ) maxis[i][1]=y;
//     ptr[i]+=strides[i];
//   }
// 
//   // check in X and Y-axis
//   for ( i=0; i<2; ++i )
//   {
//     a = minis[0][i]; b=maxis[0][i];
//     c = minis[1][i]; d=maxis[1][i];
//     if ( c>=a && c<=b || d>=a && d<=b ) 
//       ++isect;
//   }
//   if ( !isect ) 
//     return 0;
  
  // if isect==1 || isect==2, possible intersection. 
  // todo: place a good algorithm here, this is O(NM), probably good enough for small polygons
  // but for ones with huge no. of vertices isn't good.
  ptr[1]=(char*)pts1;
  for ( i=0;i<npts1;++i,ptr[1]+=strides[1])
  {
    x=*(float*)(ptr[1]); y=*(float*)(ptr[1]+sizeof(float));
    if ( cpoly_cv_point_inside(pts0,npts0,stride0,x,y) )
      return 1;
  }
  return 0;
}


#endif
