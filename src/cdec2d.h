/*

Usage:

int N = 20;
float polygonPoints[N]={...};
int** parts=0;
int* psizes=0;
int count=0;
int i,j,offs;

count = cdec2d_decomp_cw( polygonPoints, N, sizeof(float)*2, &parts, &psizes);
if ( count > 0 )
{
  for (i=0;i<count;++i)
  {
    printf( "Partition %d points:\n", i );
    for (j=0;j<psizes[i];++j)
    {
      offs = parts[i][j]*2;
      printf( "\tp%03d: %f, %f\n", polygonPoints[offs], polygonPoints[offs+1] );
    }
    printf ( "\n" );
  }
  cdec2d_free_parts(&parts, &psizes, count);
}

 */

#ifndef CDEC2D_H
#define CDEC2D_H

#ifdef __cplusplus
extern "C" {
#endif

  // Decomposes a polygon into a set of convex polygons (ClockWise)
  // 'pts' is the array of polygon points, it should point to a pair of float x,y.
  // 'npts' is the number of points
  // 'stride' is the no. of bytes until next point
  // 'parts' will return an array of arrays, indicating the indices of points.
  // 'psizes' size for every array in 'parts'
  // return the number of parts
  int cdec2d_decomp_cw(void* pts, int npts, int stride, int*** parts, int** psizes);

  // Deallocates memory previously allocated by cdec2d_decomp* functions
  void cdec2d_free_parts(int*** parts, int** psizes, int count);

#ifdef __cplusplus
};
#endif

#endif // CDEC2D_H

#ifdef CDEC2D_IMPLEMENTATION


#ifdef _MSC_VER
	#pragma warning (disable: 4996) // Switch off security warnings
	#pragma warning (disable: 4100) // Switch off unreferenced formal parameter warnings
#endif


int cdec2d_decomp_cw(void* pts, int npts, int stride, int*** parts, int** psizes)
{
  // sample test 
  int nparts=3;

  *parts = (int**)malloc(sizeof(int)*nparts);
  *psizes= (int*)malloc(sizeof(int)*nparts);

  (*parts)[0] = (int*)malloc(sizeof(int)*5); (*psizes)[0]=5;
  (*parts)[0][0] = 0;
  (*parts)[0][1] = 1;
  (*parts)[0][2] = 4;
  (*parts)[0][3] = 5;
  (*parts)[0][4] = 6;
  (*parts)[1] = (int*)malloc(sizeof(int)*4); (*psizes)[1]=4;
  (*parts)[1][0] = 1;
  (*parts)[1][1] = 2;
  (*parts)[1][2] = 3;
  (*parts)[1][3] = 4;
  (*parts)[2] = (int*)malloc(sizeof(int)*3); (*psizes)[2]=3;
  (*parts)[2][0] = 6;
  (*parts)[2][1] = 7;
  (*parts)[2][2] = 0;

  return nparts;
  /*
   LOOP1
   if first partition, 
    select starting point and next one (sp=current, np=next)
   else
    select prior and next one (sp=prior, np=current)
   P += {sp,np}
   
     LOOP2
     compute plane of last edge (sp-->np)
     p = only consider next points positives to plane 
     if no intersect (sp-->p) with any edge, p
          P+={p}
     last edge=np->p
     Goto LOOP2 while points
   
   add partition P
   add adjacent points to remaining points
   Goto LOOP1 while points (if 1+2 remaining points is last partition though)

   */
}

void cdec2d_free_parts(int*** parts, int** psizes, int count)
{
  int i=0; 

  if ( !parts || !*parts || !psizes || !*psizes )
    return;
  for ( ; i < count; ++i )
    free((*parts)[i]);

  free(*parts);
  free(*psizes);
  *parts=0;
  *psizes=0;
}


#endif
