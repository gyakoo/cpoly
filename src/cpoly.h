/*

Usage:

int N = 20;
float polygonPoints[N]={...};
int* partIndices=0;
int* partOffsets=0;
int count=0;
int i,j,start,end;
float* offs;

count = cpoly_partitioning( polygonPoints, N, sizeof(float)*2, &partIndices, &partOffsets);
if ( count > 0 )
{
  start = 0;
  for (i=0;i<count;++i)
  {
    printf( "Partition %d points:\n", i );
    end = partOffsets[i];
    for ( j=start; j<end; ++j )
    {
      offs = polygonPoints + partIndices[j];
      printf( "\tPoint %d: %f %f\n", j-start, *offs, *(offs+1) );
    }
    start = end;
    printf ( "\n" );
  }
  cpoly_free_parts(&partIndices, &partOffsets, count);
}

 */

#ifndef CPOLY_H
#define CPOLY_H

#ifdef __cplusplus
extern "C" {
#endif

  // Decomposes a polygon into a set of convex polygons (ClockWise)
  // 'pts' is the array of polygon points, it should point to a pair of float x,y.
  // 'npts' is the number of points
  // 'stride' is the no. of bytes until next point
  // 'partndx' will allocate an array with all partition indices. see remarks
  // 'poffsets' will allocate an array with offsets in 'partndxs' with partitions. length = return count - 1. see remarks
  // return the number of parts
  // Remarks:
  // 
  int cpoly_partitioning_cw(void* pts, int npts, int stride, int** partndxs, int** poffsets);

  // Deallocates memory previously allocated by cpoly_decomp* functions
  void cpoly_free_parts(int** parts, int** psizes);

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
  int i=0; 

  if ( !partndxs || !*partndxs || !poffsets || !*poffsets )
    return;
  free(*partndxs);
  free(*poffsets);
  *partndxs = *poffsets = 0;
}


#endif
