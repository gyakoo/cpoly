# cdec2d
Some tests of convex Decomposition of 2D Polygons in C<br/>
It generates convex polygons out of any arbitrary polygon<br/>
Not based in any particular algorithm, just fooling with ideas.<br/>

OpenGL viewer provided. Depends on GLFW which is included but also in https://github.com/glfw/glfw

# How it works (work in progress)
the idea is the following:
<pre>
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
</pre>

# Usage (work in progress)
```c
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
```
