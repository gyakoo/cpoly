# cdec2d
Some tests of convex Decomposition of 2D Polygons in C<br/>
It generates convex polygons out of any arbitrary polygon<br/>
Not based in any particular algorithm, just fooling with ideas.<br/>

OpenGL viewer provided. Depends on GLFW which is included but also in https://github.com/glfw/glfw

# How it works
the idea is the following, not implemented yet
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

