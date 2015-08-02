//
//
#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <GLFW/glfw3.h>
#pragma comment(lib, "glfw3.lib")

#define CPOLY_IMPLEMENTATION
#include "cpoly.h"

#define EXCOMMON_IMPLEMENTATION
#include "excommon.h"

#define SCREENWIDTH 800
#define SCREENHEIGHT 600
#define SCREENINVR (1.f/((float)SCREENHEIGHT/SCREENWIDTH))
#define MINRND 5
#define RANGRND 30.0f
#define POLYPOINTS 12

// original polygon to be decomposed
float g_polygon[POLYPOINTS*2]={0};
const int g_polycount= sizeof(g_polygon)/(sizeof(float)*2);
float g_polygon2[12*2]={0};
float g_poly2pos[2]={0};
const int g_polycount2= sizeof(g_polygon2)/(sizeof(float)*2);
int fillpoly=0;
const float STRIDE=sizeof(float)*2;
// ortho proj
const double viewbounds[]={-100, 100, -80, 80, -1, 1};

void randompoly()
{
  int i,j;
  float step;
  float a;
  float x,y;
  float l;
  int count[2]={g_polycount,g_polycount2};
  float* verts[2]={g_polygon,g_polygon2};

  for (j=0;j<2;++j)
  {
    step=3.14159f*2.0f/count[j];
    a=3.14159f*2.0f;
    for ( i=0; i < count[j];++i, a-=step)
    {
      l = MINRND + RANGRND*((float)rand()/RAND_MAX);
      x=cosf(a)*l;
      y=sinf(a)*l;
      verts[j][i*2]=x;
      verts[j][i*2+1]=y;
    }
  }
  g_poly2pos[0]=g_poly2pos[1]=0.0f;
}

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  int i,j,k,l=0;
  float x,y;
  static float globalTime=0.0f;
//  float x,y;

  globalTime += 1.0f/60.0f;
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	glClearColor(20.0f/255.0f, 20.0f/255.0f, 90.0f/255.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_TEXTURE_2D);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(viewbounds[0], viewbounds[1], viewbounds[2], viewbounds[3], viewbounds[4], viewbounds[5]);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDisable(GL_DEPTH_TEST);
	glColor4ub(255,255,255,255);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

  glLineWidth(2.0f);
  glBegin(GL_LINE_LOOP);
  for ( i=0;i<g_polycount;++i ) glVertex2f(g_polygon[i*2], g_polygon[i*2+1]);
  glEnd();

  glPointSize(6.0f);
  glColor4f(0,1,0,1);
  glBegin(GL_POINTS);
  for ( i=0;i<g_polycount;++i ) glVertex2f(g_polygon[i*2], g_polygon[i*2+1]);
  glEnd();

  glPointSize(8.0f);
  glColor4f(1,0,0,1);
  glBegin(GL_POINTS);
  glVertex2f(g_polygon[0], g_polygon[1]);
  glEnd();

  if ( glfwGetKey(window,GLFW_KEY_2)==GLFW_PRESS )
  {
    cpoly_triangulate_EC(g_polygon,g_polycount,STRIDE,NULL);
    if ( fillpoly )
    {
      glBegin(GL_TRIANGLES);
      for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_0]/3;++i) // for all partitions
      {
        j=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3);
        k=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3+1);
        l=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3+2);
        glVertex2f( g_polygon[j*2], g_polygon[j*2+1] );
        glVertex2f( g_polygon[k*2], g_polygon[k*2+1] );
        glVertex2f( g_polygon[l*2], g_polygon[l*2+1] );
      }
      glEnd();
    }
    else
    {
      for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_0]/3;++i) // for all partitions
      {
        glBegin(GL_LINE_LOOP);
        j=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3);
        k=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3+1);
        l=cpoly_pool_get_index(CPOLY_IPOOL_0,i*3+2);
        glVertex2f( g_polygon[j*2], g_polygon[j*2+1] );
        glVertex2f( g_polygon[k*2], g_polygon[k*2+1] );
        glVertex2f( g_polygon[l*2], g_polygon[l*2+1] );
        glEnd();
      }
    }
  }
  else if ( glfwGetKey(window,GLFW_KEY_3)==GLFW_PRESS )
  {
    cpoly_partition_HM(g_polygon, g_polycount, STRIDE);
    k=0;
    for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_1];++i)
    {
      glBegin(GL_LINE_LOOP);
      for (j=k;j<cpoly_pool_get_index(CPOLY_IPOOL_1,i);++j)
      {
        l=cpoly_pool_get_index(CPOLY_IPOOL_0,j);
        glVertex2f( g_polygon[l*2], g_polygon[l*2+1] );
      }
      k=j;
      glEnd();
    }
  }
  else if ( glfwGetKey(window,GLFW_KEY_1)==GLFW_PRESS )
  {
    cpolyBitPool pvs;
    cpoly_pvs_create(&pvs,g_polycount);

    cpoly_pvs(g_polygon,g_polycount,STRIDE,&pvs);
    glBegin(GL_LINES);
    for (i=0;i<g_polycount;++i)
    {
      for(j=i+2;j<g_polycount;++j)
      {
        if ( cpoly_pvs_get(&pvs,i,j) )
        {
          glVertex2f(g_polygon[i*2],g_polygon[i*2+1]);
          glVertex2f(g_polygon[j*2],g_polygon[j*2+1]);
        }
      }
    }
    glEnd();
    cpoly_pvs_destroy(&pvs);
  }

  if ( glfwGetKey(window, GLFW_KEY_A)==GLFW_PRESS )
    cpoly_transform_rotate(g_polygon2, g_polycount2, STRIDE, 0.016f, g_poly2pos, g_poly2pos+1);
  else if ( glfwGetKey(window, GLFW_KEY_Z)==GLFW_PRESS )
    cpoly_transform_rotate(g_polygon2, g_polycount2, STRIDE, -0.016f,NULL,NULL);
  else if ( glfwGetKey(window, GLFW_KEY_LEFT)==GLFW_PRESS )
  {
    g_poly2pos[0]-=0.016f*20.0f;
    cpoly_transform_translate(g_polygon2, g_polycount2, STRIDE, g_poly2pos[0], g_poly2pos[1], NULL, NULL);
  }
  else if ( glfwGetKey(window, GLFW_KEY_RIGHT)==GLFW_PRESS )
  {
    g_poly2pos[0]+=0.016f*20.0f;
    cpoly_transform_translate(g_polygon2, g_polycount2, STRIDE, g_poly2pos[0], g_poly2pos[1], NULL, NULL);
  }

  if ( glfwGetKey(window, GLFW_KEY_Q)==GLFW_PRESS )
  {
    glBegin(GL_LINE_LOOP);
    for ( i=0;i<g_polycount2;++i ) glVertex2f(g_polygon2[i*2], g_polygon2[i*2+1]);
    glEnd();
//     cpoly_cv_clip_union(g_polygon, g_polycount, STRIDE, g_polygon2, g_polycount2, STRIDE);
//     glBegin(GL_LINE_LOOP);
//     for ( i = 0; i < cpoly_pool_vcount; ++i )
//     {
//       cpoly_pool_get_vertex(i,&x,&y);
//       glVertex2f(x,y);
//     }
//     glEnd();
  }

  	
  glfwSwapBuffers(window);
}

void resizecb(GLFWwindow* window, int width, int height)
{
	frame(window);
}

void keycallback(GLFWwindow* w, int key, int scancode, int action, int mods)
{
  if ( key==GLFW_KEY_ESCAPE && action == GLFW_PRESS )
    glfwSetWindowShouldClose(w, GL_TRUE);
  if ( key==GLFW_KEY_SPACE && action == GLFW_PRESS )
    randompoly();
  if (key==GLFW_KEY_TAB && action==GLFW_PRESS)
    fillpoly=1-fillpoly;
}


int main()
{
  int i;
	GLFWwindow* window;
	const GLFWvidmode* mode;
  double t0,t1,acum=0.0;
  char timestr[32];
  char* algo;

#ifdef _DEBUG
  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

	if (!glfwInit())
		return -1;

	mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
  window = glfwCreateWindow(SCREENWIDTH, SCREENHEIGHT, "cpoly", NULL, NULL);
	if (!window)
	{
		printf("Could not open window\n");
		glfwTerminate();
		return -1;
	}

  randompoly();
  for (i=0;i<6;++i) randompoly();
	glfwSetFramebufferSizeCallback(window, resizecb);
	glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, keycallback);
  glfwSwapInterval(1);

  glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);

	while (!glfwWindowShouldClose(window))
	{
    t0=glfwGetTime();
		frame(window);
		glfwPollEvents();
    t1=(glfwGetTime() - t0);
    acum += t1;
    if ( acum >= 1.0 )
    {
      sprintf(timestr,"%.2g", 1.0/t1);
      glfwSetWindowTitle(window,timestr);
      acum-=1.0;
    }
	}

	glfwTerminate();
	return 0;
}
