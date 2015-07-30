//
//

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

#define SCREENWIDTH 1600
#define SCREENHEIGHT 1200
#define SCREENINVR (1.f/((float)SCREENHEIGHT/SCREENWIDTH))

// original polygon to be decomposed
float g_polygon[12*2]={0};
const int g_polycount= sizeof(g_polygon)/(sizeof(float)*2);
int fillpoly=0;

// ortho proj
const double viewbounds[]={-100, 100, -80, 80, -1, 1};

void randompoly()
{
  int i;
  float step=3.14159f*2.0f/g_polycount;
  float a=0.0f;
  float x,y;
  float l;

  for ( i=0; i < g_polycount;++i, a-=step)
  {
    l = 2.0f + 30.0f*((float)rand()/RAND_MAX);
    x=cos(a)*l;
    y=sin(a)*l;
    g_polygon[i*2]=x;
    g_polygon[i*2+1]=y;
  }
}

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  int i,j,k;
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
    cpoly_convex_partition(g_polygon,g_polycount,sizeof(float)*2);

    k=0;
    for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_1];++i) // for all partitions
    {
      glBegin(GL_LINE_LOOP);     
      for ( j=k; j<cpoly_pool_get_index(CPOLY_IPOOL_1,i);++j)
      {
        k=cpoly_pool_get_index(CPOLY_IPOOL_0,j);
        glVertex2f(g_polygon[k*2],g_polygon[k*2+1]);
      }
      k=j;
      glEnd();
    }

    if ( fillpoly )
    {
      k=0;
      for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_1];++i) // for all partitions
      {
        glBegin(GL_TRIANGLE_FAN);
        for ( j=k; j<cpoly_pool_get_index(CPOLY_IPOOL_1,i);++j)
        {
          k=cpoly_pool_get_index(CPOLY_IPOOL_0,j);
          glVertex2f(g_polygon[k*2],g_polygon[k*2+1]);
        }
        k=j;
        glEnd();
      }
    }
  }
  else if ( glfwGetKey(window,GLFW_KEY_1)==GLFW_PRESS )
  {
    cpolyBitPool pvs;
    cpoly_pvs_create(&pvs,g_polycount);

    cpoly_pvs(g_polygon,g_polycount,sizeof(float)*2,&pvs);
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
		frame(window);

		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}
