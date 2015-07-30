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
float g_polygon[]={ 
  -10.0f,20.0f,  // 0
  0.0f, 10.0f,   // 1
  10.0f, 20.0f,  // 2
  20.0f, -5.0f,  // 3
  5.0f, 0.0f,    // 4
  -5.0f, -40.0f,  // 5
  -5.0f, -10.0f, // 6
  -20.0f, -15.0f // 7
};
const int g_polycount= sizeof(g_polygon)/(sizeof(float)*2);

// ortho proj
const double viewbounds[]={-100, 100, -80, 80, -1, 1};

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  int i,j,k;
  static float globalTime=0.0f;
  float x,y;

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

  if ( glfwGetKey(window,GLFW_KEY_SPACE)==GLFW_PRESS )
  {
    cpoly_convex_partition(g_polygon,g_polycount,sizeof(float)*2);
    glBegin(GL_LINE_LOOP);
    for (i=0;i<cpoly_pool_icount;++i)
    {
      k = cpoly_pool_get_index(i);
      glVertex2f(g_polygon[k*2],g_polygon[k*2+1]);
    }
    glEnd();
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
}


int main()
{
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
