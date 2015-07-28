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

#define SCREENWIDTH 1024
#define SCREENHEIGHT 768
#define SCREENINVR (1.f/((float)SCREENHEIGHT/SCREENWIDTH))

// ortho proj
const double viewbounds[]={-100, 100, -80, 80, -1, 1};
const float WHITE[4]={1,1,1,1};

typedef struct sCircle
{
  float x, y, r;
}sCircle;

#define MAXCIRCLES 10
sCircle g_circles[MAXCIRCLES]={0};

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  int i;
  float* fc;
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

  
  for (i=0;i<MAXCIRCLES;++i)
  {
    drawCircle(g_circles[i].x, g_circles[i].y, g_circles[i].r, WHITE, 0, 0, 16.0f);
  }

  if ( glfwGetKey(window,GLFW_KEY_1)==GLFW_PRESS )
  {
    glPointSize(2.0f);
    cpoly_marching_sq(g_circles,MAXCIRCLES,sizeof(sCircle),5.0f);
    glBegin(GL_POINTS);
    for ( i=0;i<cpoly_pool_vcount;++i )
    {
      cpoly_pool_get_vertex(i,&x,&y);
      glVertex2f(x,y);
    }
    glEnd();
  }

	
  glPointSize(15.0f);
  glColor4ub(255,255,0,255);
  glBegin(GL_POINTS);
  glVertex2f(0,0);
  glEnd();

  glfwSwapBuffers(window);
}

void resizecb(GLFWwindow* window, int width, int height)
{
	frame(window);
}


float randrange(float a, float b)
{
  return a + ((float)rand()/RAND_MAX*fabs(b-a));
}

void randomcircles()
{
  int i;

  // init circles
  for (i=0;i<MAXCIRCLES;++i)
  {
    g_circles[i].x = randrange(-40.0f,40.0f);
    g_circles[i].y = randrange(-40.0f,40.0f);
    g_circles[i].r = randrange(2.0f, 12.0f);
  }
}


void keycallback(GLFWwindow* w, int key, int scancode, int action, int mods)
{
  if ( key==GLFW_KEY_ESCAPE && action == GLFW_PRESS )
    glfwSetWindowShouldClose(w, GL_TRUE);
  if ( key==GLFW_KEY_SPACE && action==GLFW_PRESS)
    randomcircles();
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

  randomcircles();

	while (!glfwWindowShouldClose(window))
	{
		frame(window);


		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}
