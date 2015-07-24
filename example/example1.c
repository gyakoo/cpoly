//
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <GLFW/glfw3.h>
#pragma comment(lib, "glfw3.lib")

#define CPOLY_IMPLEMENTATION
#include "cpoly.h"

#define EXCOMMON_IMPLEMENTATION
#include "excommon.h"
// ortho proj
const double viewbounds[]={-100, 100, -100, 100, -1, 1};

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

float g_convexpoly0[] ={
  -10.0f,20.0f,  // 0
  0.0f, 25.0f,   // 1
  10.0f, 20.0f,  // 2
  20.0f, -5.0f,  // 3
  25.0f, -25.0f,    // 4
  -5.0f, -40.0f,  // 5
  -20.0f, -15.0f // 6
};
const int g_convexpolycount0= sizeof(g_convexpoly0)/(sizeof(float)*2);

// decomp result
int* g_parts=0;
int* g_psizes=0;
int g_partscount=0;

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  float rgba[]={1,1,1,1};

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

	// Draw bounds
  drawPolygon(1, g_polygon, g_polycount, -40.0f, 40.0f, rgba);
  drawPolygon(0, g_polygon, g_polycount, +40.0f, 40.0f, rgba);	

  drawPolygon(1, g_convexpoly0, g_convexpolycount0, -40.0f, -40.0f, rgba);
  drawPolygon(0, g_convexpoly0, g_convexpolycount0, +40.0f, -40.0f, rgba);	

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

  cpoly_is_convex(g_convexpoly0,g_convexpolycount0,sizeof(float)*2);

	if (!glfwInit())
		return -1;

	mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
  window = glfwCreateWindow(1024, 768, "cpoly", NULL, NULL);
	if (!window)
	{
		printf("Could not open window\n");
		glfwTerminate();
		return -1;
	}

	glfwSetFramebufferSizeCallback(window, resizecb);
	glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, keycallback);

  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
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
