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

// ortho proj
const double viewbounds[]={-100, 100, -100, 100, -1, 1};

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

float g_convexpoly1[] ={
  5.0f,25.0f,  // 0
  15.0f, 30.0f,   // 1
  25.0f, 25.0f,  // 2
  35.0f, 0.0f,  // 3
  40.0f, -20.0f,    // 4
  30.0f, -35.0f,    // 5
  10.0f, -35.0f,  // 6
  -5.0f, -10.0f // 7
};
const int g_convexpolycount1= sizeof(g_convexpoly1)/(sizeof(float)*2);
float g_convexTransform[sizeof(g_convexpoly1)]={0};
const float C_REDISH[4]={50/255.0f,0,10/255.0f,1};
// decomp result
int* g_parts=0;
int* g_psizes=0;
int g_partscount=0;


void drawPolygon(char fill, float* poly, int count, float x, float y, float* color)
{
  int i;

  glLineWidth(2.0f);
  glPointSize(7.0f);

  glPushMatrix();
  glTranslatef(x,y,0.0f);
  glColor4fv(color);
  if ( fill )
  {
    glBegin(GL_TRIANGLE_FAN);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();
  }
  else
  {
    glColor4ub(255,255,255,255);
    glBegin(GL_LINE_LOOP);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();

    glColor4ub(255,255,0,255);
    glBegin(GL_POINTS);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();
  }
  glPopMatrix();
}

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  float x0,y0,x1,y1,d0,d1;
  float a,step;
  int i;
  float* fc;
  static float globalTime=0.0f;
  float othc[4]={120/255.0f,0,10/255.0f,1};

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

	// Draw 
  //drawPolygon(1, g_convexpoly0, g_convexpolycount0, 0.0f, 0.0f,C_REDISH);
  drawPolygon(1, g_convexpoly1, g_convexpolycount1, 0.0f, 0.0f,C_REDISH);	
  //drawPolygon(0, g_convexpoly0, g_convexpolycount0, 0.0f, 0.0f,C_REDISH);
  drawPolygon(0, g_convexpoly1, g_convexpolycount1, 0.0f, 0.0f,C_REDISH);


  glColor4ub(255,0,255,255);    
  glBegin(GL_POINTS);
  for ( i = 0; i < g_convexpolycount0; ++i )
  {
    x0 = g_convexpoly0[i*2]; y0= g_convexpoly0[i*2+1];
    if ( cpoly_cv_point_inside(g_convexpoly1,g_convexpolycount1,sizeof(float)*2,x0,y0) )
      glVertex2f(x0,y0);
  }
  glEnd();

  for ( i=0; i< g_convexpolycount1; ++i )
  {
    g_convexTransform[i*2] = g_convexpoly1[i*2]+cos(globalTime*0.5f)*60.0f;
    g_convexTransform[i*2+1] = g_convexpoly1[i*2+1];
  }
  fc = cpoly_cv_intersects_SAT(g_convexpoly1, g_convexpolycount1, sizeof(float)*2, g_convexTransform, g_convexpolycount1, -1) ? othc : C_REDISH;
  drawPolygon(1, g_convexTransform, g_convexpolycount1, 0.0f, 0.0f,fc);
  drawPolygon(0, g_convexTransform, g_convexpolycount1, 0.0f, 0.0f,fc);
  

  /*
  x0=cos(ta)*20.0f,y0=sin(ta)*20.0f;
  glBegin(GL_LINES);
    glColor4ub(255,255,255,255);
    glVertex2f(0,0); glVertex2f(x0,y0);
    glColor4ub(255,255,255,100);
    glVertex2f(0,0); glVertex2f(-x0,-y0);
  glEnd();

  glPointSize(4.0f);
  glBegin(GL_POINTS);
  a=.0f; step=3.14159f*2.0f/32;
  for (i=0;i<32;++i,a+=step)
  {
    x1=cos(a)*15.0f; y1=sin(a)*15.0f;
    d0 = cpoly_zcross(0,0,x0,y0,x1,y1);
    if ( d0<=.0f) glColor4ub(255,255,255,255); else glColor4ub(0,255,0,255);
    glVertex2f(x1,y1);
  }
  glEnd();
  */
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

  if ( !cpoly_is_convex(g_convexpoly0,g_convexpolycount0,sizeof(float)*2) ||
       !cpoly_is_convex(g_convexpoly1,g_convexpolycount1,sizeof(float)*2) )
  {
    printf( "Polygons should be convex for this example\n" );
    exit(EXIT_FAILURE);
  }

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
  glfwSwapInterval(1);

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
