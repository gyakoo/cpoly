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

// ortho proj
const double viewbounds[]={-80, 80, -80, 80, -1, 1};

float STRIDE=sizeof(float)*2;
float g_convexpoly0[] ={
  -10.0f  , 35.0f,  // 0
  0.0f    , 40.0f,   // 1
  10.0f   , 35.0f,  // 2
  20.0f   , 10.0f,  // 3
  25.0f   , -10.0f, // 4
  -5.0f   , -25.0f, // 5
  -10.0f  , 0.0f // 6
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
float g_pos1[2];
float g_poly[] = 
{
  -10.0f  , 35.0f,  // 0
  0.0f    , 40.0f,   // 1
  10.0f   , 5.0f,  // 2
  20.0f   , 10.0f,  // 3
  15.0f   , -10.0f, // 4
  0.0f   , 0.0f, // 5
  -20.0f   , -25.0f, // 6
  -5.0f  , 10.0f // 7
};
const int g_polycount= sizeof(g_poly)/(sizeof(float)*2);

float g_convexTransform[sizeof(g_convexpoly1)]={0};
float C_REDISH[4]={50/255.0f,0,10/255.0f,0.5f};
float C_GREENISH[4]={0,50/255.0f,10/255.0f,1};
float C_GREEN[4]={0,1,0,1};

void frame(GLFWwindow* window)
{
	int width = 0, height = 0;
  float x0,y0,x1,y1,x2,y2,x3,y3;
  float a,step;
  int i,j;
  float* fc;
  static float globalTime=0.0f;
  float othc[4]={1,0,10/255.0f,0.3f};
  float dslow=globalTime*0.01f;

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
  if ( glfwGetKey(window,GLFW_KEY_SPACE)==GLFW_PRESS )
  {
    cpoly_cv_union(g_convexpoly0, g_convexpolycount0, STRIDE, g_convexpoly1, g_convexpolycount1, STRIDE);
    glBegin(GL_LINE_LOOP);
    for ( i = 0; i < cpoly_pool_count; ++i )
    {
      cpoly_pool_get_vertex(i,&x0,&y0);
      glVertex2f(x0,y0);
    }
    glEnd();
  }
  else if ( glfwGetKey(window,GLFW_KEY_BACKSPACE)==GLFW_PRESS )
  {
    glBegin(GL_LINE_LOOP);
    for ( i = 0; i < g_polycount; ++i )
      glVertex2f( g_poly[i*2], g_poly[i*2+1]);
    glEnd();

    cpoly_convex_hull(g_poly, g_polycount,STRIDE);
    glColor4ub(0,255,0,255);
    glBegin(GL_LINE_LOOP);
    for ( i=0; i < cpoly_pool_count; ++i )
    {
      j = cpoly_pool_get_index(i);
      glVertex2f( g_poly[j*2], g_poly[j*2+1] );
    }
    glEnd();

  }
  else if ( glfwGetKey(window,GLFW_KEY_1)==GLFW_PRESS )
  {
    cpoly_cv_diff(g_convexpoly0,g_convexpolycount0, STRIDE, g_convexpoly1, g_convexpolycount1, STRIDE);
    glBegin(GL_LINE_LOOP);
    for ( i = 0; i < cpoly_pool_count; ++i )
    {
      cpoly_pool_get_vertex(i,&x0,&y0);
      glVertex2f( x0, y0);
    }
    glEnd();
  }
  else
  {
    drawPolygon(1, g_convexpoly0, g_convexpolycount0, 0.0f, 0.0f, C_GREENISH);
    drawPolygon(1, g_convexpoly1, g_convexpolycount1, 0.0f, 0.0f, C_REDISH);	
    drawPolygon(0, g_convexpoly0, g_convexpolycount0, 0.0f, 0.0f, C_GREENISH);
    drawPolygon(0, g_convexpoly1, g_convexpolycount1, 0.0f, 0.0f, C_REDISH);

    glColor4ub(255,0,255,255);    
    glBegin(GL_POINTS);
    for ( i = 0; i < g_convexpolycount0; ++i )
    {
      x0 = g_convexpoly0[i*2]; y0= g_convexpoly0[i*2+1];
      if ( cpoly_cv_point_inside(g_convexpoly1,g_convexpolycount1,STRIDE,x0,y0) )
        glVertex2f(x0,y0);
    }
    glEnd();

    // cpu transform of polygon g_convexTransform
//     for ( i=0; i< g_convexpolycount1; ++i )
//     {
//       g_convexTransform[i*2] = g_convexpoly1[i*2]+cos(globalTime*0.5f)*60.0f;
//       g_convexTransform[i*2+1] = g_convexpoly1[i*2+1]+sin(globalTime*2.0f)*10.0f;
//     }
// 
//     // coloring if both collide
//     fc = cpoly_cv_intersects_SAT(g_convexpoly1, g_convexpolycount1, STRIDE, g_convexTransform, g_convexpolycount1, -1) ? othc : C_REDISH;
//     drawPolygon(1, g_convexTransform, g_convexpolycount1, 0.0f, 0.0f,fc);
//     drawPolygon(0, g_convexTransform, g_convexpolycount1, 0.0f, 0.0f,fc);
  }

  if ( glfwGetKey(window, GLFW_KEY_A)==GLFW_PRESS )
  {
    //cpoly_transform_translate(g_convexpoly1, g_convexpolycount1, STRIDE, 0, 0);
    cpoly_transform_rotate(g_convexpoly1, g_convexpolycount1, STRIDE, 0.016f,NULL,NULL);
  }

  if ( glfwGetKey(window, GLFW_KEY_LEFT)==GLFW_PRESS )
  {
    g_pos1[0] -= 1.16f;
    cpoly_transform_translate(g_convexpoly1, g_convexpolycount1, STRIDE, g_pos1[0], g_pos1[1]);
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

void keycallback(GLFWwindow* w, int key, int scancode, int action, int mods)
{
  if ( key==GLFW_KEY_ESCAPE && action == GLFW_PRESS )
    glfwSetWindowShouldClose(w, GL_TRUE);
}

int main()
{
	GLFWwindow* window;
	const GLFWvidmode* mode;

  int v;

  if ( !cpoly_is_convex(g_convexpoly0,g_convexpolycount0,STRIDE) ||
       !cpoly_is_convex(g_convexpoly1,g_convexpolycount1,STRIDE) )
  {
    printf( "Polygons should be convex for this example\n" );
    exit(EXIT_FAILURE);
  }

  //cpoly_cv_diff(g_convexpoly0, g_convexpolycount0, STRIDE, g_convexpoly1, g_convexpolycount1, STRIDE);
  cpoly_poly_centroid(g_convexpoly1,g_convexpolycount1,STRIDE,g_pos1,g_pos1+1);

  cpoly_transform_scale(g_poly, g_polycount, STRIDE, 1.5f, 1.5f,0,0);

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
