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
#define SQSIDE 2.0f

// ortho proj
const double viewbounds[]={-100, 100, -80, 80, -1, 1};
const float WHITE[4]={1,1,1,1};

typedef struct sCircle
{
  float x, y, r;
}sCircle;

#define MAXCIRCLES 10
sCircle g_circles[MAXCIRCLES]={0};

void drawgrid()
{
  int i,j,stepsx,stepsy;
  float x,y;
  float minis[2]={FLT_MAX,FLT_MAX};
  float maxis[2]={-FLT_MAX, -FLT_MAX};
  sCircle* cir;

  // computes AABB of circles
  for (i=0;i<MAXCIRCLES;++i)
  {
    cir = g_circles+i;
    x=cir->x-cir->r; y=cir->y-cir->r; if ( x < minis[0] ) minis[0]=x; if ( y < minis[1] ) minis[1]=y;
    x=cir->x+cir->r; y=cir->y+cir->r; if ( x > maxis[0] ) maxis[0]=x; if ( y > maxis[1] ) maxis[1]=y;
  }
  minis[0]-=SQSIDE; maxis[0]+=SQSIDE;
  maxis[1]+=SQSIDE; minis[1]-=SQSIDE;

  stepsx = (int)ceilf((maxis[0]-minis[0])/SQSIDE);
  stepsy = (int)ceilf((maxis[1]-minis[1])/SQSIDE);

  maxis[0] = minis[0]+SQSIDE*stepsx;
  minis[1] = maxis[1]-SQSIDE*stepsy;

  glLineWidth(1.0f);
  glColor4f(0.5f,0.5f,0.5f,0.6f);
  glBegin(GL_LINES);  
  y=maxis[1];
  for ( j=0;j<=stepsy;++j,y-=SQSIDE)
  {
    glVertex2f( minis[0], y );
    glVertex2f( maxis[0], y );
  }
  x=minis[0];
  for (i=0;i<=stepsx;++i,x+=SQSIDE)
  {
    glVertex2f( x, maxis[1] );
    glVertex2f( x, minis[1] );
  }
  glEnd();
}

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

  
  if ( glfwGetKey(window,GLFW_KEY_1)==GLFW_PRESS )
  {
    glPointSize(4.0f);
    cpoly_marchingsq_nointerp(g_circles,MAXCIRCLES,sizeof(sCircle),SQSIDE);    
    glLineWidth(1.0f);
    glColor4f(1,1,0,1);

    k=0;
    for (i=0;i<cpoly_pool_icount[CPOLY_IPOOL_0];++i)
    {
      glBegin(GL_LINE_LOOP);
      for (j=k;j<cpoly_pool_get_index(CPOLY_IPOOL_0,i);++j )
      {
        cpoly_pool_get_vertex(j,&x,&y);
        glVertex2f(x,y);
      }
      glEnd();
      k=j;
    }

//     cpoly_convex_hull(cpoly_pool_v,cpoly_pool_get_index(CPOLY_IPOOL_0,0),sizeof(float)*2);
//     glColor4ub(0,255,0,255);
//     glBegin(GL_LINE_LOOP);
//     for ( i=0; i < cpoly_pool_icount[CPOLY_IPOOL_0]; ++i )
//     {
//       j = cpoly_pool_get_index(CPOLY_IPOOL_0,i);
//       cpoly_pool_get_vertex(j,&x,&y);
//       glVertex2f( x,y );
//     }
//     glEnd();
  }
  //else
  {
    for (i=0;i<MAXCIRCLES;++i)
    {
      drawCircle(g_circles[i].x, g_circles[i].y, g_circles[i].r, WHITE, 0, 0, 16.0f);
    }
  }
  drawgrid();
	
//   glPointSize(15.0f);
//   glColor4ub(255,255,0,255);
//   glBegin(GL_POINTS);
//   glVertex2f(0,0);
//   glEnd();

  glfwSwapBuffers(window);
}

void resizecb(GLFWwindow* window, int width, int height)
{
	frame(window);
}


float randrange(float a, float b)
{
  return a + (float)((float)rand()/RAND_MAX*fabs(b-a));
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
