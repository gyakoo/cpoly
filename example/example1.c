//
//

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <GLFW/glfw3.h>
#pragma comment(lib, "glfw3.lib")

#define CDEC2D_IMPLEMENTATION
#include "cdec2d.h"

void drawframe(GLFWwindow* window)
{
	int width = 0, height = 0;

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	glClearColor(20.0f/255.0f, 20.0f/255.0f, 90.0f/255.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	//glEnable(GL_BLEND);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_TEXTURE_2D);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(-100.0f, 100.0f, -100.0f, 100.0f, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDisable(GL_DEPTH_TEST);
	glColor4ub(255,255,255,255);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	// Draw bounds
	glColor4ub(255,0,0,255);
	glBegin(GL_LINE_LOOP);
	glVertex2f(0, 0);
	glVertex2f(50.0f, 0);
	glVertex2f(50.0f, 50.0f);
	glVertex2f(0, 50.0f);
	glEnd();

	glfwSwapBuffers(window);
}

void resizecb(GLFWwindow* window, int width, int height)
{
	drawframe(window);
}

int main()
{
	GLFWwindow* window;
	const GLFWvidmode* mode;

	if (!glfwInit())
		return -1;

	mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
  window = glfwCreateWindow(1024, 768, "cdec2d", NULL, NULL);
	if (!window)
	{
		printf("Could not open window\n");
		glfwTerminate();
		return -1;
	}

	glfwSetFramebufferSizeCallback(window, resizecb);
	glfwMakeContextCurrent(window);
	//glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_SMOOTH);

	while (!glfwWindowShouldClose(window))
	{
		drawframe(window);
		glfwPollEvents();
	}

	glfwTerminate();
	return 0;
}
