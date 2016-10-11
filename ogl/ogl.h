/* Zachary Shore
 * 2014-01-03
 * OpenGL Functions and Callbacks
 */

#ifndef __OGL_H__
#define __OGL_H__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "FrameGrabber.h"
#include "Matrix.h"
#include "Camera.h"

#define PI180 M_PI/180.0

#define FACE_FILL	0
#define WIREFRAME	1
#define POINTS		2

using namespace std;

namespace ogl
{

int windowWidth, windowHeight;
int windowXPosition, windowYPosition;

/* ---------------------------------- */
bool update_display = false;
int update_timestep = 25;
int drawMode = 0;
int mouseInWindow = 0;

OGL_Camera* cam;
/* ---------------------------------- */

void draw();
//void user_cleanup();

void perspectiveDisplay()
{
   cam->PerspectiveDisplay(windowWidth, windowHeight);
}

void setRenderingProperties(int w, int h, int x, int y, int update_dt)
{ 
   windowWidth = w;
   windowHeight = h;
   windowXPosition = x;
   windowYPosition = y;
   update_timestep = update_dt;
   update_display = true;
}

void setCamera(Vector3d eye, Vector3d view, Vector3d up)
{
   cam = new OGL_Camera(eye, view, up);
}

void ogl_cleanup() { delete cam; }

void IdleCB(){}

void timerCB(int millisec);
/*
{
   glutTimerFunc(millisec, timerCB, millisec);

   glutPostRedisplay();
}*/

void reshapeCB(int w, int h)
{
   windowWidth = w;
   windowHeight = h;
   glViewport(0, 0, (GLsizei)w, (GLsizei)h);

   //float aspectRatio = (float)w / h;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glMatrixMode(GL_MODELVIEW);
}

void keyboardCB(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 27: // ESCAPE
        //user_cleanup();
        ogl_cleanup();
        exit(0);
        break;

    case 'r': // switch rendering modes (fill -> wire -> point)
    case 'R':
        drawMode = ++drawMode % 3;
        if(drawMode == FACE_FILL)        // fill mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
        }
        else if(drawMode == WIREFRAME)  // wireframe mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_DEPTH_TEST);
        }
        else                    // point mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
            glDisable(GL_DEPTH_TEST);
        }
        break;

    case 'p':
    case 'P':
        update_display = !update_display;
        break;

    default:
        ;
    }

    glutPostRedisplay();
}

void mouseEntryCB(int state)
{
   if(state == GLUT_ENTERED)
   {
      /*printf("------------------------------\n");
      printf("!! MOUSE ENTERED WINDOW %d !!\n",state);
      printf("------------------------------\n");*/
      mouseInWindow = state;
   } else if(state == GLUT_LEFT) {
      /*printf("------------------------------\n");
      printf("!! MOUSE LEFT WINDOW %d !!\n", state);
      printf("------------------------------\n");*/
      mouseInWindow = state;
   }

}

void mouseMotionCB(int x, int y)
{
   cam->HandleMouseMotion(x, y);
   glutPostRedisplay();

}

void mouseCB(int button, int state, int x, int y)
{
   cam->HandleMouseEvent(button, state, x, y);
}

int initGLUT(int argc, char **argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
   glutInitWindowSize(windowWidth, windowHeight);
   glutInitWindowPosition(windowXPosition, windowYPosition);
   int handle = glutCreateWindow(argv[0]);

   glutDisplayFunc(draw);
   glutReshapeFunc(reshapeCB);
   glutKeyboardFunc(keyboardCB);
   glutMouseFunc(mouseCB);
   glutMotionFunc(mouseMotionCB);
   glutEntryFunc(mouseEntryCB);
   glutIdleFunc(IdleCB);
   glutTimerFunc(update_timestep, timerCB, update_timestep);

return handle;
}

void initLights()
{
   glLightModelf(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

   GLfloat light0_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light0_Kd[] = {0.7f, 0.0f, 0.0f, 1.0f};  // diffuse light
   GLfloat light0_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT0, GL_AMBIENT, light0_Ka);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_Kd);
   glLightfv(GL_LIGHT0, GL_SPECULAR, light0_Ks);
   float light0_Pos[4] = {15, 0, 0, 1}; // positional light
   glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);

   GLfloat light1_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light1_Kd[] = {0.0f, 0.7f, 0.0f, 1.0f};  // diffuse light
   GLfloat light1_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT1, GL_AMBIENT, light1_Ka);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_Kd);
   glLightfv(GL_LIGHT1, GL_SPECULAR, light1_Ks);
   float light1_Pos[4] = {-15, 0, 0, 1}; // positional light
   glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);

   GLfloat light2_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light2_Kd[] = {0.0f, 0.0f, 0.7f, 1.0f};  // diffuse light
   GLfloat light2_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT2, GL_AMBIENT, light2_Ka);
   glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_Kd);
   glLightfv(GL_LIGHT2, GL_SPECULAR, light2_Ks);
   float light2_Pos[4] = {0, 0, 15, 1}; // positional light
   glLightfv(GL_LIGHT2, GL_POSITION, light2_Pos);

   GLfloat light3_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light3_Kd[] = {0.7f, 0.7f, 0.7f, 1.0f};  // diffuse light
   GLfloat light3_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT3, GL_AMBIENT, light3_Ka);
   glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_Kd);
   glLightfv(GL_LIGHT3, GL_SPECULAR, light3_Ks);
   float light3_Pos[4] = {0, 0, -15, 1}; // positional light
   glLightfv(GL_LIGHT3, GL_POSITION, light3_Pos);

   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_LIGHT1);
   glEnable(GL_LIGHT2);
   glEnable(GL_LIGHT3);
}

void initGL()
{
   glShadeModel(GL_SMOOTH);

   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
   glEnable(GL_DEPTH_TEST);

   //glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);

   glClearColor(0, 0, 0, 0);                   // background color
   glClearDepth(1.0f);                         // 0 is near, 1 is far
   glDepthFunc(GL_LEQUAL);

   //initLights();
}

}

#endif
