#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>
#include <set>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
using namespace std;

/* --------------------------------- */
void displayCB();
void reshapeCB(int w, int h);
void keyboardCB(unsigned char key, int x, int y);
void mouseCB(int button, int stat, int x, int y);
void mouseMotionCB(int x, int y);
void mouseEntryCB(int state);
void IdleCB();

void initGL();
int  initGLUT(int argc, char **argv);
void initLights();

void editMenu(int value);
void setCamera();
/* ---------------------------------- */
