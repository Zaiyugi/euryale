#include <iostream>
#include <fstream>
#include <cstdlib>

#include <cfloat>
#include <stdint.h>
#include <climits>

#include <random>
#include <string>
#include <vector>

#include "ogl/ogl_include.h"
#include <AntTweakBar.h>

#include "Vector.h"
#include "Matrix.h"
#include "LinearAlgebra.h"
#include "ColorGradient.h"

#include "SPHEngine.h"

#include "Bin.h"
#include "Array.h"

#include "Particle.h"

#define VIEW 0
#define EDIT 1

int initGLUT(int argc, char **argv);
void initGL();
void switchLightingState(bool);

// System variables
// START

/* ---------------------------------- */
int windowWidth, windowHeight;
int windowXPosition, windowYPosition;

GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];
GLfloat winX, winY, winZ;
GLdouble posX, posY, posZ;

bool update_display = true;
int update_timestep = 0;
int drawMode = 0;
int mouseInWindow = 0;

ogl::OGL_Camera* cam;
TwBar *main_bar, *tensor_bar;

size_t glyph_size;
lux::Vector* glyph_2d;
double glyph_scale = 1.0;

ColorGradient heatMapGradient;

uint64_t simulation_counter;
int steps_before_draw = 1;
int steps_until_resort = 5;
using FpSeconds = std::chrono::duration<float, std::chrono::seconds::period>;
FpSeconds elapsed_simulation_time(0);
/* ---------------------------------- */

// END

bool converged = false;
bool drawBins = false;
bool displayParticles = true;
bool displayGlyphs = false;
bool displayAnisotropic = false;
bool displayHighEnergy = false;
bool displaySingleGlyph = false;
bool use_superquadrics = false;

bool enableDebug = false;

double E_max;
double E_avg;
double edge_sharpness = 6.0;

int newParticles = 50;
size_t particle_cnt = 0;

double high, low;
double max_step;

// Math Constants
static const double robbins = 0.66170718;

Engine *engine;
std::string IO_Model_Name = "test";
std::string IO_Path_Name = "/home/zshore/Research/Academic/models";

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

lux::Vector OGLToLux(ogl::Vector3d V)
{ return lux::Vector(V[0], V[1], V[2]); }

lux::Vector unproject(int x, int y)
{
   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
   glGetDoublev( GL_PROJECTION_MATRIX, projection );
   glGetIntegerv( GL_VIEWPORT, viewport );

   winX = static_cast<float>(x);
   winY = static_cast<float>(viewport[3] - y);
   winZ = 1.0;

   gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ );

   return lux::Vector(posX, posY, posZ);
}

double sq_pow(double x, double a)
{
   return util::sgn(x) * std::pow(std::abs(x), a);
}

void draw()
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glPushMatrix();

   glPointSize(2.0);

   // Draw Particles
   float r, g, b;
   for(int it = 0; it < engine->_particles.size(); it++)
   {
      if(engine->_particles[it]->is_deleted)
         continue;

      lux::Vector p = engine->_particles[it]->p;
      double particleEnergy = engine->_particles[it]->E;

      glPushMatrix();
      glTranslatef(p[0], p[1], p[2]);

      double C = engine->getMappedParticleEnergy(it);
      heatMapGradient.getColorAtValue(C, r, g, b);

      int id_in_bin = engine->_particles[it]->id_in_bin;

      bool draw_particle = true;
      if( displayAnisotropic && T.is_isotropic() )
         draw_particle = false;
      if( displayHighEnergy && particleEnergy < engine->_E_avg)
         draw_particle = false;

      glColor4f(r, g, b, 1);
      if(displayGlyphs && draw_particle)
      {
   //       lux::Vector normal = engine->gradVolume(p).unitvector();
   //       lux::Vector rot_dir = lux::Vector(0, 0, 1) ^ normal;
   //       double rot_mag = std::acos(normal * lux::Vector(0, 0, 1));
   //       lux::Matrix rot = lux::rotation(rot_dir, rot_mag);
			// lux::Matrix full_transform = rot;

   //       glDisable(GL_LIGHTING);
   //       glBegin(GL_POLYGON);
   //          for(size_t k = 0; k < glyph_size; ++k)
   //          {
   //             lux::Vector glyph_trans = full_transform * glyph_2d[k] * engine->_repulsion;
   //             glVertex3dv( &(glyph_trans[0]) );
   //          }
   //       glEnd();
   //       glEnable(GL_LIGHTING);
         if( use_superquadrics )
            drawSuperquadric(engine->_particles[it]->T, engine->_repulsion * glyph_scale);
         else
            drawEllipsoid(engine->_particles[it]->T, engine->_repulsion * glyph_scale);
      }

      if(draw_particle)
      {
         switchLightingState(false);
         glColor4f(r, g, b, 1);
         glBegin(GL_POINTS);
            glVertex3f(0, 0, 0);
         glEnd();
         switchLightingState(true);
      }

      glTranslatef(-p[0], -p[1], -p[2]);
      glPopMatrix();
   }

   // Draw Tensors
   // for(int it = 0; it < engine->getTensorField()->size(); ++it)
   // {
   //    size_t n = engine->getTensorField()->tensors[it].size();
   //    for(size_t k = 0; k < n; ++k)
   //    {
   //       lux::Tensor* T = engine->getTensorField()->tensors[it][k];
   //       lux::Vector trans = T->pos;

   //       glPushMatrix();
   //       glTranslatef( trans[0], trans[1], trans[2] );

   //       glColor4f(1,1,1,1);

   //       glBegin(GL_LINE_LOOP);
   //          for(size_t k = 0; k < glyph_size; ++k)
   //          {
   //             lux::Vector glyph_trans = T->transform() * glyph_2d[k] * engine->_repulsion;
   //             glVertex3dv( &(glyph_trans[0]) );
   //          }
   //       glEnd();

   //       glTranslatef( -trans[0], -trans[1], -trans[2] );
   //       glPopMatrix();
   //    }

   // }

   // Draw Bins
   if(drawBins)
   {
      glColor4f(1.0, 1.0, 1.0, 1.0);
      for(int k = 0; k < engine->getBins()->getNz(); ++k)
      {
         for(int j = 0; j < engine->getBins()->getNy(); ++j)
         {
            for(int i = 0; i < engine->getBins()->getNx(); ++i)
            {
               if(engine->getBins()->eval(lux::Coord(i,j,k)).size() == 0)
                  continue;

               lux::Vector res = engine->getBins()->getRes().x();
               lux::Vector urc = engine->getBins()->evalP(lux::Coord(i+1, j+1, k+1));
               lux::Vector llc = engine->getBins()->evalP(lux::Coord(i, j, k));

               // Left
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &llc[0] );
                  glVertex3d( llc[0], llc[1], urc[2] );
                  glVertex3d( llc[0], urc[1], urc[2] );
                  glVertex3d( llc[0], urc[1], llc[2] );
               glEnd();

               // Bottom
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &llc[0] );
                  glVertex3d( urc[0], llc[1], llc[2] );
                  glVertex3d( urc[0], llc[1], urc[2] );
                  glVertex3d( llc[0], llc[1], urc[2] );
               glEnd();

               // Back
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &llc[0] );
                  glVertex3d( llc[0], urc[1], llc[2] );
                  glVertex3d( urc[0], urc[1], llc[2] );
                  glVertex3d( urc[0], llc[1], llc[2] );
               glEnd();

               // Right
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &urc[0] );
                  glVertex3d( urc[0], llc[1], urc[2] );
                  glVertex3d( urc[0], llc[1], llc[2] );
                  glVertex3d( urc[0], urc[1], llc[2] );
               glEnd();

               // Top
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &urc[0] );
                  glVertex3d( urc[0], urc[1], llc[2] );
                  glVertex3d( llc[0], urc[1], llc[2] );
                  glVertex3d( llc[0], urc[1], urc[2] );
               glEnd();

               // Front
               glBegin(GL_LINE_LOOP);
                  glVertex3dv( &urc[0] );
                  glVertex3d( urc[0], llc[1], urc[2] );
                  glVertex3d( llc[0], llc[1], urc[2] );
                  glVertex3d( llc[0], urc[1], urc[2] );
               glEnd();

            }
         }
      }
   }

   glPopMatrix();

   TwDraw();

   // engine->clearColors();

   glutSwapBuffers();

}

void generateGlyphs()
{
   glyph_size = 20;

   // 2D Glyph
   glyph_2d = new lux::Vector[glyph_size];
   for(size_t i = 0; i < glyph_size; ++i)
   {
      double t = i / static_cast<double>(glyph_size-1) * 2.0 * M_PI;
      lux::Vector p;
      p[0] = std::cos(t);
      p[1] = std::sin(t);
      p[2] = 0.0;
      glyph_2d[i] = p;
   }

}

// TwBar Callbacks

void TW_CALL SetMyStdStringCB(const void *value, void *)
{
   const std::string *srcPtr = static_cast<const std::string *>(value);
   IO_Model_Name = *srcPtr;
}

void TW_CALL GetMyStdStringCB(void *value, void *)
{
   std::string *destPtr = static_cast<std::string *>(value);
   TwCopyStdStringToLibrary(*destPtr, IO_Model_Name);
}

void TW_CALL ExportModel(void *)
{
   std::cout << "Exporting to Karl and Obj files\n";
   std::string karl_path = IO_Path_Name + "/" + IO_Model_Name + ".karl";
   std::string obj_path = IO_Path_Name + "/" + IO_Model_Name + ".obj";

   std::vector<lux::Vector> data;
   data.resize(engine->_particles.size());
   for(int i = 0; i < engine->_particles.size(); ++i)
      data[i] = engine->_particles[i]->p;

   util::writeKarl(karl_path, data);
   util::writeObj(obj_path, data);
   std::cout << "Export done\n";
}

void TW_CALL toggle_Bins(void *)
{
   drawBins = !drawBins;
}

void TW_CALL toggle_Particles(void *)
{
   displayParticles = !displayParticles;
}

void TW_CALL toggle_Glyphs(void *)
{
   displayGlyphs = !displayGlyphs;
}

void TW_CALL toggle_Anisotropic(void *)
{
   displayAnisotropic = !displayAnisotropic;
}

void TW_CALL toggle_HighEnergy(void *)
{
   displayHighEnergy = !displayHighEnergy;
}

void TW_CALL toggle_SingleGlyph(void *)
{
   displaySingleGlyph = !displaySingleGlyph;
}

void TW_CALL toggle_Superquadrics(void *)
{
   use_superquadrics = !use_superquadrics;
}

void TW_CALL addExtraParticles(void *)
{
   engine->injectParticles(newParticles);

   particle_cnt = engine->_particles.size();
}

void TW_CALL setRepulsion(const void *value, void *clientData)
{
   double repulsion = *static_cast<const double*>(value);
   engine->_repulsion = repulsion;
   engine->updateEnergyFunction();
}

void TW_CALL getRepulsion(void *value, void *clientData)
{
   *static_cast<double*>(value) = engine->_repulsion;
}

void TW_CALL setGlyphScaling(const void *value, void *clientData)
{
   glyph_scale = *static_cast<const double*>(value);
}

void TW_CALL getGlyphScaling(void *value, void *clientData)
{
   *static_cast<double*>(value) = glyph_scale;
}

void TW_CALL setEdgeSharpness(const void *value, void *clientData)
{
   edge_sharpness = *static_cast<const double*>(value);
}

void TW_CALL getEdgeSharpness(void *value, void *clientData)
{
   *static_cast<double*>(value) = edge_sharpness;
}

/* --------------------- */

// OpenGL Callbacks
// START

void ogl_cleanup()
{
   delete [] glyph_2d;

   if(benchmarking)
      if(benchmarking_ofs.is_open())
         benchmarking_ofs.close();

   TwTerminate();
}

void timerCB(int millisec)
{
   glutTimerFunc(millisec, timerCB, millisec);

   if(update_display && !converged)
   {
		auto clock_start = std::chrono::high_resolution_clock::now();

      for(int i = 0; i < steps_before_draw; ++i)
      {
         engine->distribute();
         simulation_counter++;

         if(simulation_counter % steps_until_resort == 0)
         {
            engine->updateBins();
         }
      }

		auto clock_end = std::chrono::high_resolution_clock::now();

		FpSeconds elapsed_sec = FpSeconds(clock_end - clock_start);
      elapsed_simulation_time += elapsed_sec;

      if( benchmarking )
      {
         std::streambuf* backup;
         if( benchmarking_ofs.is_open() )
         {
            backup = std::cout.rdbuf();
            std::cout.rdbuf(benchmarking_ofs.rdbuf());
         }

         std::cout << simulation_counter << ", " << engine->_E_avg << ", " << engine->_E_max << ", ";
         std::cout << engine->_avg_stencil_size << ", " << engine->_avg_neighbors << ", ";
         std::cout << elapsed_sec.count() << ", " << elapsed_simulation_time.count() << "\n";

         if( benchmarking_ofs.is_open() )
            std::cout.rdbuf(backup);

         if( simulation_counter > benchmarking_iterations && benchmarking_iterations > 0)
         {
            std::cout << flush;
            ogl_cleanup();
            exit(0);
         }

      } else {
         std::cout << "simulation_counter: " << simulation_counter << "\n";
         std::cout << "\t>>> processed " << steps_before_draw << " iterations, resorting every " << steps_until_resort << " steps\n";
         std::cout << "\t    >>> Avg. Energy: " << engine->_E_avg << " | Max Energy: " << engine->_E_max << "\n";
         std::cout << "\t    >>> Avg. Stencil Size: " << engine->_avg_stencil_size << " | Avg. Neighbors: " << engine->_avg_neighbors << "\n";
         std::cout << "\t    >>> elapsed time for last call: " << elapsed_sec.count() << " seconds\n";
         std::cout << "\t    >>> elapsed simulation time: " << elapsed_simulation_time.count() << " seconds\n";
               
      }
      std::cout << flush;
   }

   glutPostRedisplay();
}

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

   TwWindowSize(windowWidth, windowHeight);
}

void keyboardCB(unsigned char key, int x, int y)
{
   if( !TwEventKeyboardGLUT(key, x, y) )
   {
      switch(key)
      {
         case 27: // ESCAPE
            ogl_cleanup();
            exit(0);
            break;

         case 'r': // switch rendering modes (fill -> wire -> point)
         case 'R':
            drawMode = ++drawMode % 3;
            if(drawMode == 0)        // fill mode
            {
               glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
               glEnable(GL_DEPTH_TEST);
            }
            else if(drawMode == 1)  // wireframe mode
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
   }

   glutPostRedisplay();
}

void mouseEntryCB(int state)
{
   if(state == GLUT_ENTERED)
   {
      mouseInWindow = state;
   } else if(state == GLUT_LEFT) {
      mouseInWindow = state;
   }

}

void mouseMotionCB(int x, int y)
{
   if( !TwEventMouseMotionGLUT(x, y) )
   {
      cam->HandleMouseMotion(x, y);
   }

   glutPostRedisplay();

}

void mouseCB(int button, int state, int x, int y)
{
   if( !TwEventMouseButtonGLUT(button, state, x, y) )
   {
      cam->HandleMouseEvent(button, state, x, y);
   }

   glutPostRedisplay();
}

// MAIN

int main(int argc, char** argv)
{
   std::cout.sync_with_stdio(false);

   std::string parameter_file, benchmarking_file("cout");
   switch(argc)
   {
      case 5:
         std::cout << "processing argument 4\n";
         benchmarking_file = std::string(argv[4]);
      case 4:
         std::cout << "processing argument 3\n";
         benchmarking_iterations = std::atoi(argv[3]);
      case 3:
         std::cout << "processing argument 2\n";
         benchmarking = std::atoi(argv[2]);
      case 2:
         std::cout << "processing argument 1\n";
         parameter_file = std::string(argv[1]);
         break;

      default:
         std::cerr << "sample <parameter file> [benchmark flag]\n";
         exit(0);
   }
   std::string research_path("/home/zshore/Research/Academic");
   std::string project_base_path = research_path + std::string("/repos/aniso-meshing.sampling");

   std::cout << "Got parameter_file: " << parameter_file << "\n\n";

   simulation_counter = 0;

   if( benchmarking )
      if( benchmarking_file != "cout" )
         benchmarking_ofs.open(benchmarking_file.c_str(), std::ofstream::out);

   /* Start Parameters */

   int initial_size;
   double bins_per_side;
   double repulsion, lambda_max, lambda_min, initial_lambda, ideal_energy;

   std::ifstream parameter_stream(parameter_file, std::ifstream::in);

   std::string junk;
   std::string nrrd_file;

   parameter_stream >> junk >> nrrd_file;
   parameter_stream >> junk >> repulsion;
   parameter_stream >> junk >> lambda_max;
   parameter_stream >> junk >> lambda_min;
   parameter_stream >> junk >> initial_lambda;
   parameter_stream >> junk >> initial_size;
   parameter_stream >> junk >> ideal_energy;
   parameter_stream >> junk >> steps_before_draw;
   parameter_stream >> junk >> steps_until_resort;
   parameter_stream >> junk >> bins_per_side;

   std::cout << "dataset: " << nrrd_file << "\n";
   std::cout << "repulsion: " << repulsion << "\n";
   std::cout << "lambda_max: " << lambda_max << "\n";
   std::cout << "lambda_min: " << lambda_min << "\n";
   std::cout << "initial_lambda: " << initial_lambda << "\n";
   std::cout << "initial_size: " << initial_size << "\n";
   std::cout << "ideal_energy: " << ideal_energy << "\n";
   std::cout << "steps_before_draw: " << steps_before_draw << "\n";
   std::cout << "steps_until_resort: " << steps_until_resort << "\n";
   std::cout << "bins_per_side: " << bins_per_side << "\n";

   parameter_stream.close();

   /* End Parameters */

   double high = 0.5;
   double low = -0.5;

   lux::Vector URC(high, high, high);
   lux::Vector LLC(low, low, low);

   heatMapGradient.loadFromCSV(project_base_path + std::string("/colormaps/cubeyf1.csv"), 16);

   /* ~ Tensors ~ */
   Nrrd *nin; nin = nrrdNew();

   nrrd_file = research_path + "/datasets/" + nrrd_file;

   char *err;
   if(nrrdLoad(nin, nrrd_file.c_str(), NULL))
   {
      err = biffGetDone(NRRD);
      std::cerr << "Trouble reading file \"" << nrrd_file << "\" : " << err << "\n";
      delete err;
      exit(0);
   }

   std::cout << nrrd_file << " is a " << nin->dim << "-dimensional nrrd of type " << nin->type << " (" << airEnumStr(nrrdType, nin->type) << ")\n";

   size_t *axis_size = new size_t[nin->dim-1];
   double *axis_spacing = new double[nin->dim-1];
   for(size_t i = 1; i < nin->dim; i++)
   {
      axis_size[i-1] = nin->axis[i].size;
      axis_spacing[i-1] = (high - low) / static_cast<double>(axis_size[i-1]);

      std::cout << "Axis " << i << " > size: " << axis_size[i-1] << " | spacing: " << axis_spacing[i-1] << "\n";
   }

   lux::TensorField* tensor_field = new lux::TensorField(URC, LLC, lux::Vector(axis_spacing[0], axis_spacing[1], axis_spacing[2]));

   // Index scheme: 
   // size_t ndx = m + a2_size * (k + a1_size * (j + a0_size * i));
   lux::Tensor* t;
   size_t ndx = 0;
   for(size_t k = 0; k < axis_size[2]; ++k)
      for(size_t j = 0; j < axis_size[1]; ++j)
         for(size_t i = 0; i < axis_size[0]; ++i)
         {
            float confidence = static_cast<float*>(nin->data)[ndx];

            lux::Vector pos;
            pos[0] = i * axis_spacing[0] + low;
            pos[1] = j * axis_spacing[1] + low;
            pos[2] = k * axis_spacing[2] + low;
            if(confidence < 1.0)
            {
               ndx += 7;
               
               t = new lux::Tensor;
               t->pos = pos;
               t->setValues(1, 1, 1);
               t->setVectors(lux::EUCL3D_X_AXIS, lux::EUCL3D_Y_AXIS, lux::EUCL3D_Z_AXIS);
               tensor_field->addTensor(t);

               continue;
            }

            float Dxx = static_cast<float*>(nin->data)[ndx + 1];
            float Dxy = static_cast<float*>(nin->data)[ndx + 2];
            float Dxz = static_cast<float*>(nin->data)[ndx + 3];
            float Dyy = static_cast<float*>(nin->data)[ndx + 4];
            float Dyz = static_cast<float*>(nin->data)[ndx + 5];
            float Dzz = static_cast<float*>(nin->data)[ndx + 6];

            lux::Vector v1(Dxx, Dxy, Dxz);
            lux::Vector v2(Dxy, Dyy, Dyz);
            lux::Vector v3(Dxz, Dyz, Dzz);

            double e1 = v1.magnitude();
            double e2 = v2.magnitude();
            double e3 = v3.magnitude();

            t = new lux::Tensor;
            t->pos = pos;
            t->setValues(e1, e2, e3);
            t->setVectors(v1.unitvector(), v2.unitvector(), v3.unitvector());
            tensor_field->addTensor(t);

            ndx += 7;
         }

         
   nrrdNuke(nin);
   delete [] axis_size;
   delete [] axis_spacing;

   /* Create Engine */
   // To create an array with 1 cell, 
   // use URC - LLC as the resolution
   // 3-D Bin setup
   lux::Array<Bin>* space = new lux::Array<Bin>(URC, LLC, (URC - LLC)/bins_per_side, Bin());

	engine = new Engine(URC, LLC, repulsion, lambda_max, lambda_min, initial_lambda, ideal_energy, tensor_field, space);
   std::cout.flush();

   /* ~ Particles ~ */
   engine->injectParticles(initial_size);
   engine->computeEminEmax();
   E_max = engine->_E_max;
   E_avg = engine->_E_avg;
   particle_cnt = engine->_particles.size();

   std::cout << "Initial E_max: " << E_max << "\n";
   std::cout << "Initial E_avg: " << E_avg << "\n";

   for(size_t i = 0; i < engine->_particles.size(); i++)
   {
      if(engine->_particles[i]->E < 5e-8)
         engine->_particles[i]->color = lux::Vector(1, 1, 1);
      else {
         double C = engine->getMappedParticleEnergy(i);
         float r, g, b;
         heatMapGradient.getColorAtValue(C, r, g, b);
         engine->_particles[i]->color = lux::Vector(r, g, b);
      }
   }

   /* ~ User Interface ~ */
   windowWidth = 960;
   windowHeight = 540;
   windowXPosition = windowYPosition = 0;
   update_timestep = 25;
   update_display = true;

   cam = new ogl::OGL_Camera(ogl::Vector3d(0, 0, 1), ogl::Vector3d(0, 0, 0), ogl::Vector3d(0, 1, 0));
   cam->SetFOV(60.0);

   initGLUT(argc, argv);
   initGL();

   // Setup TwBar
   TwInit(TW_OPENGL, NULL);
   TwWindowSize(windowWidth, windowHeight);

   /* ~ MAIN CONTROLS ~ */
   main_bar = TwNewBar("General");
   TwDefine(" 'General' size='300 280' valueswidth=fit contained=true ");

   TwAddVarCB(main_bar, "Filename", TW_TYPE_STDSTRING, SetMyStdStringCB, GetMyStdStringCB, NULL, " group='Main' ");
   TwAddButton(main_bar, "Export", ExportModel, NULL, " group='Main' ");

   TwAddSeparator(main_bar, "DEBUG_Sep", " group='Main' ");

   TwAddButton(main_bar, "Display particles", toggle_Particles, NULL, " group='Main' ");
   TwAddButton(main_bar, "Display bins", toggle_Bins, NULL, " group='Main' ");
   TwAddButton(main_bar, "Dispaly glyphs", toggle_Glyphs, NULL, " group='Main' ");
   TwAddButton(main_bar, "Show only anisotropic", toggle_Anisotropic, NULL, " group='Main' ");
   TwAddButton(main_bar, "Show only high energy (>E_avg)", toggle_HighEnergy, NULL, " group='Main' ");
   TwAddButton(main_bar, "Toggle superquadrics", toggle_Superquadrics, NULL, " group='Main' ");
   TwAddVarCB(main_bar, "Glyph Scale", TW_TYPE_DOUBLE, setGlyphScaling, getGlyphScaling, NULL, " group='Main' ");
   TwAddVarCB(main_bar, "Edge Sharpness", TW_TYPE_DOUBLE, setEdgeSharpness, getEdgeSharpness, NULL, " group='Main' ");

   TwAddSeparator(main_bar, "Display_Sep", " group='Main' ");

   TwAddVarRO(main_bar, "Particle Count", TW_TYPE_INT32, &particle_cnt, " group='Main' ");
   TwAddButton(main_bar, "Add particles", addExtraParticles, NULL, " group='Main' ");
   TwAddVarRW(main_bar, "Additional Particles", TW_TYPE_INT32, &newParticles, " group='Main' ");

   TwAddSeparator(main_bar, "Edit_Sep", " group='Main' ");

   TwAddVarCB(main_bar, "Repulsion", TW_TYPE_DOUBLE, setRepulsion, getRepulsion, NULL, " group='Main' ");
   TwAddVarRW(main_bar, "Max. Lambda", TW_TYPE_DOUBLE, &engine->_lambda_max, " group='Main' ");
   TwAddVarRW(main_bar, "Min. Lambda", TW_TYPE_DOUBLE, &engine->_lambda_min, " group='Main' ");
   TwAddVarRW(main_bar, "Ideal Energy (per particle)", TW_TYPE_DOUBLE, &engine->_E_ideal, " group='Main' ");
   TwAddVarRW(main_bar, "Max step size", TW_TYPE_DOUBLE, &max_step, " group='Main' ");

   // Create Glyph
   generateGlyphs();

   glutMainLoop();

return 0;
}

// INITIALIZE

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
   glutTimerFunc(update_timestep, timerCB, update_timestep);

   glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
   TwGLUTModifiersFunc(glutGetModifiers);

return handle;
}

void initLights()
{
   glLightModelf(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

   float highf = 4.0f * static_cast<float>(engine->_high[0]);
   float lowf = 4.0f * static_cast<float>(engine->_low[0]);

   GLfloat light_Ka[] = {0.5f, 0.5f, 0.5f, 1.0f};


   float light0_Pos[4] = {highf, 0.f, 0.f, 1};
   float light0_Dir[3] = {-1.f, 0.f, 0.f};

   float light1_Pos[4] = {0.f, highf, 0.f, 1};
   float light1_Dir[3] = {0.f, -1.f, 0.f};

   float light2_Pos[4] = {0.f, 0.f, highf, 1};
   float light2_Dir[3] = {0.f, 0.f, -1.f};


   float light3_Pos[4] = {lowf, 0.f, 0.f, 1};
   float light3_Dir[3] = {1.f, 0.f, 0.f};

   float light4_Pos[4] = {0.f, lowf, 0.f, 1};
   float light4_Dir[3] = {0.f, 1.f, 0.f};

   float light5_Pos[4] = {0.f, 0.f, lowf, 1};
   float light5_Dir[3] = {0.f, 0.f, 1.f};


   glLightfv(GL_LIGHT0, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);
   glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light0_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);

   glLightfv(GL_LIGHT1, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);
   glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, light1_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);

   glLightfv(GL_LIGHT2, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT2, GL_POSITION, light2_Pos);
   glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light2_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);


   glLightfv(GL_LIGHT3, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT3, GL_POSITION, light3_Pos);
   glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, light3_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);

   glLightfv(GL_LIGHT4, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT4, GL_POSITION, light4_Pos);
   glLightfv(GL_LIGHT4, GL_SPOT_DIRECTION, light4_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);

   glLightfv(GL_LIGHT5, GL_DIFFUSE, light_Ka);
   glLightfv(GL_LIGHT5, GL_POSITION, light5_Pos);
   glLightfv(GL_LIGHT5, GL_SPOT_DIRECTION, light5_Dir);
   glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.0f);


   switchLightingState(true);
}

void switchLightingState(bool state)
{
   if(state)
   {
      glEnable(GL_LIGHTING);
      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      glEnable(GL_LIGHT2);
      // glEnable(GL_LIGHT3);
      // glEnable(GL_LIGHT4);
      // glEnable(GL_LIGHT5);

   } else {
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
      glDisable(GL_LIGHT1);
      glDisable(GL_LIGHT2);
      // glDisable(GL_LIGHT3);
      // glDisable(GL_LIGHT4);
      // glDisable(GL_LIGHT5);
   }

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

   // BG Color
   glClearColor(0.1, 0.1, 0.1, 0);

   // 0 is near, 1 is far
   glClearDepth(1.0f);
   glDepthFunc(GL_LEQUAL);

   initLights();
}
