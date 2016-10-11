/*
* Final SPH2D Project
* SPH with material editing and fluid vis.
*/
#include <cstdlib>
#include <cstdio>
#include <cfloat>
#include <climits>
#include <cstring>
#include <omp.h>
#include <string>
#include <vector>
#include <utility>

#include "ogl/ogl_include.h"
#include <OpenImageIO/imageio.h>

#include "Image.h"
#include "OIIOFiles.h"

#include "Vector.h"
#include "ColorGradient.h"
#include "CmdLineFind.h"
#include "SPHEngine.h"

#define OFFSET_BUFFER(offset) ((char*)NULL + offset)

OIIO_NAMESPACE_USING;

int initGLUT(int argc, char **argv);
void initGL();

// OpenGL variables
// START

/* ---------------------------------- */
int window_width, window_height;
int window_x_position, window_y_position;

GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];
GLfloat win_x, win_y, win_z;
GLdouble pos_x, pos_y, pos_z;

int draw_mode = 0;
int mouse_in_window = 0;

size_t img_width, img_height, img_depth;
float* master_image = nullptr;
float scale_factor = 1.0f;

GLuint texture_id;
GLuint pbo_id;

int prev_mx = 0, prev_my = 0;
float mv_x, mv_y;

ColorGradient heatmapGradient;
bool use_velocity_for_color = false;
bool draw_bins = false;
bool draw_glyphs = false;
bool draw_density = true;
bool draw_particles = true;
bool paused = true;
std::vector<lux::Vector2d> glyph;

// Frame recording
bool recording = false;
size_t frame_number = 1;
size_t frame_max;
std::string frame_store;
std::string frame_id;

// SPH
float timestep = 1.0 / 24.0;
int nparticles;
int particles_to_inject;
int particles_in_system = 0;
int system_step = 0;
int steps_til_inject = 1;

sim::SPHEngine* engine = nullptr;
/* ---------------------------------- */

// END

std::random_device global_rd;
std::mt19937 global_gen(global_rd());
std::uniform_real_distribution<> global_dis(0, 1);

void draw()
{
   glClear(GL_COLOR_BUFFER_BIT);

   glPushMatrix();

   lux::Vector2d urc = engine->_urc;
   lux::Vector2d llc = engine->_llc;

   if(draw_density)
   {
      glEnable(GL_TEXTURE_2D);
      glBindTexture(GL_TEXTURE_2D, texture_id);
      glBegin(GL_TRIANGLES);
         // 1st tri
         glTexCoord2f(0.0, 0.0);
         glVertex3f(llc[0], llc[1], 0.0);

         glTexCoord2f(0.0, 1.0);
         glVertex3f(llc[0], urc[1], 0.0);

         glTexCoord2f(1.0, 1.0);
         glVertex3f(urc[0], urc[1], 0.0);

         // 2nd tri
         glTexCoord2f(1.0, 1.0);
         glVertex3f(urc[0], urc[1], 0.0);

         glTexCoord2f(1.0, 0.0);
         glVertex3f(urc[0], llc[1], 0.0);

         glTexCoord2f(0.0, 0.0);
         glVertex3f(llc[0], llc[1], 0.0);
      glEnd();
      glBindTexture(GL_TEXTURE_2D, 0);
      glDisable(GL_TEXTURE_2D);
   }

   if(draw_particles)
   {
      glPointSize(2.0);
      for(auto particle : engine->_particles)
      {
         lux::Vector2d x = particle->_p;
         lux::Vector rgb = particle->_mat->_color;
         if( use_velocity_for_color )
         {
            float C = particle->_v.magnitude() / engine->getStats().max_velocity;
            heatmapGradient.getColorAtValue(C, rgb[0], rgb[1], rgb[2]);
         }

         glPushMatrix();
         glTranslatef(x[0], x[1], 0);

         glColor3dv( &(rgb[0]) );
         if( draw_glyphs )
         {
            glBegin(GL_LINE_LOOP);
            for(size_t i = 0; i < glyph.size(); ++i)
               glVertex3dv( &(glyph[i][0]) );
            glEnd();
         }

         glBegin(GL_POINTS);
            glVertex3f(0.0, 0.0, 0.0);
         glEnd();

         glTranslatef(-x[0], -x[1], 0);
         glPopMatrix();
      }
   }

   if( draw_bins )
   {
      for(int i = 0; i < engine->getVolume()->getNx(); ++i)
      {
         for(int j = 0; j < engine->getVolume()->getNy(); ++j)
         {
            if( engine->getVolume()->eval(i, j).size() != 0)
               glColor3d(1.0, 0.0, 0.0);
            else
               glColor3d(1.0, 1.0, 1.0);

            lux::Vector2d urc = engine->getVolume()->evalP(i+1, j+1);
            lux::Vector2d llc = engine->getVolume()->evalP(i, j);

            glBegin(GL_LINE_LOOP);
               glVertex3f(llc[0], llc[1], 0);
               glVertex3f(llc[0], urc[1], 0);
               glVertex3f(urc[0], urc[1], 0);
               glVertex3f(urc[0], llc[1], 0);
            glEnd();
         }
      }
   }

   glPopMatrix();

   glColor3d(1.0, 1.0, 1.0);

   glutSwapBuffers();
}

// Utility

void readImage(const char* fname)
{
   ImageInput *in = ImageInput::create (fname);
   if( !in ) { return; }

   ImageSpec spec;
   in->open (fname, spec);
   img_width = spec.width;
   img_height = spec.height;
   img_depth = spec.nchannels;

   master_image = new float[img_width*img_height*img_depth];
   float* pixels = new float[img_width*img_height*img_depth];
   in->read_image (TypeDesc::FLOAT, pixels);

   for(size_t i = 0; i < img_width; ++i)
      for(size_t j = 0; j < img_height; ++j)
         for(size_t k = 0; k < img_depth; ++k)
         {
            int ndx = (i + img_width * j) * img_depth + k;
            int inv = (i + img_width * (img_height - j - 1)) * img_depth + k;
            master_image[ndx] = pixels[inv];
         }

   in->close ();
   delete in;
   delete pixels;

}

void createGlyph(const int sections)
{
   glyph.clear();
   glyph.resize(sections);
   float angle = 2.0 * M_PI / static_cast<float>(sections);
   for(int i = 0; i < sections; ++i)
   {
      glyph[i][0] = engine->getRadius() * std::cos(i * angle);
      glyph[i][1] = engine->getRadius() * std::sin(i * angle);
   }
}

void recordCurrentFrame()
{
   float* img = new float[window_width * window_height * 3];
   glReadPixels(0, 0, window_width, window_height, GL_RGB, GL_FLOAT, img);

   for(int i = 0; i < window_width; i++)
      for(int j = 0, k = window_height-1; j < window_height/2 && k >= window_height/2; j++, k--)
      {
         for(int c = 0; c < 3; ++c)
            std::swap(img[3 * (j * window_width + i) + c], img[3 * (k * window_width + i) + c]);
      }

   // Build filename
   std::string fnum = std::to_string(frame_number);
   if(frame_number < 1000)
      fnum = "0" + fnum;
   if(frame_number < 100)
      fnum = "0" + fnum;
   if(frame_number < 10)
      fnum = "0" + fnum;
   std::string full_frame_id = frame_store + "/" + frame_id + "." + fnum + ".png";

   // Write Image to file
   lux::writeOIIOImage(full_frame_id.c_str(), img, window_width, window_height, 3);

   delete [] img;
}

void cleanUp()
{
   if( engine )
      delete engine;

   if( master_image )
      delete [] master_image;
}

// Textures

void loadTextures()
{
   // Create Texture
   glGenTextures(1, &texture_id);
   glBindTexture(GL_TEXTURE_2D, texture_id);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,     GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,     GL_CLAMP_TO_EDGE);

   // Setup PBO
   glGenBuffers(1, &pbo_id);
   glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id);

   // Stream texture data to GPU
   size_t width = engine->_Nx;
   size_t height = engine->_Ny;
   size_t buffer_size = width * height * 3 * sizeof(float);
   glBufferData(GL_PIXEL_UNPACK_BUFFER, buffer_size, NULL, GL_STREAM_DRAW);
   glBufferSubData(GL_PIXEL_UNPACK_BUFFER, 0, buffer_size, engine->getColorField());

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_FLOAT, OFFSET_BUFFER(0));

   glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
   glBindTexture(GL_TEXTURE_2D, 0);

}

void updateTextures()
{
   // Stream updated texture to GPU
   glBindTexture(GL_TEXTURE_2D, texture_id);
   glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_id);

   size_t width = engine->_Nx;
   size_t height = engine->_Ny;
   size_t buffer_size = width * height * 3 * sizeof(float);
   glBufferSubData(GL_PIXEL_UNPACK_BUFFER, 0, buffer_size, engine->getColorField());

   glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGB, GL_FLOAT, OFFSET_BUFFER(0));

   glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
   glBindTexture(GL_TEXTURE_2D, 0);

}

// Materials

void handleMaterialEditing()
{
   std::cout << std::endl
      << "### Materials Menu ###" << std::endl
      << "========================================" << std::endl;

   size_t option = 0;
   std::cout << "(0) Select or (1) Create" << std::endl << "~> ";
   std::cin >> option;
   std::cout << "========================================" << std::endl;

   if( option == 0 )
   {
      do
      {
         std::cout << "Select Material: " << std::endl;
         std::cout << "----------------------------------------" << std::endl;
         for(size_t i = 1; i < engine->_materials.size(); i++)
            std::cout << "(" << i << ") " << engine->_materials[i]->_name << std::endl;
         std::cout << std::endl << "~> ";
         std::cin >> option;
      } while( option < 1 || option >= engine->_materials.size() );

      engine->_current_material = option;

   } else if( option == 1 )
   {
      std::string name;
      lux::Vector C;
      float m, bd, ps, pg, v, ve;

      std::cout << "Create Material: " << std::endl;
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Name: ";
      std::cin.ignore();
      std::getline(cin, name);
      std::cout << "Color: "; std::cin >> C[0] >> C[1] >> C[2];
      std::cout << "Mass: "; std::cin >> m;
      std::cout << "Reference Density: "; std::cin >> bd;
      std::cout << "Pressure Scale: "; std::cin >> ps;
      std::cout << "Pressure Gamma: "; std::cin >> pg;
      std::cout << "Viscosity: "; std::cin >> v;
      std::cout << "Viscosity Epsilon: "; std::cin >> ve;

      engine->genMaterial(C, m, bd, ps, pg, v, ve);
      engine->_materials.back()->_name = name;
      engine->_current_material = engine->_materials.size() - 1;
   }

   std::cout << "----------------------------------------" << std::endl;
   engine->printCurrentMaterialAttributes();
   std::cout << "========================================" << std::endl;

}

// OpenGL Callbacks
// START

void idleCallback()
{
   if( !paused )
   {
      if(!(system_step % 24))
         std::cout << "Iter: " << system_step << std::endl;

      // Kazam!
      engine->update();
      updateTextures();

      if(recording)
      {
         if( frame_number > frame_max )
         {
            std::cout << "Exceeded max frame count. Exiting..." << std::endl;
            exit(0);
         }

         recordCurrentFrame();
         frame_number++;
      }

      system_step += 1;
   }

   glutPostRedisplay();
}

void reshapeCallback(int w, int h)
{
   window_width = w;
   window_height = h;
   glViewport(0, 0, (GLsizei)w, (GLsizei)h);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   glOrtho(
      engine->_llc[0], engine->_urc[0],
      engine->_llc[1], engine->_urc[1],
      -1.0, 1.0
   );

   glMatrixMode(GL_MODELVIEW);

}

void keyboardCallback(unsigned char key, int x, int y)
{
   switch(key)
   {
      case 27: // ESCAPE
         exit(0);
         break;

      case 'r':
         // Reset some stuff
         engine->reset();
         break;

      case '1':
         use_velocity_for_color = false;
         std::cout << "Toggled: use_velocity_for_color -> " << (use_velocity_for_color ? "ON" : "OFF") << std::endl;
         break;

      case '2':
         use_velocity_for_color = true;
         std::cout << "Toggled: use_velocity_for_color -> " << (use_velocity_for_color ? "ON" : "OFF") << std::endl;
         break;

      case '4':
         draw_density = !draw_density;
         std::cout << "Toggled: draw_density -> " << (draw_density ? "ON" : "OFF") << std::endl;
         break;

      case '5':
         draw_particles = !draw_particles;
         std::cout << "Toggled: draw_particles -> " << (draw_particles ? "ON" : "OFF") << std::endl;
         break;

      case 'P':
         paused = !paused;
         std::cout << "Toggled: pause -> " << (paused ? "ON" : "OFF") << std::endl;
         break;

      case 'B':
         draw_bins = !draw_bins;
         std::cout << "Toggled: draw bins -> " << (draw_bins ? "ON" : "OFF") << std::endl;
         break;

      case 'G':
         draw_glyphs = !draw_glyphs;
         std::cout << "Toggled: draw glyphs -> " << (draw_glyphs ? "ON" : "OFF") << std::endl;
         break;

      case 'R':
         // Reset some stuff
         recording = !recording;
         std::cout << "Recording: " << ((recording) ? "ON" : "OFF") << std::endl;
         break;

      case 'S':
         engine->useSixth();
         std::cout << "Using Sixth integration scheme" << std::endl;
         break;

      case 'L':
         engine->useLeapfrog();
         std::cout << "Using Leapfrog integration scheme" << std::endl;
         break;

      case 'M':
         paused = true;
         std::cout << "Switching to material editing" << std::endl;
         handleMaterialEditing();
         break;

      case 'I':
         std::cout << "Injecting " << nparticles << " with current material" << std::endl;
         engine->injectParticles(nparticles);
         break;

      case 'F':
         std::cout << "FLIP THA HOUSE!" << std::endl;
         engine->flip();
         break;

      case 'N':
         std::cout << "MAKE SOME NOISE!" << std::endl;
         engine->swapCoords();
         break;

      default:
         ;
   }

   glutPostRedisplay();
}

void mouseEntryCallback(int state)
{
   if(state == GLUT_ENTERED)
   {
      /*printf("------------------------------\n");
      printf("!! MOUSE ENTERED WINDOW %d !!\n",state);
      printf("------------------------------\n");*/
      mouse_in_window = state;
   } else if(state == GLUT_LEFT) {
      /*printf("-------1----------------------\n");
      printf("!! MOUSE LEFT WINDOW %d !!\n", state);
      printf("------------------------------\n");*/
      mouse_in_window = state;
   }

}

// Passive CB: Mouse moving, no clicks
void mousePassiveMotionCallback(int x, int y)
{
   y = window_height - y;
   mv_x = x - prev_mx;
   mv_y = y - prev_my;
   prev_mx = x;
   prev_my = y;
}

// Active CB: Mouse moving while clicked
void mouseMotionCallback(int x, int y)
{
   if( !(x > window_width || x < 0 || y > window_height || y < 0) )
   {
      y = window_height - y;
      float x_s = x / static_cast<float>(window_width);
      float y_s = y / static_cast<float>(window_height);

      engine->injectParticle(x_s, y_s);
   }

   // Skip if painting off-screen
   glutPostRedisplay();

}

void mouseCallback(int button, int state, int x, int y)
{
   if( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
   {
      if( !(x > window_width || x < 0 || y > window_height || y < 0) )
      {
         y = window_height - y;
         float x_s = x / static_cast<float>(window_width);
         float y_s = y / static_cast<float>(window_height);

         engine->injectParticle(x_s, y_s);
      }

  }

   glutPostRedisplay();
}

// MAIN

int main(int argc, char** argv)
{
   lux::CmdLineFind clf(argc, argv);

   lux::Vector2d llc(-1, -1), urc(1, 1);

   window_width = clf.find("-w", 1000, "Window Width");
   window_height = clf.find("-h", 1000, "Window Height");

   int grid_width = clf.find("-Nx", 300, "Grid Width");
   int grid_height = clf.find("-Ny", 300, "Grid Height");

   timestep = clf.find("-dt", 0.01f, "Timestep");
   nparticles = clf.find("-npart", 50, "Particles to inject");
   frame_id = clf.find("-frame", "sph2d", "Frame id for recording");
   frame_store = clf.find("-fstore", "./frames", "Where to output frames");
   frame_max = clf.find("-fmax", 9999, "Max frame count");

   // Engine Parameters
   float roi = clf.find("-roi", 0.1f, "Radius of influence");

   float baseD = 1.0 * 10.0 / (M_PI * roi * roi);
   baseD = clf.find("-density", baseD, "Base density");

   float viscosity = clf.find("-viscosity", 1.0f, "Viscosity");
   float vE = clf.find("-veps", 0.01f, "Viscosity Epsilon");
   float pG = clf.find("-pgamma", 3.0f, "Pressure gamma");
   float pS = clf.find("-pscale", 10.0f, "Pressure scale");
   float cor = clf.find("-cor", 1.0f, "Coefficient of Restitution");
   float cof = clf.find("-cof", 1.0f, "Coefficient of Friction");

   std::string image_file = clf.find("-image", "", "Image for color");

   clf.usage("--help");
   clf.printFinds();

   std::cout << std::endl
      << "========================================" << std::endl
      << "Keyboard shortcuts: " << std::endl
      << std::endl
      << "Esc    :\t Exit" << std::endl
      << "r      :\t Reset particle storage (doesn't reset materials)" << std::endl
      << "P      :\t Start/Stop simulation" << std::endl
      << "R      :\t Toggle recording" << std::endl
      << "B      :\t Toggle display of occupancy volume" << std::endl
      << "G      :\t Toggle display of particle glyphs" << std::endl
      << "========================================" << std::endl
      << "1      :\t Use material color for particles" << std::endl
      << "2      :\t Use velocity color for particles" << std::endl
      << "4      :\t Toggle display of color grid" << std::endl
      << "5      :\t Toggle display of particles" << std::endl
      << "========================================" << std::endl
      << "L      :\t Use Leapfrog integration" << std::endl
      << "S      :\t Use Sixth integration" << std::endl
      << "M      :\t Open Material Editing menu in console" << std::endl
      << "I      :\t Inject -npart # of particles with current material " << std::endl
      << "========================================" << std::endl
      << "F      :\t Flip the tank" << std::endl
      << "N      :\t Rotate the tank" << std::endl
      << std::endl;

   if( image_file != "" )
   {
      readImage(image_file.c_str());
      window_width = img_width;
      window_height = img_height;

      std::cout <<
         "Image Width: " << img_width << std::endl <<
         "Image Height: " << img_height << std::endl;
   } else {
      size_t pixels = window_width * window_height * 3;
      master_image = new float[pixels];
      std::memset(master_image, 0, sizeof(float) * pixels);
   }
   window_x_position = window_y_position = 0;

   heatmapGradient.loadFromCSV(std::string("./colormaps/cubeyf1.csv"), 16);

   // Create and init SPHEngine
   engine = new sim::SPHEngine(urc, llc, timestep, roi, grid_width, grid_height);

   // Setup config for engine
   engine->initializeRandomDevice(0.0, 0.5);

   // Create some default materials
   engine->genMaterial(lux::Vector(0.0, 1.0, 0.5), 1.0, baseD, pS, pG, viscosity, vE);
   engine->_materials.back()->_name = "Regular";

   engine->genMaterial(lux::Vector(1.0, 0.0, 0.5), 1.0, baseD / 2.0, pS, pG, viscosity, vE);
   engine->_materials.back()->_name = "0.5x Density";

   engine->genMaterial(lux::Vector(0.0, 0.5, 1.0), 5.0, 2.0 * baseD, pS, pG, 2.0 * viscosity, vE);
   engine->_materials.back()->_name = "2x Density / 2x Viscosity";

   engine->genMaterial(lux::Vector(0.407, 0.271, 0.718), 10.0, 3.5 * baseD, pS, pG, 4.0 * viscosity, vE);
   engine->_materials.back()->_name = "Goop";

   engine->_current_material = 1;

   // Set Force parameters
   engine->getForce()->_gravity = lux::Vector2d(0, -9.8);

   // Collision
   engine->setCollisionCoeffs(cor, cof);

   engine->updateOccupancyVolume();

   // Create a cirle glyph for particles
   createGlyph(16);

   // Init GL and Extensions
   initGLUT(argc, argv);
   initGL();
   GLenum err = glewInit();
   if( GLEW_OK != err )
   {
      std::cerr << "Oops, it broke: " << glewGetString(err) << std::endl;
      exit(1);
   }

   loadTextures();

   /* Be sure to cleanup resources
    * Don't do this: threads get mixed up and try to double delete items
    * atexit(cleanUp);
    */

   glutMainLoop();

return 0;
}

// INITIALIZE

int initGLUT(int argc, char **argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
   glutInitWindowSize(window_width, window_height);
   glutInitWindowPosition(window_x_position, window_y_position);
   int handle = glutCreateWindow(argv[0]);

   glutDisplayFunc(draw);
   glutReshapeFunc(reshapeCallback);
   glutKeyboardFunc(keyboardCallback);
   glutMouseFunc(mouseCallback);
   glutPassiveMotionFunc(mousePassiveMotionCallback);
   glutMotionFunc(mouseMotionCallback);
   glutEntryFunc(mouseEntryCallback);
   glutIdleFunc(idleCallback);

   // Allow for us to cleanup after glut ends
   glutCloseFunc(cleanUp);

return handle;
}

void initGL()
{
   glShadeModel(GL_SMOOTH);

   glEnable(GL_TEXTURE_2D);

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   //glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
   //glEnable(GL_COLOR_MATERIAL);

   glClearColor(0, 0, 0, 0);                   // background color
   glClearDepth(1.0f);                         // 0 is near, 1 is far
   glDepthFunc(GL_LEQUAL);
}
