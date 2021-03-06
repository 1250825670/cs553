#include <stdio.h>
	// yes, I know stdio.h is not good C++, but I like the *printf( )
#include <stdlib.h>
#include <ctype.h>

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef WIN32
#include <windows.h>
#pragma warning(disable:4996)
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include "glut.h"
#include "glui.h"

//	This is a sample OpenGL / GLUT / GLUI program
//
//	The objective is to draw a 3d object and change the color of the axes
//		with radio buttons
//
//	The left mouse button allows rotation
//	The middle mouse button allows scaling
//	The glui window allows:
//		1. The 3d object to be transformed
//		2. The projection to be changed
//		3. The color of the axes to be changed
//		4. The axes to be turned on and off
//		5. The transformations to be reset
//		6. The program to quit
//
//	Author:			Joe Graphics
//
//  Latest update:	March 26, 2015

// constants:
//
// NOTE: There are a bunch of good reasons to use const variables instead
// of #define's.  However, Visual C++ does not allow a const variable
// to be used as an array size or as the case in a switch( ) statement.  So in
// the following, all constants are const variables except those which need to
// be array sizes or cases in switch( ) statements.  Those are #defines.
//


// title of these windows:

const char *WINDOWTITLE = { "cs553 project8 -- Li Li" };
const char *GLUITITLE   = { "User Interface Window" };


// what the glui package defines as true and false:

const int GLUITRUE  = { true  };
const int GLUIFALSE = { false };


// the escape key:

#define ESCAPE		0x1b


// initial window size:

const int INIT_WINDOW_SIZE = { 600 };


// size of the box:

const float BOXSIZE = { 2.f };



// multiplication factors for input interaction:
//  (these are known from previous experience)

const float ANGFACT = { 1. };
const float SCLFACT = { 0.005f };


// able to use the left mouse for either rotation or scaling,
// in case have only a 2-button mouse:

enum LeftButton
{
	ROTATE,
	SCALE
};


// minimum allowable scale factor:

const float MINSCALE = { 0.05f };


// active mouse buttons (or them together):

const int LEFT   = { 4 };
const int MIDDLE = { 2 };
const int RIGHT  = { 1 };


// which projection:

enum Projections
{
	ORTHO,
	PERSP
};


// which button:

enum ButtonVals
{
	RESET,
	QUIT
};


// window background color (rgba):

const float BACKCOLOR[ ] = { 0., 0., 0., 0. };


// line width for the axes:

const GLfloat AXES_WIDTH   = { 3. };


// the color numbers:
// this order must match the radio button order

enum Colors
{
	RED,
	YELLOW,
	GREEN,
	CYAN,
	BLUE,
	MAGENTA,
	WHITE,
	BLACK
};


// the color definitions:
// this order must match the radio button order

const GLfloat Colors[ ][3] = 
{
	{ 1., 0., 0. },		// red
	{ 1., 1., 0. },		// yellow
	{ 0., 1., 0. },		// green
	{ 0., 1., 1. },		// cyan
	{ 0., 0., 1. },		// blue
	{ 1., 0., 1. },		// magenta
	{ 1., 1., 1. },		// white
	{ 0., 0., 0. },		// black
};


// fog parameters:

const GLfloat FOGCOLOR[4] = { .0, .0, .0, 1. };
const GLenum  FOGMODE     = { GL_LINEAR };
const GLfloat FOGDENSITY  = { 0.30f };
const GLfloat FOGSTART    = { 1.5 };
const GLfloat FOGEND      = { 4. };


// non-constant global variables:

int	ActiveButton;			// current button that is down
GLuint	AxesList;			// list to hold the axes
int	AxesOn;					// != 0 means to draw the axes
int	DebugOn;				// != 0 means to print debugging info
int	DepthCueOn;				// != 0 means to use intensity depth cueing
GLUI *	Glui;				// instance of glui window
int	GluiWindow;				// the glut id for the glui window
int	LeftButton;				// either ROTATE or SCALE
GLuint	BoxList;			// object display list
int	MainWindow;				// window id for main graphics window
GLfloat	RotMatrix[4][4];	// set by glui rotation widget
float	Scale, Scale2;		// scaling factors
int	WhichColor;				// index into Colors[ ]
int	WhichProjection;		// ORTHO or PERSP
int	Xmouse, Ymouse;			// mouse values
float	Xrot, Yrot;			// rotation angles in degrees
float	TransXYZ[3];		// set by glui translation widgets


inline float SQR(float x)
{
	return x * x;
}

struct centers
{
	float xc, yc, zc;       // center location
	float a;                // amplitude
} Centers[] =
{
	{ 1.00f, 0.00f, 0.00f, 90.00f },
	{ -1.00f, 0.30f, 0.00f, 120.00f },
	{ 0.00f, 1.00f, 0.00f, 120.00f },
	{ 0.00f, 0.40f, 1.00f, 170.00f },
};

//functions and variables for proj8
#define CLOUD 1
#define EXAGGERATE 2
#define TEXTURE 3

GLUI_Checkbox *Cloud;
int Cloud_live = 1;

GLUI_RadioGroup *Texture;
GLUI_RadioButton *Notex;
GLUI_RadioButton *Replace;
GLUI_RadioButton *Modulate;

float Exagg;
GLUI_Spinner *ExaggHeight;

bool Paused = 0;

GLuint TerrainTex;
GLuint TerrainTex2;

#define NUMLNGS 201
#define NUMLATS 105

const float LNGMIN = { -289.6f };
const float LNGMAX = { 289.6f };
const float LATMIN = { -197.5f };
const float LATMAX = { 211.2f };
const float HGTMIN = { 0.0f };
const float HGTMAX = { 3.429f };

float Time = 0.;
float x=20., y = 20.;
float delta = 0.;

#define Pi 3.1416
int CloudWindow;
int FunctionWindow;

struct LngLatHgt
{
	float lng, lat, hgt;
}Points[NUMLATS][NUMLNGS];

#define Rad 100.0
#define rad 60.0

extern unsigned char *BmpToTexture(char *filename, int *width, int *height);

void InitGraphics2();
void InitGraphics3();

void Display2();
void Display3();


void Checkbox(int id)
{
	glutSetWindow(MainWindow);
	glutPostRedisplay();
	glutSetWindow(CloudWindow);
	glutPostRedisplay();
	glutSetWindow(FunctionWindow);
	glutPostRedisplay();
}

void Spinner(int id)
{
	glutSetWindow(MainWindow);
	glutPostRedisplay();
}

void NormalColor(int z, int x, float HeightExag);

// function prototypes:

void	Animate( );
void	Buttons( int );
void	Display( );
void	DoRasterString( float, float, float, char * );
void	DoStrokeString( float, float, float, float, char * );
float	ElapsedSeconds( );
void	InitGlui( );
void	InitGraphics( );
void	InitLists( );
void	Keyboard( unsigned char, int, int );
void	MouseButton( int, int, int, int );
void	MouseMotion( int, int );
void	Reset( );
void	Resize( int, int );
void	Visibility( int );

void	Arrow( float [3], float [3] );
void	Cross( float [3], float [3], float [3] );
float	Dot( float [3], float [3] );
float	Unit( float [3], float [3] );
void	Axes( float );
void	HsvRgb( float[3], float [3] );

// main program:

int
main( int argc, char *argv[ ] )
{
	// turn on the glut package:
	// (do this before checking argc and argv since it might
	// pull some command line arguments out)

	glutInit( &argc, argv );


	// setup all the graphics stuff:

	InitGraphics( );
	InitGraphics2();
	InitGraphics3();

	// create the display structures that will not change:

	InitLists( );


	// init all the global variables used by Display( ):
	// this will also post a redisplay
	// it is important to call this before InitGlui( )
	// so that the variables that glui will control are correct
	// when each glui widget is created

	Reset( );


	// setup all the user interface stuff:

	InitGlui( );


	// draw the scene once and wait for some interaction:
	// (this will never return)

	glutMainLoop( );


	// this is here to make the compiler happy:

	return 0;
}


// this is where one would put code that is to be called
// everytime the glut main loop has nothing to do
//
// this is typically where animation parameters are set
//
// do not call Display( ) from here -- let glutMainLoop( ) do it

void
Animate( )
{
	// put animation stuff in here -- change some global variables
	// for Display( ) to find:
	// make sure every system it is 10 second animate
	int ms = glutGet(GLUT_ELAPSED_TIME);
	ms = ms % (10 * 1000);
	Time = (float)ms / (float)(10 * 1000 - 1);
	delta = Time;

	// force a call to Display( ) next time it is convenient:

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
	glutSetWindow(CloudWindow);
	glutPostRedisplay();
	glutSetWindow(FunctionWindow);
	glutPostRedisplay();
}



// glui buttons callback:

void
Buttons( int id )
{
	switch( id )
	{
		case RESET:
			Reset( );
			Glui->sync_live( );
			glutSetWindow( MainWindow );
			glutPostRedisplay( );
			break;

		case QUIT:
			// gracefully close the glui window:
			// gracefully close out the graphics:
			// gracefully close the graphics window:
			// gracefully exit the program:

			Glui->close( );
			glutSetWindow( MainWindow );
			glFinish( );
			glutDestroyWindow( MainWindow );
			exit( 0 );
			break;

		default:
			fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}

}



// draw the complete scene:

void
Display()
{
	if (DebugOn != 0)
	{
		fprintf(stderr, "Display\n");
	}


	// set which window we want to do the graphics into:

	glutSetWindow(MainWindow);


	// erase the background:

	glDrawBuffer(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);


	// specify shading to be flat:

	glShadeModel(GL_SMOOTH);


	// set the viewport to a square centered in the window:

	GLsizei vx = glutGet(GLUT_WINDOW_WIDTH);
	GLsizei vy = glutGet(GLUT_WINDOW_HEIGHT);
	GLsizei v = vx < vy ? vx : vy;			// minimum dimension
	GLint xl = (vx - v) / 2;
	GLint yb = (vy - v) / 2;
	glViewport(xl, yb, v, v);


	// set the viewing volume:
	// remember that the Z clipping  values are actually
	// given as DISTANCES IN FRONT OF THE EYE
	// USE gluOrtho2D( ) IF YOU ARE DOING 2D !

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (WhichProjection == ORTHO)
		glOrtho(-3., 3., -3., 3., 0.1, 1000.);
	else
		gluPerspective(90., 1., 0.1, 1000.);


	// place the objects into the scene:

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();


	// set the eye position, look-at position, and up-vector:
	// IF DOING 2D, REMOVE THIS -- OTHERWISE ALL YOUR 2D WILL DISAPPEAR !

	gluLookAt(0., 300., 0., 0., 0., 0., 0., 0., -1.);


	// translate the objects in the scene:
	// note the minus sign on the z value
	// this is to make the appearance of the glui z translate
	// widget more intuitively match the translate behavior
	// DO NOT TRANSLATE IN Z IF YOU ARE DOING 2D !

	glTranslatef((GLfloat)TransXYZ[0], (GLfloat)TransXYZ[1], -(GLfloat)TransXYZ[2]);


	// rotate the scene:
	// DO NOT ROTATE (EXCEPT ABOUT Z) IF YOU ARE DOING 2D !

	glRotatef((GLfloat)Yrot, 0., 1., 0.);
	glRotatef((GLfloat)Xrot, 1., 0., 0.);
	glMultMatrixf((const GLfloat *)RotMatrix);


	// uniformly scale the scene:

	glScalef((GLfloat)Scale, (GLfloat)Scale, (GLfloat)Scale);
	GLfloat scale2 = 1. + Scale2;		// because glui translation starts at 0.
	if (scale2 < MINSCALE)
		scale2 = MINSCALE;
	glScalef((GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2);


	// set the fog parameters:
	// DON'T NEED THIS IF DOING 2D !

	if (DepthCueOn != 0)
	{
		glFogi(GL_FOG_MODE, FOGMODE);
		glFogfv(GL_FOG_COLOR, FOGCOLOR);
		glFogf(GL_FOG_DENSITY, FOGDENSITY);
		glFogf(GL_FOG_START, FOGSTART);
		glFogf(GL_FOG_END, FOGEND);
		glEnable(GL_FOG);
	}
	else
	{
		glDisable(GL_FOG);
	}


	// possibly draw the axes:

	if (AxesOn != 0)
	{
		glColor3fv(&Colors[WhichColor][0]);
		glCallList(AxesList);
	}


	// set the color of the object:

	glColor3fv(Colors[WhichColor]);


	// draw the current object:

	float dlng = (LNGMAX - LNGMIN) / (NUMLNGS - 1.);
	float dlat = (LATMAX - LATMIN) / (NUMLATS - 1.);
	float s, t;
	int z, x;
	float HeightExag = ExaggHeight->get_float_val();

	
	if (Texture->get_int_val() != 0)
	{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, TerrainTex);
		if (Texture->get_int_val() == 1)
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		else
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glBegin(GL_TRIANGLES);
		for (z = 0; z < NUMLATS; z++)
		{
			for (x = 0; x < NUMLNGS; x++)
			{

				if (z + 1 < NUMLATS && x + 1 < NUMLNGS)
				{
					s = (Points[z][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					glTexCoord2f(s, t);
					NormalColor(z, x, HeightExag);
					glVertex3f(Points[z][x].lng, Points[z][x].hgt * HeightExag, -Points[z][x].lat);
					s = (Points[z][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					NormalColor(z, x + 1, HeightExag);
					glTexCoord2f(s, t);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					s = (Points[z + 1][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z + 1][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					NormalColor(z + 1, x, HeightExag);
					glTexCoord2f(s, t);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
					s = (Points[z][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					glTexCoord2f(s, t);
					NormalColor(z, x + 1, HeightExag);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					s = (Points[z + 1][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z + 1][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					NormalColor(z + 1, x + 1, HeightExag);
					glTexCoord2f(s, t);
					glVertex3f(Points[z + 1][x + 1].lng, Points[z + 1][x + 1].hgt * HeightExag, -Points[z + 1][x + 1].lat);
					s = (Points[z + 1][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
					t = (Points[z + 1][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
					NormalColor(z + 1, x, HeightExag);
					glTexCoord2f(s, t);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
				}
			}
		}
		glEnd();
		glDisable(GL_TEXTURE_2D);
	}
	else
	{
		glBegin(GL_TRIANGLES);
		for (z = 0; z < NUMLATS; z++)
		{
			for (x = 0; x < NUMLNGS; x++)
			{
				if (z + 1 < NUMLATS && x + 1 < NUMLNGS)
				{
					NormalColor(z, x, HeightExag);
					glVertex3f(Points[z][x].lng, Points[z][x].hgt * HeightExag, -Points[z][x].lat);
					NormalColor(z, x + 1, HeightExag);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					NormalColor(z + 1, x, HeightExag);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
					NormalColor(z, x + 1, HeightExag);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					NormalColor(z + 1, x + 1, HeightExag);
					glVertex3f(Points[z + 1][x + 1].lng, Points[z + 1][x + 1].hgt * HeightExag, -Points[z + 1][x + 1].lat);
					NormalColor(z + 1, x, HeightExag);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
				}
			}
		}
		glEnd();
	}
	if (Cloud->get_int_val() != 0)
	{
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		
			float a = Rad * sin(2. * Pi*delta);
			float b = Rad * cos(2. * Pi*delta);

		for (float localrad = rad; localrad > 0.; localrad -= 3.)
		{
			glBegin(GL_TRIANGLE_STRIP);
			glColor4f(1., 0., 1., .4);
			int n = 12;
		
			for (int i = 0; i < n; i++)
				glVertex3f(a + localrad*sin(2. * Pi / n*i), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n*i)));
			glVertex3f(a + localrad*sin(2. * Pi / n*0), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n*0)));
			glVertex3f(a + localrad*sin(2. * Pi / n * 1), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n * 1)));
			glEnd();
		}
		
		
		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
	}

	// draw some gratuitous text that just rotates on top of the scene:

	glDisable( GL_DEPTH_TEST );
	glColor3f( 0., 1., 1. );
	//DoRasterString( 0., 1., 0., "Text That Moves" );


	// draw some gratuitous text that is fixed on the screen:
	//
	// the projection matrix is reset to define a scene whose
	// world coordinate system goes from 0-100 in each axis
	//
	// this is called "percent units", and is just a convenience
	//
	// the modelview matrix is reset to identity as we don't
	// want to transform these coordinates

	glDisable( GL_DEPTH_TEST );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	gluOrtho2D( 0., 100.,     0., 100. );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
	glColor3f( 1., 1., 1. );
	//DoRasterString( 5., 5., 0., "Text That Doesn't" );


	// swap the double-buffered framebuffers:

	glutSwapBuffers( );


	// be sure the graphics buffer has been sent:
	// note: be sure to use glFlush( ) here, not glFinish( ) !

	glFlush( );
}


void Display2()
{
	glutSetWindow(CloudWindow);
	glDrawBuffer(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	GLsizei vx = glutGet(GLUT_WINDOW_WIDTH);
	GLsizei vy = glutGet(GLUT_WINDOW_HEIGHT);
	GLsizei v = vx < vy ? vx : vy;			// minimum dimension
	GLint xl = (vx - v) / 2;
	GLint yb = (vy - v) / 2;
	float a = Rad * sin(2 * Pi*delta);
	float b = Rad * cos(2 * Pi*delta);
	float HeightExag = ExaggHeight->get_float_val();
	glViewport(xl, yb, v, v);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (WhichProjection == ORTHO)
		glOrtho(-3., 3., -3., 3., 0.1, 1000.);
	else
		gluPerspective(90., 1., 0.1, 1000.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if (Cloud->get_int_val() == 0)
		gluLookAt(0., 300., 0., 0., 0., 0., 0., 0., -1.);
	else
		gluLookAt(a, HGTMAX * HeightExag + 20, -b, a + rad*sin(2 * Pi*delta + Pi / 2), 0., -(b + rad * cos(2 * Pi*delta + Pi / 2)), 0., 1., 0.);
	float dlng = (LNGMAX - LNGMIN) / (NUMLNGS - 1.);
	float dlat = (LATMAX - LATMIN) / (NUMLATS - 1.);
	float s, t;
	int z, x;

	if (Texture->get_int_val() != 0)
	{
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, TerrainTex2);
			if (Texture->get_int_val() == 1)
				glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
			else
				glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glBegin(GL_TRIANGLES);
			for (z = 0; z < NUMLATS; z++)
			{
				for (x = 0; x < NUMLNGS; x++)
				{

					if (z + 1 < NUMLATS && x + 1 < NUMLNGS)
					{
						s = (Points[z][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						glTexCoord2f(s, t);
						NormalColor(z, x, HeightExag);
						glVertex3f(Points[z][x].lng, Points[z][x].hgt * HeightExag, -Points[z][x].lat);
						s = (Points[z][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						NormalColor(z, x + 1, HeightExag);
						glTexCoord2f(s, t);
						glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
						s = (Points[z + 1][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z + 1][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						NormalColor(z + 1, x, HeightExag);
						glTexCoord2f(s, t);
						glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
						s = (Points[z][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						glTexCoord2f(s, t);
						NormalColor(z, x + 1, HeightExag);
						glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
						s = (Points[z + 1][x + 1].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z + 1][x + 1].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						NormalColor(z + 1, x + 1, HeightExag);
						glTexCoord2f(s, t);
						glVertex3f(Points[z + 1][x + 1].lng, Points[z + 1][x + 1].hgt * HeightExag, -Points[z + 1][x + 1].lat);
						s = (Points[z + 1][x].lng - LNGMIN) / (LNGMAX - LNGMIN);
						t = (Points[z + 1][x].lat - LATMIN) / (LATMAX - LATMIN) * 0.718f;
						NormalColor(z + 1, x, HeightExag);
						glTexCoord2f(s, t);
						glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
					}
				}
			}
			glEnd();
			glDisable(GL_TEXTURE_2D);
		}
	else
	{
		glBegin(GL_TRIANGLES);
		for (z = 0; z < NUMLATS; z++)
		{
			for (x = 0; x < NUMLNGS; x++)
			{
				if (z + 1 < NUMLATS && x + 1 < NUMLNGS)
				{
					NormalColor(z, x, HeightExag);
					glVertex3f(Points[z][x].lng, Points[z][x].hgt * HeightExag, -Points[z][x].lat);
					NormalColor(z, x + 1, HeightExag);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					NormalColor(z + 1, x, HeightExag);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
					NormalColor(z, x + 1, HeightExag);
					glVertex3f(Points[z][x + 1].lng, Points[z][x + 1].hgt * HeightExag, -Points[z][x + 1].lat);
					NormalColor(z + 1, x + 1, HeightExag);
					glVertex3f(Points[z + 1][x + 1].lng, Points[z + 1][x + 1].hgt * HeightExag, -Points[z + 1][x + 1].lat);
					NormalColor(z + 1, x, HeightExag);
					glVertex3f(Points[z + 1][x].lng, Points[z + 1][x].hgt * HeightExag, -Points[z + 1][x].lat);
				}
			}
		}
		glEnd();
	}
	if (Cloud->get_int_val() != 0)
	{
		glEnable(GL_BLEND);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


		for (float localrad = rad; localrad > 0.; localrad -= 3.)
		{
			glBegin(GL_TRIANGLE_STRIP);
			glColor4f(1., 0., 1., .4);
			int n = 12;

			for (int i = 0; i < n; i++)
				glVertex3f(a + localrad*sin(2. * Pi / n*i), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n*i)));
			glVertex3f(a + localrad*sin(2. * Pi / n * 0), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n * 0)));
			glVertex3f(a + localrad*sin(2. * Pi / n * 1), HGTMAX*HeightExag, -(b + localrad*cos(2. * Pi / n * 1)));
			glEnd();
		}

		glDepthMask(GL_TRUE);
		glDisable(GL_BLEND);
	}
	glutSwapBuffers();
	glFlush();
}

void Display3()
{
	glutSetWindow(FunctionWindow);
	glDrawBuffer(GL_BACK);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_FLAT);
	GLsizei vx = glutGet(GLUT_WINDOW_WIDTH);
	GLsizei vy = glutGet(GLUT_WINDOW_HEIGHT);
	GLsizei v = vx < vy ? vx : vy;			// minimum dimension
	GLint xl = (vx - v) / 2;
	GLint yb = (vy - v) / 2;
	float a = Rad * sin(2 * Pi*delta);
	float b = Rad * cos(2 * Pi*delta);
	glViewport(xl, yb, v, v);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0., 2., 0., 2.);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glBegin(GL_LINES);
	glColor3f(1., 1., 0.);
	glLineWidth(2.);
	glVertex2f(0.2, 1.2);		//x
	glVertex2f(0.2, 1.8);
	glVertex2f(0.2, 1.2);
	glVertex2f(1.8, 1.2);
	glVertex2f(0.2, 1.5);
	glVertex2f(1.8, 1.5);
	glVertex2f(0.2, 0.2);		//y
	glVertex2f(0.2, 0.8);
	glVertex2f(0.2, 0.2);
	glVertex2f(1.8, 0.2);
	glVertex2f(0.2, 0.5);
	glVertex2f(1.8, 0.5);
	glEnd();
	glBegin(GL_QUADS);
	glColor3f(1., 0., 0.);
	glVertex2f(0.2 + 1.6*delta, 1.5 + (0.6*a / (2. * Rad)));
	glVertex2f(0.2 + 1.6*delta + 0.05, 1.5 + (0.6*a / (2. * Rad)));
	glVertex2f(0.2 + 1.6*delta + 0.05, 1.5 + (0.6*a / (2. * Rad)) - 0.05);
	glVertex2f(0.2 + 1.6*delta, 1.5 + (0.6*a / (2. * Rad)) - 0.05);
	glVertex2f(0.2 + 1.6*delta, 0.5 + (0.6*b / (2. * Rad)));
	glVertex2f(0.2 + 1.6*delta + 0.05, 0.5 + (0.6*b / (2. * Rad)));
	glVertex2f(0.2 + 1.6*delta + 0.05, 0.5 + (0.6*b / (2. * Rad)) - 0.05);
	glVertex2f(0.2 + 1.6*delta, 0.5 + (0.6*b / (2. * Rad)) - 0.05);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glLineWidth(2.);
	glColor3f(0., 1., 0.);
	for (float x = 0.; x <= 2 * Pi; x += 0.01)
		glVertex2f(0.2 + (1.6*x / (2. * Pi)), 1.5 + (0.6*sin(x) / 2.));
	glEnd();
	glBegin(GL_LINE_STRIP);
	glLineWidth(2.);
	glColor3f(0., 1., 0.);
	for (float x = 0.; x <= 2 * Pi; x += 0.01)
		glVertex2f(0.2 + (1.6*x / (2. * Pi)), 0.5 + (0.6*cos(x) / 2.));
	glEnd();
	glDisable(GL_DEPTH_TEST);
	glColor3f(1., 1., 1.);
	DoRasterString(1., 1., 0., "a");
	DoRasterString(1., 0.05, 0., "b");
	glutSwapBuffers();
	glFlush();
}

// use glut to display a string of characters using a raster font:

void
DoRasterString( float x, float y, float z, char *s )
{
	char c;			// one character to print

	glRasterPos3f( (GLfloat)x, (GLfloat)y, (GLfloat)z );
	for( ; ( c = *s ) != '\0'; s++ )
	{
		glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, c );
	}
}


// use glut to display a string of characters using a stroke font:

void
DoStrokeString( float x, float y, float z, float ht, char *s )
{
	char c;			// one character to print

	glPushMatrix( );
		glTranslatef( (GLfloat)x, (GLfloat)y, (GLfloat)z );
		float sf = ht / ( 119.05 + 33.33 );
		glScalef( (GLfloat)sf, (GLfloat)sf, (GLfloat)sf );
		for( ; ( c = *s ) != '\0'; s++ )
		{
			glutStrokeCharacter( GLUT_STROKE_ROMAN, c );
		}
	glPopMatrix( );
}


// return the number of seconds since the start of the program:

float
ElapsedSeconds( )
{
	// get # of milliseconds since the start of the program:

	int ms = glutGet( GLUT_ELAPSED_TIME );

	// convert it to seconds:

	return (float)ms / 1000.;
}


// initialize the glui window:

void
InitGlui( )
{
	GLUI_Panel *panel;
	GLUI_RadioGroup *group;
	GLUI_Rotation *rot;
	GLUI_Translation *trans, *scale;


	// setup the glui window:

	glutInitWindowPosition( INIT_WINDOW_SIZE + 50, 0 );
	Glui = GLUI_Master.create_glui( (char *) GLUITITLE );


	Glui->add_statictext( (char *) GLUITITLE );
	Glui->add_separator( );

	Glui->add_checkbox( "Axes", &AxesOn );

	Glui->add_checkbox( "Perspective", &WhichProjection );

	Glui->add_checkbox( "Intensity Depth Cue", &DepthCueOn );


	//panel = Glui->add_panel(  "Axes Color" );
	//	group = Glui->add_radiogroup_to_panel( panel, &WhichColor );
	//		Glui->add_radiobutton_to_group( group, "Red" );
	//		Glui->add_radiobutton_to_group( group, "Yellow" );
	//		Glui->add_radiobutton_to_group( group, "Green" );
	//		Glui->add_radiobutton_to_group( group, "Cyan" );
	//		Glui->add_radiobutton_to_group( group, "Blue" );
	//		Glui->add_radiobutton_to_group( group, "Magenta" );
	//		Glui->add_radiobutton_to_group( group, "White" );
	//		Glui->add_radiobutton_to_group( group, "Black" );

	//add GLuis for proj8
	
	panel = Glui->add_panel("Cloud and Exaggeration");
	
		ExaggHeight = Glui->add_spinner_to_panel(panel, "Exaggerate Height", GLUI_SPINNER_FLOAT, 0, EXAGGERATE, (GLUI_Update_CB)Spinner);
		ExaggHeight->set_float_limits(1.0f, 20.0f, GLUI_LIMIT_CLAMP);
		ExaggHeight->set_speed(1.0f);
		
		Glui->add_column_to_panel(panel, GLUIFALSE);
		Cloud = Glui->add_checkbox_to_panel(panel, "Cloud", &Cloud_live, CLOUD, (GLUI_Update_CB)Checkbox);

		
	panel = Glui->add_panel("Texture");
		Texture = Glui->add_radiogroup_to_panel(panel, 0, TEXTURE, (GLUI_Update_CB)Checkbox);
		Glui->add_column_to_panel(panel, GLUIFALSE);
		Notex = Glui->add_radiobutton_to_group(Texture, "No Texture");
		Glui->add_column_to_panel(panel, GLUIFALSE);
		Replace = Glui->add_radiobutton_to_group(Texture, "GL_REPLACE");
		Glui->add_column_to_panel(panel, GLUIFALSE);
		Modulate = Glui->add_radiobutton_to_group(Texture, "GL_MODULATE");
		Texture->set_int_val(1);
	
	
	panel = Glui->add_panel( "Object Transformation" );

		rot = Glui->add_rotation_to_panel( panel, "Rotation", (float *) RotMatrix );

		// allow the object to be spun via the glui rotation widget:

		rot->set_spin( 1.0 );


		Glui->add_column_to_panel( panel, GLUIFALSE );
		scale = Glui->add_translation_to_panel( panel, "Scale",  GLUI_TRANSLATION_Y , &Scale2 );
		scale->set_speed( 0.005f );

		Glui->add_column_to_panel( panel, GLUIFALSE );
		trans = Glui->add_translation_to_panel( panel, "Trans XY", GLUI_TRANSLATION_XY, &TransXYZ[0] );
		trans->set_speed( 0.05f );

		Glui->add_column_to_panel( panel, GLUIFALSE );
		trans = Glui->add_translation_to_panel( panel, "Trans Z",  GLUI_TRANSLATION_Z , &TransXYZ[2] );
		trans->set_speed( 0.05f );

	Glui->add_checkbox( "Debug", &DebugOn );


	panel = Glui->add_panel( "", GLUIFALSE );

	Glui->add_button_to_panel( panel, "Reset", RESET, (GLUI_Update_CB) Buttons );

	Glui->add_column_to_panel( panel, GLUIFALSE );

	Glui->add_button_to_panel( panel, "Quit", QUIT, (GLUI_Update_CB) Buttons );


	// tell glui what graphics window it needs to post a redisplay to:

	Glui->set_main_gfx_window( MainWindow );


	// set the graphics window's idle function if needed:

	GLUI_Master.set_glutIdleFunc( Animate );
}


// initialize the glut and OpenGL libraries:
//	also setup display lists and callback functions

void
InitGraphics( )
{
	

	// setup the display mode:
	// ( *must* be done before call to glutCreateWindow( ) )
	// ask for color, double-buffering, and z-buffering:

	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );


	// set the initial window configuration:

	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize( INIT_WINDOW_SIZE, INIT_WINDOW_SIZE );


	// open the window and set its title:

	MainWindow = glutCreateWindow( WINDOWTITLE );
	glutSetWindowTitle( WINDOWTITLE );


	// setup the clear values:

	glClearColor( BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3] );


	// setup the callback functions:

	// DisplayFunc -- redraw the window
	// ReshapeFunc -- handle the user resizing the window
	// KeyboardFunc -- handle a keyboard input
	// MouseFunc -- handle the mouse button going down or up
	// MotionFunc -- handle the mouse moving with a button down
	// PassiveMotionFunc -- handle the mouse moving with a button up
	// VisibilityFunc -- handle a change in window visibility
	// EntryFunc	-- handle the cursor entering or leaving the window
	// SpecialFunc -- handle special keys on the keyboard
	// SpaceballMotionFunc -- handle spaceball translation
	// SpaceballRotateFunc -- handle spaceball rotation
	// SpaceballButtonFunc -- handle spaceball button hits
	// ButtonBoxFunc -- handle button box hits
	// DialsFunc -- handle dial rotations
	// TabletMotionFunc -- handle digitizing tablet motion
	// TabletButtonFunc -- handle digitizing tablet button hits
	// MenuStateFunc -- declare when a pop-up menu is in use
	// TimerFunc -- trigger something to happen a certain time from now
	// IdleFunc -- what to do when nothing else is going on

	glutSetWindow( MainWindow );
	glutDisplayFunc( Display );
	glutReshapeFunc( Resize );
	glutKeyboardFunc( Keyboard );
	glutMouseFunc( MouseButton );
	glutMotionFunc( MouseMotion );
	glutPassiveMotionFunc( NULL );
	glutVisibilityFunc( Visibility );
	glutEntryFunc( NULL );
	glutSpecialFunc( NULL );
	glutSpaceballMotionFunc( NULL );
	glutSpaceballRotateFunc( NULL );
	glutSpaceballButtonFunc( NULL );
	glutButtonBoxFunc( NULL );
	glutDialsFunc( NULL );
	glutTabletMotionFunc( NULL );
	glutTabletButtonFunc( NULL );
	glutMenuStateFunc( NULL );
	glutTimerFunc( 0, NULL, 0 );

	// DO NOT SET THE GLUT IDLE FUNCTION HERE !!
	// glutIdleFunc( NULL );
	// let glui take care of it in InitGlui( )
	
	//init proj8 data
	unsigned char *Texture;
	int Width, Height;
	
	Texture = BmpToTexture("Oregon2048.bmp", &Width, &Height);
	if (Texture == NULL)
		fprintf(stderr, "Texture is NULL\n");
	if (Width != 2048 && Height != 2048)
		fprintf(stderr, "Texture size not right\n");
	
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &TerrainTex);
	glBindTexture(GL_TEXTURE_2D, TerrainTex);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, Width, Height, 0, GL_RGB, GL_UNSIGNED_BYTE, Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	FILE *fp = fopen("or.hgt", "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Cannot open .hgt file\n");
		exit(1);
	}
	float dlng = (LNGMAX - LNGMIN) / (NUMLNGS - 1);
	float dlat = (LATMAX - LATMIN) / (NUMLATS - 1);
	float lat, lng;
	int z, x;
	for (z = 0, lat = LATMAX; z < NUMLATS; z++, lat -= dlat)
	{
		for (x = 0, lng = LNGMIN; x < NUMLNGS; x++, lng += dlng)
		{
			Points[z][x].lat = lat;
			Points[z][x].lng = lng;
			fscanf(fp, "%f", &Points[z][x].hgt);
		}
	}
}


void InitGraphics2()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(100, 500);
	glutInitWindowSize(200, 200);
	CloudWindow = glutCreateWindow(WINDOWTITLE);
	glutSetWindowTitle(WINDOWTITLE);
	glClearColor(BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3]);
	glutSetWindow(CloudWindow);
	glutKeyboardFunc(Keyboard);
	glutDisplayFunc(Display2);

	unsigned char *Texture;
	int Width, Height;
	Texture = BmpToTexture("Oregon2048.bmp", &Width, &Height);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glGenTextures(1, &TerrainTex2);
	glBindTexture(GL_TEXTURE_2D, TerrainTex2);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, Width, Height, 0, GL_RGB, GL_UNSIGNED_BYTE, Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

void InitGraphics3()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(400, 500);
	glutInitWindowSize(200, 200);
	FunctionWindow = glutCreateWindow(WINDOWTITLE);
	glutSetWindowTitle(WINDOWTITLE);
	glClearColor(BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3]);
	glutSetWindow(FunctionWindow);
	glutKeyboardFunc(Keyboard);
	glutDisplayFunc(Display3);
}

// initialize the display lists that will not change:
// (a display list is a way to store opengl commands in
//  memory so that they can be played back efficiently at a later time
//  with a call to glCallList( )

void
InitLists( )
{
	
	// create the axes:

	AxesList = glGenLists( 1 );
	glNewList( AxesList, GL_COMPILE );
		glLineWidth( AXES_WIDTH );
			Axes( 1.5 );
		glLineWidth( 1. );
	glEndList( );
}


// the keyboard callback:

void
Keyboard( unsigned char c, int x, int y )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Keyboard: '%c' (0x%0x)\n", c, c );

	switch( c )
	{
		case 'o':
		case 'O':
			WhichProjection = !WhichProjection;
			break;

		case 'p':
		case 'P':
			Paused = !Paused;
			if (Paused)
				GLUI_Master.set_glutIdleFunc(NULL);
			else
				GLUI_Master.set_glutIdleFunc(Animate);
			break;

		case 'q':
		case 'Q':
		case ESCAPE:
			Buttons( QUIT );	// will not return here
			break;			// happy compiler

		case 'r':
		case 'R':
			LeftButton = ROTATE;
			break;

		case 's':
		case 'S':
			LeftButton = SCALE;
			break;

		default:
			fprintf( stderr, "Don't know what to do with keyboard hit: '%c' (0x%0x)\n", c, c );
	}


	// synchronize the GLUI display with the variables:

	Glui->sync_live( );


	// force a call to Display( ):

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// called when the mouse button transitions down or up:

void
MouseButton( int button, int state, int x, int y )
{
	int b = 0;			// LEFT, MIDDLE, or RIGHT

	if( DebugOn != 0 )
		fprintf( stderr, "MouseButton: %d, %d, %d, %d\n", button, state, x, y );

	
	// get the proper button bit mask:

	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			b = LEFT;		break;

		case GLUT_MIDDLE_BUTTON:
			b = MIDDLE;		break;

		case GLUT_RIGHT_BUTTON:
			b = RIGHT;		break;

		default:
			b = 0;
			fprintf( stderr, "Unknown mouse button: %d\n", button );
	}


	// button down sets the bit, up clears the bit:

	if( state == GLUT_DOWN )
	{
		Xmouse = x;
		Ymouse = y;
		ActiveButton |= b;		// set the proper bit
	}
	else
	{
		ActiveButton &= ~b;		// clear the proper bit
	}
}


// called when the mouse moves while a button is down:

void
MouseMotion( int x, int y )
{
	if( DebugOn != 0 )
		fprintf( stderr, "MouseMotion: %d, %d\n", x, y );


	int dx = x - Xmouse;		// change in mouse coords
	int dy = y - Ymouse;

	if( ( ActiveButton & LEFT ) != 0 )
	{
		switch( LeftButton )
		{
			case ROTATE:
				Xrot += ( ANGFACT*dy );
				Yrot += ( ANGFACT*dx );
				break;

			case SCALE:
				Scale += SCLFACT * (float) ( dx - dy );
				if( Scale < MINSCALE )
					Scale = MINSCALE;
				break;
		}
	}


	if( ( ActiveButton & MIDDLE ) != 0 )
	{
		Scale += SCLFACT * (float) ( dx - dy );

		// keep object from turning inside-out or disappearing:

		if( Scale < MINSCALE )
			Scale = MINSCALE;
	}

	Xmouse = x;			// new current position
	Ymouse = y;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}

void NormalColor(int z, int x, float HeightExag)
{
	float gradLong[3], gradLat[3], norm[3];
	float LightingHeight = 1.0f;

	if (x - 1 < 0)
	{
		gradLong[0] = Points[z][x + 1].lng - Points[z][x].lng;
		gradLong[1] = HeightExag*LightingHeight*(Points[z][x + 1].hgt - Points[z][x].hgt);
	}
	else if (x + 1 > NUMLNGS)
	{
		gradLong[0] = Points[z][x].lng - Points[z][x - 1].lng;
		gradLong[1] = HeightExag*LightingHeight*(Points[z][x].hgt - Points[z][x - 1].hgt);
	}
	else
	{
		gradLong[0] = Points[z][x + 1].lng - Points[z][x - 1].lng;
		gradLong[1] = HeightExag*LightingHeight*(Points[z][x + 1].hgt - Points[z][x - 1].hgt);
	}
	gradLong[2] = 0.;
	if (z - 1 < 0)
	{
		gradLat[1] = HeightExag*LightingHeight*(Points[z + 1][x].hgt - Points[z][x].hgt);
		gradLat[2] = Points[z + 1][x].lat - Points[z][x].lat;
	}
	else if (z + 1 > NUMLATS)
	{
		gradLat[1] = HeightExag*LightingHeight*(Points[z][x].hgt - Points[z - 1][x].hgt);
		gradLat[2] = Points[z][x].lat - Points[z - 1][x].lat;
	}
	else
	{
		gradLat[1] = HeightExag*LightingHeight*(Points[z + 1][x].hgt - Points[z - 1][x].hgt);
		gradLat[2] = Points[z + 1][x].lat - Points[z - 1][x].lat;
	}
	gradLat[0] = 0;
	Cross(gradLong, gradLat, norm);
	Unit(norm, norm);
	float intern = fabs(norm[1]);
	glColor3f(intern, intern, intern);
}

// reset the transformations and the colors:
// this only sets the global variables --
// the glut main loop is responsible for redrawing the scene

void
Reset( )
{
	
	ActiveButton = 0;
	AxesOn = GLUITRUE;
	DebugOn = GLUIFALSE;
	DepthCueOn = GLUIFALSE;
	LeftButton = ROTATE;
	Scale  = 1.0;
	Scale2 = 0.0;		// because we add 1. to it in Display( )
	WhichColor = WHITE;
	WhichProjection = PERSP;
	Xrot = Yrot = 0.;
	TransXYZ[0] = TransXYZ[1] = TransXYZ[2] = 0.;

	                  RotMatrix[0][1] = RotMatrix[0][2] = RotMatrix[0][3] = 0.;
	RotMatrix[1][0]                   = RotMatrix[1][2] = RotMatrix[1][3] = 0.;
	RotMatrix[2][0] = RotMatrix[2][1]                   = RotMatrix[2][3] = 0.;
	RotMatrix[3][0] = RotMatrix[3][1] = RotMatrix[3][3]                   = 0.;
	RotMatrix[0][0] = RotMatrix[1][1] = RotMatrix[2][2] = RotMatrix[3][3] = 1.;

}


// called when user resizes the window:

void
Resize( int width, int height )
{
	if( DebugOn != 0 )
		fprintf( stderr, "ReSize: %d, %d\n", width, height );

	// don't really need to do anything since window size is
	// checked each time in Display( ):

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}


// handle a change to the window's visibility:

void
Visibility ( int state )
{
	if( DebugOn != 0 )
		fprintf( stderr, "Visibility: %d\n", state );

	if( state == GLUT_VISIBLE )
	{
		glutSetWindow( MainWindow );
		glutPostRedisplay( );
	}
	else
	{
		// could optimize by keeping track of the fact
		// that the window is not visible and avoid
		// animating or redrawing it ...
	}
}




///////////////////////////////////////   HANDY UTILITIES:  //////////////////////////

// size of wings as fraction of length:

#define WINGS	0.10


// axes:

#define X	1
#define Y	2
#define Z	3


// x, y, z, axes:

static float axx[3] = { 1., 0., 0. };
static float ayy[3] = { 0., 1., 0. };
static float azz[3] = { 0., 0., 1. };


void
Arrow( float tail[3], float head[3] )
{
	float u[3], v[3], w[3];		// arrow coordinate system

	// set w direction in u-v-w coordinate system:

	w[0] = head[0] - tail[0];
	w[1] = head[1] - tail[1];
	w[2] = head[2] - tail[2];


	// determine major direction:

	int axis = X;
	float mag = fabs( w[0] );
	if(  fabs( w[1] )  > mag  )
	{
		axis = Y;
		mag = fabs( w[1] );
	}
	if(  fabs( w[2] )  > mag  )
	{
		axis = Z;
		mag = fabs( w[2] );
	}


	// set size of wings and turn w into a Unit vector:

	float d = WINGS * Unit( w, w );


	// draw the shaft of the arrow:

	glBegin( GL_LINE_STRIP );
		glVertex3fv( tail );
		glVertex3fv( head );
	glEnd( );

	// draw two sets of wings in the non-major directions:

	float x, y, z;

	if( axis != X )
	{
		Cross( w, axx, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
	}


	if( axis != Y )
	{
		Cross( w, ayy, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
	}



	if( axis != Z )
	{
		Cross( w, azz, v );
		(void) Unit( v, v );
		Cross( v, w, u  );
		x = head[0] + d * ( u[0] - w[0] );
		y = head[1] + d * ( u[1] - w[1] );
		z = head[2] + d * ( u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
		x = head[0] + d * ( -u[0] - w[0] );
		y = head[1] + d * ( -u[1] - w[1] );
		z = head[2] + d * ( -u[2] - w[2] );
		glBegin( GL_LINE_STRIP );
			glVertex3fv( head );
			glVertex3f( x, y, z );
		glEnd( );
	}
}



float
Dot( float v1[3], float v2[3] )
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}



void
Cross( float v1[3], float v2[3], float vout[3] )
{
	float tmp[3];

	tmp[0] = v1[1]*v2[2] - v2[1]*v1[2];
	tmp[1] = v2[0]*v1[2] - v1[0]*v2[2];
	tmp[2] = v1[0]*v2[1] - v2[0]*v1[1];

	vout[0] = tmp[0];
	vout[1] = tmp[1];
	vout[2] = tmp[2];
}



float
Unit( float vin[3], float vout[3] )
{
	float dist = vin[0]*vin[0] + vin[1]*vin[1] + vin[2]*vin[2];

	if( dist > 0.0 )
	{
		dist = sqrt( dist );
		vout[0] = vin[0] / dist;
		vout[1] = vin[1] / dist;
		vout[2] = vin[2] / dist;
	}
	else
	{
		vout[0] = vin[0];
		vout[1] = vin[1];
		vout[2] = vin[2];
	}

	return dist;
}



// the stroke characters 'X' 'Y' 'Z' :

static float xx[ ] = {
		0.f, 1.f, 0.f, 1.f
	      };

static float xy[ ] = {
		-.5f, .5f, .5f, -.5f
	      };

static int xorder[ ] = {
		1, 2, -3, 4
		};


static float yx[ ] = {
		0.f, 0.f, -.5f, .5f
	      };

static float yy[ ] = {
		0.f, .6f, 1.f, 1.f
	      };

static int yorder[ ] = {
		1, 2, 3, -2, 4
		};


static float zx[ ] = {
		1.f, 0.f, 1.f, 0.f, .25f, .75f
	      };

static float zy[ ] = {
		.5f, .5f, -.5f, -.5f, 0.f, 0.f
	      };

static int zorder[ ] = {
		1, 2, 3, 4, -5, 6
		};


// fraction of the length to use as height of the characters:

const float LENFRAC = 0.10f;


// fraction of length to use as start location of the characters:

const float BASEFRAC = 1.10f;


//	Draw a set of 3D axes:
//	(length is the axis length in world coordinates)

void
Axes( float length )
{
	glBegin( GL_LINE_STRIP );
		glVertex3f( length, 0., 0. );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., length, 0. );
	glEnd( );
	glBegin( GL_LINE_STRIP );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., 0., length );
	glEnd( );

	float fact = LENFRAC * length;
	float base = BASEFRAC * length;

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 4; i++ )
		{
			int j = xorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( base + fact*xx[j], fact*xy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 5; i++ )
		{
			int j = yorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( fact*yx[j], base + fact*yy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 6; i++ )
		{
			int j = zorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( 0.0, fact*zy[j], base + fact*zx[j] );
		}
	glEnd( );

}


// function to convert HSV to RGB
//
// Reference:  Foley, van Dam, Feiner, Hughes,
//		"Computer Graphics Principles and Practices,"

void
HsvRgb( float hsv[3], float rgb[3] )
{
	float r, g, b;			// red, green, blue

	// guarantee valid input:

	float h = hsv[0] / 60.;
	while( h >= 6. )	h -= 6.;
	while( h <  0. ) 	h += 6.;

	float s = hsv[1];
	if( s < 0. )
		s = 0.;
	if( s > 1. )
		s = 1.;

	float v = hsv[2];
	if( v < 0. )
		v = 0.;
	if( v > 1. )
		v = 1.;


	// if sat==0, then is a gray:

	if( s == 0.0 )
	{
		rgb[0] = rgb[1] = rgb[2] = v;
		return;
	}


	// get an rgb from the hue itself:
	
	float i = floor( h );
	float f = h - i;
	float p = v * ( 1. - s );
	float q = v * ( 1. - s*f );
	float t = v * ( 1. - ( s * (1.-f) ) );

	switch( (int) i )
	{
		case 0:
			r = v;	g = t;	b = p;
			break;
	
		case 1:
			r = q;	g = v;	b = p;
			break;
	
		case 2:
			r = p;	g = v;	b = t;
			break;
	
		case 3:
			r = p;	g = q;	b = v;
			break;
	
		case 4:
			r = t;	g = p;	b = v;
			break;
	
		case 5:
			r = v;	g = p;	b = q;
			break;
	}


	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}
