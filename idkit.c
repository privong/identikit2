/*
 * idkit.c: fit Identakit N-body model to interacting galaxy XYV data.
 * Copyright (C) 2012 by Joshua Edward Barnes, Honolulu, Hawaii.  This
 * is free software, distributed under the terms of GNU General Public
 * License.  See <http://www.gnu.org/licenses/> for details.
 */

#include "stdinc.h"
#include "mathfns.h"
#include "getparam.h"
#include "vectmath.h"
#include "filestruct.h"
#include "phatbody.h"

#if defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

//  ______________________________________________________
//  Command line parameters and their description strings.

string defv[] = {		";Fit Identakit model to XYVz data.",
				";Drag mouse to adjust view; command keys:",
				";? <sp> <del> n p x y z r R l m 1 2 c s o v",
				";t Q q D L j J V M w W d C Z S b f $ <esc>",
    "jhat=???",			";Spin direction and weight data.",
				";Used to select disk particles to display.",
    "sims=???",			";Patterns for identakit simulation files.",
				";Embedded printf directive uses current",
				";value of step to generate filenames.",
    "zrot=0.0",			";Rotate simulation about initial z axis",
    "isim=0",			";Initially display this simulation",
    "step=0x200",		";Initially display this step number",
    "incr=0x020",		";Increment/decrement for step number",
    "ncenter=2048",		";Number of bodies used to track centers",
    "maxview=32768",		";Maximum number of bodies to display",
    "panelsize=400",		";Size of panels; entire window is 2x2",
    "xyimage=",			";Image for XY projection",
    "vyimage=",			";Image for VY projection",
    "xvimage=",			";Image for XV projection",
    "datacube=",		";XYV data in particle form",
    "initview=",		";Initial view file",
    "saveview=idk%03d.dat",	";Output pattern for view files",
    "saveimage=idk%03d.ppm",	";Output pattern for screen dump",
    "spindepth=4",		";Tessellation depth for spin ball",
    "viewdepth=2",		";Tessellation depth for view search",
    "xwidth=250.0",		";Image X,Y width in kpc",
    "vwidth=500.0",		";Image V width in km/sec",
    "script=",			";String of keyboard commands",
    "VERSION=2.2" VARIANT,	";Joshua Barnes  11 June 2012.",
				";This program is free software; you can",
				";redistribute it under certain conditions.",
				";This program comes with NO WARRANTY.",
    NULL,
};

//  __________________________________________________________
//  idkstate: structure for current and past identikit states.

#define NUMPARS  10			// num. of parameter-pair variables

typedef struct {
  real par[NUMPARS][2];			// values of variable parameters
  matrix pmat;				// preliminary transformation mat.
  int step;				// current value of step number
  int isim;				// index into list of simulations
} idkstate;

//  _____________________________________________________
//  XYROT, ..., JVIEW: indexes into parameter-pair array.

#define XYROT  0			// rotation about X and Y axes
#define SCALE  1			// X,Y scale and rotation about Z
#define DISK1  2			// incl. and peri. arg. for disk 1
#define DISK2  3			// incl. and peri. arg. for disk 2
#define XYOFF  4			// offsets in X and Y directions
#define VELOS  5			// V offset and scale factor
#define ITRAN  6			// image contrast and zeropoint
#define BOXXY  7			// selection box X and Y coords.
#define BOXVR  8			// selection box V coord and size
#define JVIEW  9			// rotate view of spin ball

//  ________________________________________________
//  box: structure for phase-space selection region.

typedef struct {
  vector XYbox;				// screen X, Y coords. of box
  real Vbox;				// screen V coord. of box
  real Rbox;				// size of box
  bool velflag;				// select V as well as X and Y
  int disk;				// which disk box selects
} box;

#define MAXBOX   32			// max. number of phase-space boxes

//  ___________________________________________________________________
//  facet: structure for triangular face of solid approximating a ball.

typedef struct {
  vector vert1, vert2, vert3;		// coordinates of corners
  vector midp;				// unit vector through midpoint
  real spin[2], prod[2], view[2];	// values on ball, indexed by disk
} facet;

#define NUMICOS  20			// num. of faces on icosahedron

//  ______________________________________________________________________
//  XZSHOW, ..., VWSHOW: content of lower-right panel, selected by lrshow.

#define XZSHOW 0			// show (X,Z) projection
#define BXSHOW (MAXBOX)			// show spin ball (1:BXSHOW)
#define PRSHOW (MAXBOX+1)		// show product ball
#define VWSHOW (MAXBOX+2)		// show view ball

//  ______________________________________
//  Global variables and input parameters.

idkstate idk;				// current identikit state
bool parset[NUMPARS];			// flag changed parameter-pairs

#define MAXHIST  1024			// max. number of saved states
idkstate old[MAXHIST];			// past identikit history
int ihist = 0;				// index into history array
int nhist = 0;				// num. of items in history array

box boxtab[MAXBOX+1];			// storage for phase-space boxes
int nbox = 0;				// box indices run from 1 to nbox

facet *ball;				// tessellation of icosahedron
int nfacet;				// number of facets in tessellation

bodyptr btab = NULL;			// identikit array of bodies
int nbody = 0;				// number of bodies in btab[]
real tnow;				// current time from input file

bodyptr vtab = NULL;			// array of visible bodies
int maxview, nview1, nview2;		// max size, number from each disk

bodyptr ctab = NULL;			// data cube body array
int ncube = 0;				// number of bodies in data cube

string *sims;				// patterns for input filenames
string *jhat;				// names of jhat files for above
string *zrot;				// values for initial z rot'n
int nsims;				// number of names in these lists

GLsizei pansize, imgsize = -1;		// size of panels and bkgnd images
GLubyte *xyimage = NULL, *xyimage1;	// image behind XY projection
GLubyte *vyimage = NULL, *vyimage1;	// image behind VY projection
GLubyte *xvimage = NULL, *xvimage1;	// image behind XV projection

int butbind[] = { XYROT,DISK1,DISK2 };	// map mouse buttons to parameters
int actpars = -1, actlast = XYROT;	// current and last param adjusted
int xlast, ylast;			// previous position of mouse

int ncenter;				// bodies used to find centers

real disktol = 0.01;			// sets number of bodies visible

vector cpos1, cpos2;			// screen positions of centers
real cvel1[2], cvel2[2];		// screen vel. ranges of centers
bool cposflag = FALSE;			// TRUE if cpos vectors are set
bool cvelflag = FALSE;			// TRUE if cvel ranges are set
int centflag = 0;			// if > 0, lock centers at cpos

int lrshow = XZSHOW;			// what to show in LR panel
int parflag = 1;			// display parameter values
bool cubeflag = FALSE;			// if TRUE, show data cube
int zoomflag = 0;			// if >0, zoom on one panel

int ndirect;				// number of directions to scan

char *script;				// command keys to be executed
bool animate = FALSE;			// TRUE if animation in progress

char msgbuf[256], titlebuf[256];	// buffers for messages and title

//  __________________________________________________________
//  Function and procedure prototypes, in order of appearance.

void display(void);
void seldisks(void);
void lockcent(void);
void mapdisks(void);
void adjustimage(GLubyte *img1, GLubyte *img0, real scale, real bias);
void displayimage(GLubyte *image);
void displaycube(int xind, int yind, real xfield, real vfield);
void displaydisks(int xind, int yind);
void displayball(void);
bool loadball(matrix vmat, real xbox, real ybox, real vbox,
	      real rbox, real hbox, int disk);
facet *showball(matrix tmat, real *rgb1, real *rgb2, int ind);
void showspins(matrix jmat, vector spin1, vector spin2);
void displaypars(void);
void keyboard(unsigned char key, int x, int y);
void printhelp(void);
void storehist(void);
void scanboxes(void);
void scanviews(int arg);
void findmax(void);
void printscale(void);
void mousebut(int but, int state, int x, int y);
void mousemove(int x, int y);
void special(int key, int x, int y);
void getjhat(string jhat, string zrot, bool firstcall);
void getdata(string spat, string zrot, bool firstcall);
void getcube(string cube);
void getview(void);
void putview(void);
void getimage(GLubyte **image, GLubyte **imcpy, string ifile);
string inputline(stream istr);
void putscreen(void);
void clrpars(bool all);
void getpars(matrix vmat, real *xfield, vector spin1, vector spin2,
	     vector xyoff, real *vfield, real *voff, real *scale, real *bias,
	     vector XYbox, real *Vbox, real *Rbox, matrix jmat);
void setpars(real *xyzrot, real *xfield, vector spin1, vector spin2,
	     vector xyoff, vector XYbox, real *Vbox, real *Rbox);
void showtext(string str, int x, int y);
void showcross(real x, real y, real size);
void showbox(real x, real y, real size);
void rotmatrix(matrix rmat, real xrot, real yrot, real zrot);
string *padlist(string *old, int len, string name);
void findcenters(void);
bool centbbox(real xmin[2], real xmax[2], real ymin[2], real ymax[2],
	      GLubyte *img, int size);
void buildsphere(int lev);
void buildtriangle(vector v1, vector v2, vector v3, int lev);

//  __________________
//  Macro definitions.

#define rsinD(x)  (((x) ==  0 || (x) == 180) ? 0 : rsin((x) * (PI / 180)))
#define rcosD(x)  (((x) == 90 || (x) == 270) ? 0 : rcos((x) * (PI / 180)))
#define rasinD(x)  (rasin(x) * (180 / PI))
#define racosD(x)  (racos(x) * (180 / PI))
#define ratan2D(x,y)  ((x<0 ? 360 : 0) + ratan2(x,y) * (180 / PI))
#define PERDIST  4.0			// viewing distance for ball
#define PerTrans(x,z)  (1.125 * (x) / (1.0 - (z) / PERDIST))

//  _________________________________________________________________
//  main: read initial data, initialize display, and start main loop.

int main(int argc, string argv[])
{
  string bodytags[] = { PosTag, VelTag, AuxTag, AuxVecTag, NULL, };

  glutInit(&argc, argv);
  initparam(argv, defv);
  layout_body(bodytags, Precision, NDIM);

  sims = burststring(getparam("sims"), ",; ");
  nsims = xstrlen(sims, sizeof(string)) - 1;
  if (nsims < 1)
    error("%s: input file list list empty\n", getargv0());
  jhat = padlist(burststring(getparam("jhat"), ",; "), nsims, "jhat");
  zrot = padlist(burststring(getparam("zrot"), ",; "), nsims, "zrot");
  idk.isim = getiparam("isim");
  if (idk.isim < 0 || idk.isim >= nsims)
    error("%s: isim out of range\n", getargv0());
  idk.step = getiparam("step");
  getjhat(jhat[idk.isim], zrot[idk.isim], TRUE);
  getdata(sims[idk.isim], zrot[idk.isim], TRUE);
  ncenter = getiparam("ncenter");
  maxview = getiparam("maxview");
  vtab = (bodyptr) allocate(maxview * SizeofBody);
  pansize = getiparam("panelsize");
  if (! strnull(getparam("xyimage")))
    getimage(&xyimage, &xyimage1, getparam("xyimage"));
  if (! strnull(getparam("vyimage")))
    getimage(&vyimage, &vyimage1, getparam("vyimage"));
  if (! strnull(getparam("xvimage")))
    getimage(&xvimage, &xvimage1, getparam("xvimage"));
  findcenters();
  if (! strnull(getparam("datacube")))
    getcube(getparam("datacube"));
  clrpars(TRUE);
  SETMI(idk.pmat);
  buildsphere(getiparam("spindepth"));
  if (! strnull(getparam("initview")))
    getview();
  storehist();
  ndirect = MIN(NUMICOS * (1 << (2 * getiparam("viewdepth"))), nfacet);
  script = getparam("script");
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(2 * pansize, 2 * pansize); 
  glutCreateWindow(getargv0());
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glutDisplayFunc(display); 
  glutMouseFunc(mousebut);
  glutMotionFunc(mousemove);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(special);
#if defined(BOLD_GRAPHICS)
  glLineWidth((GLfloat) BOLD_GRAPHICS);
  glPointSize((GLfloat) BOLD_GRAPHICS);
#endif
  glutMainLoop();
  return (0);
}

//  ________________________________________________
//  display: main routine called to update graphics.

void display(void)
{
  int box, i;
  real xfield, vfield, scale, bias;

  if (parset[DISK1] || parset[DISK2])
    seldisks();
  if (centflag > 0)
    lockcent();
  box = (lrshow <= BXSHOW ? lrshow : 0);	// box == 0 has no effect
  getpars(NULL, &xfield, NULL, NULL, NULL, &vfield, NULL, &scale, &bias,
	  boxtab[box].XYbox, &boxtab[box].Vbox, &boxtab[box].Rbox, NULL);
  mapdisks();
  if (parset[ITRAN] && xyimage != NULL)
    adjustimage(xyimage1, xyimage, scale, bias);
  if (parset[ITRAN] && xvimage != NULL)
    adjustimage(xvimage1, xvimage, scale, bias);
  if (parset[ITRAN] && vyimage != NULL)
    adjustimage(vyimage1, vyimage, scale, bias);
  glClear(GL_COLOR_BUFFER_BIT);
  glViewport(0, 0, 2*pansize, 2*pansize);
  if (zoomflag == 0 || zoomflag == 1) {
    if (zoomflag == 0)
      glViewport(0, pansize, pansize, pansize);
    displayimage(xyimage1);
    displaycube(0, 1, xfield, vfield);
    displaydisks(0, 1);
  }
  if (zoomflag == 0 || zoomflag == 2) {
    if (zoomflag == 0)
      glViewport(0, 0, pansize, pansize);
    displayimage(xvimage1);
    displaycube(0, 3, xfield, vfield);
    displaydisks(0, 3);
  }
  if (zoomflag == 0 || zoomflag == 3) {
    if (zoomflag == 0)
      glViewport(pansize, pansize, pansize, pansize);
    displayimage(vyimage1);
    displaycube(3, 1, xfield, vfield);
    displaydisks(3, 1);
  }
  if (zoomflag == 0 || zoomflag == 4) {
    if (zoomflag == 0)
      glViewport(pansize, 0, pansize, pansize);
    if (lrshow == XZSHOW)
      displaydisks(0, 2);			// show XZ projection
    else
      displayball();				// show function on ball
  }
  displaypars();				// show pars on last panel
  glFlush();
  glutSwapBuffers();
  glutSetWindowTitle(titlebuf);
  for (i = 0; i < NUMPARS; i++)
    parset[i] = FALSE;
  if (*script != (long) NULL && !animate)	// if script is playing
    keyboard(*script++, 0, 0);			// handle keys one-by-one
}

//  _______________________________________________________________
//  seldisks: select bodies in current disks, and locate centroids.

void seldisks(void)
{
  vector spin1, spin2, cmpos, cmvel;
  matrix zmat;
  int block, nc;
  bodyptr vp, bp, cp;

  getpars(NULL, NULL, spin1, spin2, NULL, NULL, NULL,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  vp = NthBody(vtab, 2);
  nview1 = nview2 = 0;
  block = 0;					// step thru body blocks
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    if (AuxVec(bp)[0] == 0 && AuxVec(bp)[1] == 0 && AuxVec(bp)[2] == 0) {
      if (block % 2 == 0) {			// start of sphr block?
	block++;
	cp = NthBody(vtab, block == 1 ? 0 : 1);
	CLRV(cmpos);
	CLRV(cmvel);
	nc = 0;
      }
      if (nc < ncenter) {
	nc++;
	ADDV(cmpos, cmpos, Pos(bp));
	ADDV(cmvel, cmvel, Vel(bp));
	DIVVS(Pos(cp), cmpos, nc);
	DIVVS(Vel(cp), cmvel, nc);
      }
    } else {
      if (block % 2 == 1)			// start of disk block?
	block++;
      if (dotvp(AuxVec(bp), block==2 ? spin1 : spin2) > 1 - disktol*Aux(bp)) {
	if (vp >= NthBody(vtab, maxview)) {
	  sprintf(msgbuf, "#rtoo many bodies to view");
	  break;
	}
	SETV(Pos(vp), Pos(bp));
	SETV(Vel(vp), Vel(bp));
	vp = NextBody(vp);
	(* (block==2 ? &nview1 : &nview2))++;
      }
    }
}

//  __________________________________________________________________
//  lockcent: adjust scale, z rotation, offsets to keep centers fixed.

void lockcent(void)
{
  real dxref, dyref, xfield, dxmod, dymod, sfact, dxoff, dyoff;
  matrix vmat;
  vector xyoff;
  bodyptr bp, c1 = NthBody(vtab, 0), c2 = NthBody(vtab, 1);

  dxref = cpos2[0] - cpos1[0];			// work in screen coords.
  dyref = cpos2[1] - cpos1[1];
  getpars(vmat, &xfield, NULL, NULL, xyoff, NULL, NULL,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  for (bp = vtab; bp < NthBody(vtab, 2); bp = NextBody(bp)) {
    MULMV(AuxVec(bp), vmat, Pos(bp));
    ADDV(AuxVec(bp), AuxVec(bp), xyoff);
    DIVVS(AuxVec(bp), AuxVec(bp), xfield);
  }
  dxmod = AuxVec(c2)[0] - AuxVec(c1)[0];
  dymod = AuxVec(c2)[1] - AuxVec(c1)[1];
  sfact = rsqrt((dxmod*dxmod + dymod*dymod) / (dxref*dxref + dyref*dyref));
  /********** convert following to call to setpars() **********/
  idk.par[SCALE][0] += (2/PI) * (ratan2(dxref, dyref) - ratan2(dxmod, dymod));
  idk.par[SCALE][1] += rlog10(sfact);
  parset[SCALE] = TRUE;
  /********** convert previous to call to setpars() **********/
  getpars(vmat, &xfield, NULL, NULL, xyoff, NULL, NULL,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  for (bp = vtab; bp < NthBody(vtab, 2); bp = NextBody(bp)) {
    MULMV(AuxVec(bp), vmat, Pos(bp));
    ADDV(AuxVec(bp), AuxVec(bp), xyoff);
    DIVVS(AuxVec(bp), AuxVec(bp), xfield);
  }
  dxoff = (cpos1[0] + cpos2[0])/2 - (AuxVec(c1)[0] + AuxVec(c2)[0])/2;
  dyoff = (cpos1[1] + cpos2[1])/2 - (AuxVec(c1)[1] + AuxVec(c2)[1])/2;
  xyoff[0] += xfield * dxoff;
  xyoff[1] += xfield * dyoff;
  setpars(NULL, NULL, NULL, NULL, xyoff, NULL, NULL, NULL);
}

//  _________________________________________________________________
//  mapdisks: apply viewing transformation to sampled disk particles.

void mapdisks(void)
{
  matrix vmat;
  real xfield, vfield, voff;
  vector xyoff, tmpv;
  bodyptr bp;

  getpars(vmat, &xfield, NULL, NULL, xyoff, &vfield, &voff,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  for (bp = vtab; bp < NthBody(vtab, 2+nview1+nview2); bp = NextBody(bp)) {
    MULMV(AuxVec(bp), vmat, Pos(bp));
    ADDV(AuxVec(bp), AuxVec(bp), xyoff);
    DIVVS(AuxVec(bp), AuxVec(bp), xfield);
    MULMV(tmpv, vmat, Vel(bp));
    Aux(bp) = (tmpv[2] + voff) / vfield;
  }
}

//  ____________________________________________________
//  adjustimage: perform linear transformation on image.

void adjustimage(GLubyte *img1, GLubyte *img0, real scale, real bias)
{
  int i, imgval;

  for (i = 0; i < 3 * imgsize * imgsize; i++) {
    imgval = 256.0 * (bias + scale * img0[i] / 256.0);
    img1[i] = MIN(255, MAX(0, imgval));
  }
}

//  ________________________________________
//  displayimage: fill panel with RGB image.

void displayimage(GLubyte *image)
{
  int effsize = (zoomflag > 0 ? 2 * pansize : pansize);

  if (image != NULL) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, effsize, effsize, 0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glRasterPos2f((GLfloat) 0, (GLfloat) 0);
    glPixelZoom((effsize / (real) imgsize), -(effsize / (real) imgsize));
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glDrawPixels(imgsize, imgsize, GL_RGB, GL_UNSIGNED_BYTE, image);
    glPopMatrix();
  }
}

//  _________________________________________________________
//  displaycube: display data cube, highlighting current box.

void displaycube(int xind, int yind, real xfield, real vfield)
{
  bodyptr cp;
  real XYV[3];

  if (cubeflag && ctab != NULL) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.5, 0.5, -0.5, 0.5, -100.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glBegin(GL_POINTS);
    for (cp = ctab; cp < NthBody(ctab, ncube); cp = NextBody(cp)) {
      XYV[0] = Pos(cp)[0] / xfield;
      XYV[1] = Pos(cp)[1] / xfield;
      XYV[2] = Pos(cp)[2] / vfield;
      if ((1 <= lrshow && lrshow <= BXSHOW) &&	// lrshow is box index
	  ABS(XYV[0] - boxtab[lrshow].XYbox[0]) < boxtab[lrshow].Rbox &&
	  ABS(XYV[1] - boxtab[lrshow].XYbox[1]) < boxtab[lrshow].Rbox &&
	  (ABS(XYV[2] - boxtab[lrshow].Vbox) < boxtab[lrshow].Rbox ||
	   ! boxtab[lrshow].velflag))
	glColor3f(0.0, 1.0, 0.0);
      else
	glColor3f(0.0, 0.5, 0.0);
      glVertex2f(xind<2 ? XYV[xind] : XYV[2],
		 yind<2 ? XYV[yind] : XYV[2]);
    }
    glEnd(/* GL_POINTS */);
    glPopMatrix();
  }
}

//  _____________________________________________________________
//  displaydisks: display specified projection of disk particles.

#define SelCoord(bp,ind)  ((ind) < 3 ? AuxVec(bp)[ind] : Aux(bp))

void displaydisks(int xind, int yind)
{
  int j, i;
  bodyptr bp;
  vector XYbox;
  real Rbox, Vbox, Xbox, Ybox, boxlev;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-0.5, 0.5, -0.5, 0.5, -100.0, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glPushMatrix();
  glBegin(GL_POINTS);
  for (j = 0; j < 2 * MAX(nview1, nview2); j++)
    if ((j & 1) == 0 && j/2 < nview1) {
      glColor3f(0.0, 1.0, 1.0);
      i = j/2 + 2;
      glVertex2f(SelCoord(NthBody(vtab, i), xind),
		 SelCoord(NthBody(vtab, i), yind));
    } else if (j/2 < nview2) {
      glColor3f(1.0, 0.0, 1.0);
      i = j/2 + 2 + nview1;
      glVertex2f(SelCoord(NthBody(vtab, i), xind),
		 SelCoord(NthBody(vtab, i), yind));
    }
  glEnd(/* GL_POINTS */);
  if (centflag > 0 && xind == 0 && yind == 1) {
    glColor3f(0.0, 1.0, 0.0);
    glRectf(cpos1[0]-0.01, cpos1[1]-0.01, cpos1[0]+0.01, cpos1[1]+0.01);
    glRectf(cpos2[0]-0.01, cpos2[1]-0.01, cpos2[0]+0.01, cpos2[1]+0.01);
  }
  if (centflag >= 0 || (centflag == -1 && xind == 0 && yind == 1)) {
    glColor3f(0.0, 0.0, 1.0);
    for (bp = vtab; bp < NthBody(vtab, 2); bp = NextBody(bp))
      showcross(SelCoord(bp, xind), SelCoord(bp, yind), 0.05);
  }
  for (i = 1; i <= nbox; i++)
    if ((xind == 0 && yind == 1) ||
	(boxtab[i].velflag && (xind == 3 || yind == 3))) {
      Xbox = (xind < 3 ? boxtab[i].XYbox[xind] : boxtab[i].Vbox);
      Ybox = (yind < 3 ? boxtab[i].XYbox[yind] : boxtab[i].Vbox);
      boxlev = (i == lrshow ? 1.0 : 0.7);	// highlight current box
      glColor3f((boxtab[i].disk == 1 ? 0.0 : boxlev),
		(boxtab[i].disk == 2 ? 0.0 : boxlev),
		(boxtab[i].disk == 0 ? 0.0 : boxlev));
      showbox(Xbox, Ybox, 2 * boxtab[i].Rbox);
    }
  glPopMatrix();
}

//  _________________________________________________
//  displayball: display spin, product, or view ball.

void displayball(void)
{
  matrix vmat, jmat;
  real xfield, vfield, voff, xbox, ybox, vbox, rbox, hbox;
  real rgb1[3] = { 0.0, 0.8, 0.8 }, rgb2[3] = { 0.8, 0.0, 0.8 };
  vector spin1, spin2, xyoff, jvec, jorbit = { 0.0, 0.0, 1.0 };
  facet *fp;
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1.5, 1.5, -1.5, 1.5, -100.0, 100.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glPushMatrix();
  getpars(vmat, &xfield, spin1, spin2, xyoff, &vfield, &voff,
	  NULL, NULL, NULL, NULL, NULL, jmat);
  if (1 <= lrshow && lrshow <= BXSHOW) {	// lrshow is box index
    xbox = boxtab[lrshow].XYbox[0] * xfield - xyoff[0];
    ybox = boxtab[lrshow].XYbox[1] * xfield - xyoff[1];
    vbox = boxtab[lrshow].Vbox * vfield - voff;
    rbox = boxtab[lrshow].Rbox * xfield;
    hbox = (boxtab[lrshow].velflag ? boxtab[lrshow].Rbox : 100.0) * vfield;
    (void) loadball(vmat, xbox, ybox, vbox, rbox, hbox, boxtab[lrshow].disk);
  }
  if (lrshow != VWSHOW) {			// show spin or prod ball
    (void) showball(jmat, rgb1, rgb2, lrshow <= BXSHOW ? 0 : 1);
    showspins(jmat, spin1, spin2);
    MULMV(jvec, jmat, jorbit);			// rotate orbital (z) axis
    glColor3f(1.0, 1.0, 1.0);
    glRectf(PerTrans(jvec[0] - 0.02, jvec[2]),	// mark z direction
	    PerTrans(jvec[1] - 0.02, jvec[2]),
	    PerTrans(jvec[0] + 0.02, jvec[2]),
	    PerTrans(jvec[1] + 0.02, jvec[2]));
  } else {					// show view ball
    fp = showball(vmat, rgb1, rgb2, 2);		// get nearest facet
    sprintf(msgbuf, "#yscores = %6.2f #c%6.2f #m%6.2f",
	    fp->view[0] * fp->view[1] > 0 ?
	      rlog10(fp->view[0] * fp->view[1]) / 2 : -99.0,
	    fp->view[0] > 0 ? rlog10(fp->view[0]) : -99.0,
	    fp->view[1] > 0 ? rlog10(fp->view[1]) : -99.0);
    if (centflag >= 0) {
      glColor3f(1.0, 1.0, 1.0);
      showbox(0.0, 0.0, 0.04);
    }
  }
  glPopMatrix();
}

//  _____________________________________________________________
//  loadball: find particles in box and load them onto spin ball.

bool loadball(matrix vmat, real xbox, real ybox, real vbox,
	      real rbox, real hbox, int disk)
{
  real rhosmth, rhoface, dotsmth, dotcrit, q, max1, max2;
  facet *sp;
  int block, npart[2], icf;
  bodyptr bp;

  rhosmth = racos(1 - rsqrt(2.0) * disktol);
  rhoface = racos(ball[0].midp[2]);
  dotsmth = rcos(rhosmth);
  dotcrit = rcos(rhosmth + rhoface);
  for (sp = ball; sp < &ball[nfacet]; sp++)
    sp->spin[0] = sp->spin[1] = 0.0;
  block = 0;					// count particle blocks
  npart[0] = npart[1] = 0;
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    if (AuxVec(bp)[0] != 0 || AuxVec(bp)[1] != 0 || AuxVec(bp)[2] != 0) {
      block = (block % 2 == 1 ? block + 1 : block);
      if ((2*disk == block || disk == 0) &&
	  rabs(dotvp(vmat[0], Pos(bp)) - xbox) < rbox &&
	  rabs(dotvp(vmat[1], Pos(bp)) - ybox) < rbox &&
	  rabs(dotvp(vmat[2], Vel(bp)) - vbox) < hbox) {
	npart[block/2 - 1]++;			// count particles in box
	for (icf = 0; icf < NUMICOS; icf++) {	// loop over primary facets
	  sp = &ball[icf * (nfacet / NUMICOS)];	// point to middle facet
	  if (dotvp(sp->midp, AuxVec(bp)) > dotcrit)
	    for (sp = sp; sp < &ball[(icf + 1) * (nfacet / NUMICOS)]; sp++)
	      if (dotvp(sp->midp, AuxVec(bp)) > dotsmth) {
		q = racos(dotvp(sp->midp, AuxVec(bp))) / (0.5 * rhosmth);
#if defined(NO_WEIGHTING)
		sp->spin[block/2 - 1] +=
		  (q < 1 ? 1 - 1.5*rsqr(q) + 0.75*rqbe(q) : rqbe(2 - q) / 4);
#else
		sp->spin[block/2 - 1] += Aux(bp) *
		  (q < 1 ? 1 - 1.5*rsqr(q) + 0.75*rqbe(q) : rqbe(2 - q) / 4);
#endif
	      }
	}
      }
    } else
      block = (block % 2 == 0 ? block + 1 : block);
#if !defined(NO_NORMALIZE)
  max1 = max2 = 0;
  for (sp = ball; sp < &ball[nfacet]; sp++) {
    max1 = MAX(max1, sp->spin[0]);
    max2 = MAX(max2, sp->spin[1]);
  }
  for (sp = ball; sp < &ball[nfacet]; sp++) {
    sp->spin[0] = (max1 > 0 ? sp->spin[0] / max1 : 0);
    sp->spin[1] = (max2 > 0 ? sp->spin[1] / max2 : 0);
  }
#endif
  return (npart[0] + npart[1] > 0);
}

//  _______________________________________________
//  showball: display functions on surface of ball.

facet *showball(matrix tmat, real *rgb1, real *rgb2, int ind)
{
  facet *fp, *fnear;
  real *func, max1 = 0, max2 = 0, znear = 0, f, g;
  vector mp, vp1, vp2, vp3;

  for (fp = ball; fp < &ball[nfacet]; fp++) {
    func = (ind == 0 ? fp->spin : ind == 1 ? fp->prod : fp->view);
    max1 = MAX(max1, func[0]);
    max2 = MAX(max2, func[1]);
  }
  for (fp = ball; fp < &ball[nfacet]; fp++) {
    MULMV(mp, tmat, fp->midp);
    if (mp[2] > 1.0 / PERDIST) {		// is midp visible?
      func = (ind == 0 ? fp->spin : ind == 1 ? fp->prod : fp->view);
      f = rcbrt(max1 > 0 ? func[0] / max1 : 0);
      g = rcbrt(max2 > 0 ? func[1] / max2 : 0);
      glColor3f(rgb1[0] * f + rgb2[0] * g,
		rgb1[1] * f + rgb2[1] * g,
		rgb1[2] * f + rgb2[2] * g);
      MULMV(vp1, tmat, fp->vert1);
      MULMV(vp2, tmat, fp->vert2);
      MULMV(vp3, tmat, fp->vert3);
      glBegin(GL_TRIANGLES);
      glVertex2f(PerTrans(vp1[0], vp1[2]), PerTrans(vp1[1], vp1[2]));
      glVertex2f(PerTrans(vp2[0], vp2[2]), PerTrans(vp2[1], vp2[2]));
      glVertex2f(PerTrans(vp3[0], vp3[2]), PerTrans(vp3[1], vp3[2]));
      glEnd(/* GL_TRIANGLES */);
      glColor3f(0.25, 0.25, 0.25);
      glBegin(GL_POINTS);
      glVertex2f(PerTrans(mp[0], mp[2]), PerTrans(mp[1], mp[2]));
      glEnd(/* GL_POINTS */);
      fnear = (mp[2] > znear ? fp : fnear);
      znear = (mp[2] > znear ? mp[2] : znear);
    }
  }
  return (fnear);				// return nearest facet
}

//  ________________________________________________
//  showspins: plot spins of particles in each disk.

void showspins(matrix jmat, vector spin1, vector spin2)
{
  int block = 0;
  bodyptr bp;
  vector jvec;

  glBegin(GL_POINTS);
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    if (AuxVec(bp)[0] != 0 || AuxVec(bp)[1] != 0 || AuxVec(bp)[2] != 0) {
      block = (block % 2 == 1 ? block + 1 : block);
      if (dotvp(AuxVec(bp), block==2 ? spin1 : spin2) > 1 - disktol) {
	MULMV(jvec, jmat, AuxVec(bp));
	glColor3f(block==2 ? 0.0 : 0.5, block==2 ? 0.5 : 0.0, 0.5);
	glVertex2f(PerTrans(jvec[0], jvec[2]), PerTrans(jvec[1], jvec[2]));
      }	
    } else
      block = (block % 2 == 0 ? block + 1 : block);
  glEnd(/* GL_POINTS */);
}

//  ________________________________________________________
//  displaypars: display parameters describing current view.

void displaypars(void)
{
  int effsize = (zoomflag > 0 ? 2 * pansize : pansize);
  matrix vmat;
  real xfield, vfield, voff;
  vector spin1, spin2, xyoff;
  char buf[64];

  if (parflag > 0) {
    getpars(vmat, &xfield, spin1, spin2, xyoff, &vfield, &voff,
	    NULL, NULL, NULL, NULL, NULL, NULL);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, effsize, effsize, 0, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    sprintf(buf, "#wview: %.1f %.1f %.1f", 90 * idk.par[XYROT][1],
	    90 * idk.par[XYROT][0], 90 * idk.par[SCALE][0]);
    showtext(buf, 10, 15);
    sprintf(buf, "#wtime: %.2f", tnow);
    showtext(buf, effsize - (8 * strlen(buf) + 10), 15);
    glColor3f(1.0, 1.0, 0.3);
    showtext(msgbuf, 10, 30);
    if (parflag > 1) {
      sprintf(buf, "#wxyoff: %.3f %.3f", xyoff[0], xyoff[1]);
      showtext(buf, 10, effsize - 40);
      sprintf(buf, "#wvoff: %.3f", voff);
      showtext(buf, effsize - (8 * strlen(buf) + 10), effsize - 40);
    }
    sprintf(buf, "#cdisk1: %.1f %.1f (%.1f)",
	    racosD(spin1[2]), ratan2D(spin1[0], spin1[1]),
	    racosD(dotvp(spin1, vmat[2])));
    showtext(buf, 10, effsize - 25);
    sprintf(buf, "#wxfield: %.2f", xfield);
    showtext(buf, effsize - (8 * strlen(buf) + 10), effsize - 25);
    sprintf(buf, "#mdisk2: %.1f %.1f (%.1f)",
	    racosD(spin2[2]), ratan2D(spin2[0], spin2[1]),
	    racosD(dotvp(spin2, vmat[2])));
    showtext(buf, 10, effsize - 10);
    sprintf(buf, "#wvfield: %.2f", vfield);
    showtext(buf, effsize - (8 * strlen(buf) + 10), effsize - 10);
    glPopMatrix();
  }
}

//  _______________________________
//  keyboard: process command keys.

void keyboard(unsigned char key, int x, int y)
{
  int i;
  vector XYbox;

  switch (key) {
    case '?':
      printhelp();
      break;
    case 040:
    case 010:
    case 0177:
      storehist();
      idk.step += (key == 040 ? 1 : -1) * getiparam("incr");
      getdata(sims[idk.isim], zrot[idk.isim], FALSE);
      break;
    case 'n':
    case 'p':
      storehist();
      idk.isim = (key == 'n' ? (idk.isim < nsims - 1 ? idk.isim + 1 : 0) :
		               (idk.isim > 0 ? idk.isim - 1 : nsims - 1));
      getjhat(jhat[idk.isim], zrot[idk.isim], FALSE);
      getdata(sims[idk.isim], zrot[idk.isim], FALSE);
      break;
    case 'x':
    case 'y':
    case 'z':
      storehist();
      rotmatrix(idk.pmat, key == 'y' ? 90 : 0, key == 'x' ? 90 : 0, 0);
      break;
    case 'r':
    case 'R':
      storehist();
      clrpars(key == 'R');
      break;
    case 'l':
    case 'm':
      disktol = disktol * rsqrt(rsqrt(key == 'l' ? 0.5 : 2.0));
      seldisks();
      sprintf(msgbuf, "#ydisk tol: %.4f  nview: %d,%d",
	      disktol, nview1, nview2);
      break;
    case '1':
    case '2':
      if (1 <= lrshow && lrshow <= BXSHOW) {	// current box = lrshow
	boxtab[lrshow].disk = (key == '1' ? 1 : 2);
	sprintf(msgbuf, "#ycurrent box selects disk %c", key);
      } else {					// no current box?
	butbind[key == '1' ? 1 : 2] = actlast = (key == '1' ? DISK1 : DISK2);
	sprintf(msgbuf, "#y%s button adjusts disk %c",
		key == '1' ? "2nd (opt)" : "3rd (cmd)", key);
      }
      break;
    case 'c':
      centflag = (centflag == 1 ? -2 : centflag + 1);
      if (centflag > 0 && cposflag == FALSE) {
	SETV(cpos1, AuxVec(NthBody(vtab, 0)));	// use current positions
	SETV(cpos2, AuxVec(NthBody(vtab, 1)));
      }
      sprintf(msgbuf, "#ycenters %s", centflag > 0 ? "locked" : "unlocked");
      break;
    case 's':
      butbind[1] = actlast = SCALE;
      sprintf(msgbuf, "#y2rd (opt) button sets scale & rotation");
      break;
    case 'o':
      butbind[2] = actlast = XYOFF;
      sprintf(msgbuf, "#y3rd (cmd) button sets XY offsets");
      break;
    case 'v':
      butbind[1] = actlast = VELOS;
      sprintf(msgbuf, "#y2nd (opt) button sets velocities");
      break;
    case 't':
      butbind[2] = actlast = ITRAN;
      sprintf(msgbuf, "#y3rd (cmd) button sets contrast");
      break;
    case 'Q':
    case 'q':
      if (nbox < MAXBOX) {
	lrshow = ++nbox;			// make new box current
	boxtab[lrshow].velflag = (key == 'Q');
	boxtab[lrshow].disk = 0;
	getpars(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		XYbox, NULL, NULL, NULL);
	XYbox[0] += 3 / ((real) pansize);
	XYbox[1] -= 3 / ((real) pansize);
	setpars(NULL, NULL, NULL, NULL, NULL, XYbox, NULL, NULL);
	butbind[1] = BOXXY;
	butbind[2] = BOXVR;
	sprintf(msgbuf, "#y2nd and 3rd buttons set box %d", lrshow);
      } else
	sprintf(msgbuf, "#rtoo many boxes (max = %d)", MAXBOX);
      break;
    case 'D':
      if (nbox > 0 && 1 <= lrshow && lrshow <= BXSHOW) {
	for (i = lrshow; i < nbox; i++)
	  boxtab[i] = boxtab[i + 1];
	nbox--;
	lrshow = MIN(lrshow, nbox);
	setpars(NULL, NULL, NULL, NULL, NULL, boxtab[lrshow].XYbox,
		&boxtab[lrshow].Vbox, &boxtab[lrshow].Rbox);
	sprintf(msgbuf, nbox > 0 ? "#y%d boxes left; box %d is current" :
		"#yno boxes left", nbox, lrshow);
      } else
	sprintf(msgbuf, "#rno box to delete");
      break;
    case 'L':
      printf("\n%3s  %6s  %6s  %6s  %6s\n",
	     "box", "Xbox", "Ybox", "Vbox", "Rbox");
      for (i = 1; i <= nbox; i++)
	printf("%2d%1s  %6.3f  %6.3f  %6.3f  %6.3f\n", i,
	       i == lrshow ? "*" : " ", boxtab[i].XYbox[0],
	       boxtab[i].XYbox[1], boxtab[i].Vbox, boxtab[i].Rbox);
      break;
    case 'j':
      lrshow = (lrshow < nbox ? lrshow + 1 : lrshow == nbox ? PRSHOW :
		lrshow == PRSHOW ? VWSHOW : XZSHOW);
      butbind[0] = (lrshow == XZSHOW || lrshow == VWSHOW ? XYROT : JVIEW);
      if (1 <= lrshow && lrshow <= BXSHOW) {	// is there a current box?
	butbind[1] = BOXXY;
	butbind[2] = BOXVR;
	setpars(NULL, NULL, NULL, NULL, NULL, boxtab[lrshow].XYbox,
		&boxtab[lrshow].Vbox, &boxtab[lrshow].Rbox);
      }
      sprintf(msgbuf, "#y1st button adjusts %s",
	      lrshow == VWSHOW ? "view ball" :
	      lrshow == PRSHOW ? "product ball" :
	      lrshow == XZSHOW ? "viewing direction" : "spin ball");
      break;
    case 'J':
      if (nbox > 0)
	scanboxes();
      else
	sprintf(msgbuf, "#rno boxes to scan");
      break;
    case 'V':
      if (nbox > 0 && centflag > 0) {
	animate = TRUE;
	glutTimerFunc(10, scanviews, 0);
      } else
	sprintf(msgbuf,
		nbox == 0 ? "#rno boxes to scan" : "#rcenters not locked");
      break;
    case 'M':
      findmax();
      break;
    case 'w':
      putview();
      break;
    case 'W':
      putscreen();
      break;
    case 'd':
      parflag = (parflag + 1) % 3;
      break;
    case 'C':
      cubeflag = !cubeflag;
      break;
    case 'Z':
      zoomflag = (zoomflag + 1) % 5;
      break;
    case 'b':
    case 'f':
      ihist = (key == 'b' ? MAX(ihist - 1, 0) : MIN(ihist + 1, nhist - 1));
      idk = old[ihist % MAXHIST];
      getjhat(jhat[idk.isim], zrot[idk.isim], FALSE);
      getdata(sims[idk.isim], zrot[idk.isim], FALSE);
      sprintf(msgbuf, "#yrestoring view %d", ihist);
      break;
    case 'S':
      printscale();
      break;
    case 033:
    case '$':
      exit(0);
    default:
      sprintf(msgbuf, "#rkeystroke 0%o ignored", (int) key);
  }
  for (i = 0; i < NUMPARS; i++)			// force complete recalc.
    parset[i] = TRUE;
  glutPostRedisplay();				// always redisplay screen
}

//  _________________________________________________
//  printhelp: print message describing command keys.

void printhelp(void)
{
  eprintf("\n%s (v%s) keyboard commands:\n\n", getargv0(), getversion());
  eprintf("  ?\tdisplay keyboard commands (this message)\n");
  eprintf("  <sp>\tnext time-step in current simulation\n");
  eprintf("  <del>\tprevious time-step in current simulation\n");
  eprintf("  n\tnext pattern in simulation list\n");
  eprintf("  p\tprevious pattern in simulation list\n");
  eprintf("  x\tproject along X axis\n");
  eprintf("  y\tproject along Y axis\n");
  eprintf("  z\tproject along Z axis\n");
  eprintf("  r\treset view angles\n");
  eprintf("  R\treset all view parameters\n");
  eprintf("  l\tdisplay less particles\n");
  eprintf("  m\tdisplay more particles\n");
  eprintf("  1\t2nd mouse button adjusts disk 1 or\n");
  eprintf("   \trestricts current box to disk 1\n");
  eprintf("  2\t3rd mouse button adjusts disk 2 or\n");
  eprintf("   \trestricts current box to disk 2\n");
  eprintf("  c\tlock/unlock center positions\n");
  eprintf("  s\t2nd mouse button adjusts scale & rotation\n");
  eprintf("  o\t3rd mouse button adjusts XY offsets\n");
  eprintf("  v\t2nd mouse button adjusts velocities\n");
  eprintf("  t\t3rd mouse button adjusts image contrast\n");
  eprintf("  Q\t2nd & 3rd mouse buttons set XYVz box\n");
  eprintf("  q\t2nd & 3rd mouse buttons set XY box\n");
  eprintf("  D\tdeletes current box\n");
  eprintf("  L\tlist all box parameters\n");
  eprintf("  j\tcycle through list of boxes\n");
  eprintf("  J\tscan boxes and compute product ball\n");
  eprintf("  V\tscan viewing directions and boxes\n");
  eprintf("  M\tfind maxima on view or product ball\n");
  eprintf("  w\tsave view parameters to file\n");
  eprintf("  W\tsave current display to file\n");
  eprintf("  d\ttoggle display of parameters\n");
  eprintf("  C\ttoggle display of data cube\n");
  eprintf("  Z\ttoggle zoom on XY, XV, or VY panel\n");
  eprintf("  S\tprint out scaling parameters\n");
  eprintf("  b\tgo back in view history\n");
  eprintf("  f\tgo forward in view history\n");
  eprintf("  $\texit idkit program\n");
  eprintf("  <esc>\texit idkit program\n");
}

//  _____________________________________________________
//  storehist: append current parameters to history list.

void storehist(void)
{
  ihist = nhist++;
  old[ihist % MAXHIST] = idk;
}

//  _______________________________________________________________________
//  scanboxes: cycle through boxes, display results, and form intersection.

void scanboxes(void)
{
  matrix vmat;
  real xfield, vfield, voff;
  vector xyoff;
  bool set[2];
  facet *fp;
  box *bp;
  real xbox, ybox, vbox, rbox, hbox, max1, max2, score1, score2, score;

  getpars(vmat, &xfield, NULL, NULL, xyoff, &vfield, &voff,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  for (fp = ball; fp < &ball[nfacet]; fp++)
    fp->prod[0] = fp->prod[1] = 1;
  set[0] = set[1] = FALSE;
  for (bp = &boxtab[1]; bp <= &boxtab[nbox]; bp++) {
    xbox = bp->XYbox[0] * xfield - xyoff[0];
    ybox = bp->XYbox[1] * xfield - xyoff[1];
    vbox = bp->Vbox * vfield - voff;
    rbox = bp->Rbox * xfield;
    hbox = (bp->velflag ? bp->Rbox : 100.0) * vfield;
    if (bp->disk > 0) {
      (void) loadball(vmat, xbox, ybox, vbox, rbox, hbox, bp->disk);
      for (fp = ball; fp < &ball[nfacet]; fp++)
	fp->prod[bp->disk-1] *= fp->spin[bp->disk-1];
      set[bp->disk-1] = TRUE;
    }
  }
  max1 = max2 = 0.0;
  for (fp = ball; fp < &ball[nfacet]; fp++) {
    fp->prod[0] = (set[0] ? fp->prod[0] : 0);	// zero any prod not set
    fp->prod[1] = (set[1] ? fp->prod[1] : 0);
    max1 = MAX(max1, fp->prod[0]);		// find max prod values
    max2 = MAX(max2, fp->prod[1]);
  }
  score1 = (max1 > 0.0 ? rlog10(max1) : -99);	// compute clipped scores
  score2 = (max2 > 0.0 ? rlog10(max2) : -99);
  score = (score1 > -99 && score2 > -99 ? (score1 + score2) / 2 : -99);
  sprintf(msgbuf, "#yscores = %6.2f #c%6.2f #m%6.2f", score, score1, score2);
  lrshow = PRSHOW;				// show product ball
  butbind[0] = JVIEW;
}

//  _________________________________________________
//  scanviews: scan all viewing directions and boxes.

void scanviews(int arg)
{
  facet *vp, *fp;
  real xyzrot[3], xfield, vfield, voff, v1, v2;
  matrix vmat;
  vector xyoff, tmpv;
  bool rejview, try[2], set[2];
  box *bp;
  real xbox, ybox, vbox, rbox, hbox, max1, max2, score1, score2, score;

  vp = &ball[arg * (nfacet / ndirect)];		// set viewing direction
  xyzrot[0] = rasinD(- vp->midp[1]);
  xyzrot[1] = ratan2D(vp->midp[0], vp->midp[2]);
  xyzrot[2] = 0;
  setpars(xyzrot, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  lockcent();					// lock centers in place
  getpars(vmat, &xfield, NULL, NULL, xyoff, &vfield, &voff,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  mapdisks();					// transform center vel's
  v1 = Aux(NthBody(vtab,0));
  v2 = Aux(NthBody(vtab,1));
  rejview = cvelflag && 			// test velocity constraint
            (v1<cvel1[0] || cvel1[1]<v1 || v2<cvel2[0] || cvel2[1]<v2);
  try[0] = try[1] = ! rejview;			// if OK, try both disks
  for (fp = ball; fp < &ball[nfacet]; fp++)
    fp->prod[0] = fp->prod[1] = 1;
  set[0] = set[1] = FALSE;
  for (bp = &boxtab[1]; bp <= &boxtab[nbox]; bp++) {
    xbox = bp->XYbox[0] * xfield - xyoff[0];
    ybox = bp->XYbox[1] * xfield - xyoff[1];
    vbox = bp->Vbox * vfield - voff;
    rbox = bp->Rbox * xfield;
    hbox = (bp->velflag ? bp->Rbox : 100.0) * vfield;
    if (bp->disk > 0 && try[bp->disk-1] &&	// still testing this disk?
	loadball(vmat, xbox, ybox, vbox, rbox, hbox, bp->disk)) {
      for (fp = ball; fp < &ball[nfacet]; fp++)	// multiply latest factor
	fp->prod[bp->disk-1] *= fp->spin[bp->disk-1];
      set[bp->disk-1] = TRUE;			// remember it has been set
    } else if (bp->disk > 0) {
      try[bp->disk-1] = FALSE;			// stop testing this disk
      set[bp->disk-1] = FALSE;			// forget it was ever set
    }
  }
  max1 = max2 = 0;
  for (fp = ball; fp < &ball[nfacet]; fp++) {
    fp->prod[0] = (set[0] ? fp->prod[0] : 0);	// zero any prod never set
    fp->prod[1] = (set[1] ? fp->prod[1] : 0);
    max1 = MAX(max1, fp->prod[0]);		// find max prod values
    max2 = MAX(max2, fp->prod[1]);
  }
  for (fp = vp; fp < vp + (nfacet / ndirect); fp++) {
    fp->view[0] = max1;				// store max to view ball
    fp->view[1] = max2;
  }
  score1 = (max1 > 0 ? rlog10(max1) : -99);	// compute clipped scores
  score2 = (max2 > 0 ? rlog10(max2) : -99);
  score = (score1 > -99 && score2 > -99 ? (score1+score2)/2 : -99);
  sprintf(msgbuf, rejview ? "#yfacet %4d: rejected!" :
	  "#yfacet %4d:  scores = %6.2f #c%6.2f #m%6.2f",
	  arg + 1, score, score1, score2);
  if (arg < ndirect - 1) {
    glutTimerFunc(10, scanviews, arg + 1);	// schedule next cycle
    lrshow = PRSHOW;				// display product ball
  } else {
    animate = FALSE;
    lrshow = VWSHOW;				// display view ball
    butbind[0] = XYROT;
  }
  glutPostRedisplay();
}

//  ______________________________________________________________________
//  findmax: find peaks on view or prod ball, and set corresponding params

void findmax(void)
{
  facet *fp, *mp, *sp1, *sp2;
  real xyzrot[3];

  if (lrshow == VWSHOW) {				// scan view ball?
    mp = NULL;
    for (fp = ball; fp < &ball[nfacet]; fp++)
      if (fp->view[0]*fp->view[1] > (mp != NULL ? mp->view[0]*mp->view[1] : 0))
	mp = fp;
    if (mp != NULL) {
      xyzrot[0] = rasinD(- mp->midp[1]);
      xyzrot[1] = ratan2D(mp->midp[0], mp->midp[2]);
      xyzrot[2] = 0.0;
      setpars(xyzrot, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
      sprintf(msgbuf, "#yscore = %.2f", rlog10(mp->view[0]*mp->view[1]) / 2);
    } else
      sprintf(msgbuf, "#rno peak found");
  } else if (lrshow == PRSHOW) {			// scan prod ball?
    sp1 = sp2 = NULL;
    for (fp = ball; fp < &ball[nfacet]; fp++) {
      if (fp->prod[0] > (sp1 != NULL ? sp1->prod[0] : 0))
	sp1 = fp;
      if (fp->prod[1] > (sp2 != NULL ? sp2->prod[1] : 0))
	sp2 = fp;
    }
    setpars(NULL, NULL,
	    (sp1 != NULL ? sp1->midp : NULL),
	    (sp2 != NULL ? sp2->midp : NULL),
	    NULL, NULL, NULL, NULL);
    sprintf(msgbuf, "#yfound peaks:  #c%s  #m%s",
	    (sp1 != NULL ? "yes" : " no"), (sp2 != NULL ? "yes" : " no"));
  } else
    sprintf(msgbuf, "#rnothing to maximize");
}

//  _________________________________________
//  printscale: print out scaling parameters.

void printscale(void)
{
  real xfield, vfield;
  double Lscl, Vscl, Tscl, Mscl;		// need range for CGS units

  getpars(NULL, &xfield, NULL, NULL, NULL, &vfield, NULL,
	  NULL, NULL, NULL, NULL, NULL, NULL);
  Lscl = (getdparam("xwidth") / xfield) / 3.2408e-22;
  Vscl = (getdparam("vwidth") / vfield) * 1.0e5;
  Tscl = Lscl / Vscl;
  Mscl = Vscl * Vscl * Lscl / 6.672e-8;
  printf("\n\t\tAstronomical\tCGS\n");
  printf("  Length:\t%.2f kpc\t%.3e cm\n", 3.2408e-22 * Lscl, Lscl);
  printf("Velocity:\t%.2f kpc/Gyr\t%.3e cm/s = %.2f km/s\n",
	 1.0226e-5 * Vscl, Vscl, 1.0e-5 * Vscl);
  printf("    Time:\t%.2f Myr\t%.3e s\n", 3.1689e-14 * Tscl, Tscl);
  printf("    Mass:\t%.1f GMsun\t%.3e gm\n\n", 5.025e-43 * Mscl, Mscl);
  sprintf(msgbuf, "#y%.2f kpc  %.1f km/s  %.1f Myr  %.1f GMsun",
	  3.2408e-22*Lscl, 1.0e-5*Vscl, 3.1689e-14*Tscl, 5.025e-43*Mscl); 
}

//  ______________________________________
//  mousebut: process mouse button events.

void mousebut(int but, int state, int x, int y)
{
  int idx;

  if (state == GLUT_DOWN && actpars == -1) {
    storehist();
    idx = (but == GLUT_LEFT_BUTTON ? 0 : but == GLUT_MIDDLE_BUTTON ? 1 : 2);
    actpars = butbind[idx];
    xlast = x;
    ylast = y;
  } else if (state == GLUT_UP && actpars != -1) {
    actlast = actpars;
    actpars = -1;
  } else
    eprintf("mouse button ignored");
}

//  _________________________________________________________________
//  mousemove: process mouse movement; updates params if button down.

void mousemove(int x, int y)
{
  if (actpars != -1) {
    idk.par[actpars][0] += (x - xlast) / ((real) pansize);
    idk.par[actpars][1] += (y - ylast) / ((real) pansize);
    parset[actpars] = TRUE;
    xlast = x;
    ylast = y;
    glutPostRedisplay();
  }
}

//  ________________________________________________________________
//  special: process arrow keys; adjust params last set using mouse.

void special(int key, int x, int y)
{
  int idx, dir;

  idx = (key == GLUT_KEY_LEFT || key == GLUT_KEY_RIGHT ? 0 : 1);
  dir = (key == GLUT_KEY_LEFT || key == GLUT_KEY_UP ? -1 : 1);
  idk.par[actlast][idx] += dir / ((real) pansize);
  parset[actlast] = TRUE;
  glutPostRedisplay();
}

//  __________________________________
//  getjhat: read the jhat input file.

void getjhat(string jhat, string zrot, bool firstcall)
{
  stream istr = stropen(jhat, "r");
  string itags[MaxBodyFields];
  matrix zmat;
  bodyptr bp;

  if (firstcall)
    get_history(istr);
  else
    skip_history(istr);
  if (! (get_snap(istr, &btab, &nbody, &tnow, itags, FALSE) &&
	 set_member(itags, AuxVecTag)))
    error("%s: jhat file %s lacks %s", getargv0(), jhat, AuxVecTag);
  if (! set_member(itags, AuxTag))
    for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
      Aux(bp) = 1.0;
  rotmatrix(zmat, 0.0, 0.0, strtod((zrot), (char **) NULL));
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp))
    MULMV1(AuxVec(bp), zmat, AuxVec(bp));
  strclose(istr);
}

//  ________________________________________
//  getdata: read simulation data from file.

void getdata(string spat, string zrot, bool firstcall)
{
  char namebuf[256];
  stream istr;
  string itags[MaxBodyFields];
  matrix zmat;
  bodyptr bp;

  sprintf(namebuf, spat, idk.step);
  istr = stropen(namebuf, "r?");
  if (istr == NULL) {
    if (firstcall)
      error("%s: can't open file %s\n", getargv0(), namebuf);
    else
      sprintf(msgbuf, "#rcan't open file %s", namebuf);
    return;
  }
  sprintf(msgbuf, "#yreading data file %s", namebuf);
  if (firstcall)
    get_history(istr);
  else
    skip_history(istr);
  if (! (get_snap(istr, &btab, &nbody, &tnow, itags, FALSE) &&
	 set_member(itags, PosTag) && set_member(itags, VelTag)))
    error("%s: required data missing from %s\n", getargv0(), namebuf);
  strclose(istr);
  sprintf(titlebuf, "%s %s", getargv0(), namebuf);
  rotmatrix(zmat, 0.0, 0.0, strtod((zrot), (char **) NULL));
  for (bp = btab; bp < NthBody(btab, nbody); bp = NextBody(bp)) {
    MULMV1(Pos(bp), zmat, Pos(bp));
    MULMV1(Vel(bp), zmat, Vel(bp));
  }
}

//  ______________________________________
//  getcube: read the datacube input file.

void getcube(string cube)
{
  stream istr = stropen(cube, "r");
  real tcube;
  string itags[MaxBodyFields];

  get_history(istr);
  if (! (get_snap(istr, &ctab, &ncube, &tcube, itags, FALSE) &&
	 set_member(itags, PosTag)))
    error("%s: datacube file %s lacks %s", getargv0(), cube, PosTag);
  strclose(istr);
}

//  ___________________________________________
//  getview: read viewing parameters from file.

void getview(void)
{
  stream istr;
  int *pardims, np, nf;
  int vbmask[] = { -16 * ((int) sizeof(real)), 2 * ((int) sizeof(real)), 0 };

  istr = stropen(getparam("initview"), "r");
  get_history(istr);
  if (get_tag_ok(istr, "TitleBuffer"))
    skip_item(istr);
  pardims = get_dimensions(istr, "ViewParam");
  np = MIN(NUMPARS, pardims[0]);
  if (np < NUMPARS)
    eprintf("[getview: reading %d parameter pairs]\n", np);
  get_data(istr, "ViewParam", RealType, idk.par, np, 2, 0);
  get_data(istr, "ProjMatrix", RealType, idk.pmat, NDIM, NDIM, 0);
  if (get_tag_ok(istr, "NBox")) {
    get_data(istr, "NBox", IntType, &nbox, 0);
    if (nbox > MAXBOX)
      error("%s.getview: too many boxes (%d > %d)\n",
	    getargv0(), nbox, MAXBOX);
   get_data(istr, "BoxTab", AnyType, boxtab, nbox+1, sizeof(box), 0);
  }
  if (get_tag_ok(istr, "NFacet")) {
    get_data(istr, "NFacet", IntType, &nf, 0);
    if (nf != nfacet)
      error("%s.getview: facet counts don't match (%d != %d)\n",
	    getargv0(), nf, nfacet);
    get_data_masked(istr, "ViewBall", RealType, ball, nfacet, 2, 0, vbmask);
  }
  strclose(istr);
}

//  __________________________________________
//  putview: write viewing parameters to file.

void putview()
{
  static int nview = 1;
  char namebuf[256];
  struct stat statbuf;
  stream ostr;
  int vbmask[] = { -16 * ((int) sizeof(real)), 2 * ((int) sizeof(real)), 0 };

  if (strchr(getparam("saveview"), '%') != NULL) {	// embedded format?
    do {
      sprintf(namebuf, getparam("saveview"), nview++);	// make a file name
    } while (stat(namebuf, &statbuf) == 0);		// while it exists
  } else						// fixed file name?
    sprintf(namebuf, "%s", getparam("saveview"));	// copy to namebuf
  sprintf(msgbuf, "#y%swriting view %s",
	  (stat(namebuf, &statbuf) == 0 ? "over" : ""), namebuf);
  ostr = stropen(namebuf, "w!");			// open new or old
  put_history(ostr);
  put_string(ostr, "TitleBuffer", titlebuf);
  put_data(ostr, "ViewParam", RealType, idk.par, NUMPARS, 2, 0);
  put_data(ostr, "ProjMatrix", RealType, idk.pmat, NDIM, NDIM, 0);
  if (nbox > 0) {					// write box data?
    put_data(ostr, "NBox", IntType, &nbox, 0);
    put_data(ostr, "BoxTab", AnyType, boxtab, nbox+1, sizeof(box), 0);
  }
  if (lrshow == VWSHOW) {				// write view ball?
    put_data(ostr, "NFacet", IntType, &nfacet, 0);
    put_data_masked(ostr, "ViewBall", RealType, ball, nfacet, 2, 0, vbmask);
  }
  strclose(ostr);
}

//  ___________________________________
//  getimage: read PPM image from file.

void getimage(GLubyte **image, GLubyte **imcpy, string ifile)
{
  stream istr;
  int xsize, ysize, maxv, nbyte, i, c;

  istr = stropen(ifile, "r");
  if (sscanf(inputline(istr), "P6\n") != 0 ||
      sscanf(inputline(istr), "%d %d\n", &xsize, &ysize) != 2 ||
      sscanf(inputline(istr), "%d\n", &maxv) != 1)
    error("%s: error reading image header\n", getargv0());
  if (imgsize == -1)
    imgsize = xsize;
  if (xsize != imgsize || ysize != imgsize || maxv > 255)
    error("%s: image %s not %d by %d by 3 (deep)\n",
	  getargv0(), ifile, imgsize, imgsize);
  nbyte = 3 * xsize * ysize;
  *image = (GLubyte *) allocate(nbyte * sizeof(GLubyte));
  *imcpy = (GLubyte *) allocate(nbyte * sizeof(GLubyte));
  for (i = 0; i < nbyte; i++) {
    c = fgetc(istr);
    if (c == EOF)
      error("%s: unexpected EOF reading image\n", getargv0());
    (*image)[i] = (*imcpy)[i] = (GLubyte) c;
  }
}

//  _______________________________________________________________
//  inputline: read text line from input stream, skipping comments.

string inputline(stream istr)
{
  static char buf[128];

  do {
    if (fgets(buf, 128, istr) == NULL)
      error("%s: unexpected EOF\n", getargv0());
  } while (buf[0] == '\n' || buf[0] == '#');	// skip comments & blanks
  return (buf);
}

//  ______________________________________________
//  putscreen: write image of display to PPM file.

void putscreen(void)
{
  static int nimage = 1;
  char namebuf[256];
  struct stat statbuf;
  GLsizei winsize = 2 * pansize;
  GLubyte *image;
  stream ostr;
  int i, j;

  if (strchr(getparam("saveimage"), '%') != NULL) {
    do {
      sprintf(namebuf, getparam("saveimage"), nimage++);
    } while (stat(namebuf, &statbuf) == 0);
  } else
    sprintf(namebuf, "%s", getparam("saveimage"));
  sprintf(msgbuf, "#y%swriting image %s",
	  (stat(namebuf, &statbuf) == 0 ? "over" : ""), namebuf);
  ostr = stropen(namebuf, "w!");
  sleep(1);
  glViewport(0, 0, winsize, winsize);
  image = (GLubyte *) allocate(3 * winsize * winsize * sizeof(GLubyte));
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, winsize, winsize, GL_RGB, GL_UNSIGNED_BYTE, image);
  fprintf(ostr, "P6\n");
  fprintf(ostr, "%d %d\n", winsize, winsize);
  fprintf(ostr, "%d\n", 255);
  for (j = winsize - 1; j >= 0; j--)
    for (i = 0; i < 3 * winsize; i++)
      fputc((int) image[i + 3 * winsize * j], ostr);
  free(image);
  fclose(ostr);
}

//  _______________________________________
//  clrpars: reset viewing parameter vector

void clrpars(bool allpars)
{
  int i;

  if (allpars) {
    for (i = 0; i < NUMPARS; i++) {
      idk.par[i][0] = idk.par[i][1] = 0.0;
      parset[i] = TRUE;
    }
  } else {
    idk.par[XYROT][0] = idk.par[XYROT][1] = idk.par[SCALE][0] = 0.0;
    parset[XYROT] = parset[SCALE] = TRUE;
  }
}

//  ______________________________________________________
//  getpars: get viewing paramaters from parameter vector.

void getpars(matrix vmat, real *xfield, vector spin1, vector spin2,
	     vector xyoff, real *vfield, real *voff, real *scale, real *bias,
	     vector XYbox, real *Vbox, real *Rbox, matrix jmat)
{
  matrix rmat;

  if (vmat != NULL) {
    rotmatrix(rmat,
	      90 * idk.par[XYROT][1],
	      90 * idk.par[XYROT][0],
	      90 * idk.par[SCALE][0]);
#if defined(SWAPMAT)
    MULM(vmat, idk.pmat, rmat);
#else
    MULM(vmat, rmat, idk.pmat);
#endif
  }
  if (xfield != NULL)
    *xfield = 4.0 * rpow(10.0, idk.par[SCALE][1]);
  if (spin1 != NULL) {
    spin1[0] = rsinD(90*idk.par[DISK1][1]) * rsinD(90*idk.par[DISK1][0]);
    spin1[1] = rsinD(90*idk.par[DISK1][1]) * rcosD(90*idk.par[DISK1][0]);
    spin1[2] = rcosD(90*idk.par[DISK1][1]);
  }
  if (spin2 != NULL) {
    spin2[0] = rsinD(90*idk.par[DISK2][1]) * rsinD(90*idk.par[DISK2][0]);
    spin2[1] = rsinD(90*idk.par[DISK2][1]) * rcosD(90*idk.par[DISK2][0]);
    spin2[2] = rcosD(90*idk.par[DISK2][1]);
  }
  if (xyoff != NULL) {
    xyoff[0] = idk.par[XYOFF][0];
    xyoff[1] = - idk.par[XYOFF][1];
    xyoff[2] = 0.0;
  }
  if (voff != NULL)
    *voff = idk.par[VELOS][0];
  if (vfield != NULL)
    *vfield = 4.0 * rpow(10.0, idk.par[VELOS][1]);
  if (bias != NULL)
    *bias = idk.par[ITRAN][1];
  if (scale != NULL)
    *scale = idk.par[ITRAN][0] + 1.0;
  if (XYbox != NULL) {
    XYbox[0] = idk.par[BOXXY][0];
    XYbox[1] = - idk.par[BOXXY][1];
    XYbox[2] = 0.0;
  }
  if (Vbox != NULL)
    *Vbox = idk.par[BOXVR][0];
  if (Rbox != NULL)
    *Rbox = 0.02 * rpow(10.0, -2.5 * idk.par[BOXVR][1]);
  if (jmat != NULL)
    rotmatrix(jmat, -90 * idk.par[JVIEW][1], -90 * idk.par[JVIEW][0], 0.0);
}

//  ______________________________________________________
//  setpars: set parameter vector from viewing paramaters.

void setpars(real *xyzrot, real *xfield, vector spin1, vector spin2,
	     vector xyoff, vector XYbox, real *Vbox, real *Rbox)
{
  if (xyzrot != NULL) {
    idk.par[XYROT][1] = xyzrot[0] / 90.0;
    idk.par[XYROT][0] = xyzrot[1] / 90.0;
    idk.par[SCALE][0] = xyzrot[2] / 90.0;
    parset[XYROT] = parset[SCALE] = TRUE;
  }
  if (xfield != NULL) {
    idk.par[SCALE][1] = rlog10(*xfield / 4.0);
    parset[SCALE] = TRUE;
  }
  if (spin1 != NULL) {
    idk.par[DISK1][0] = ratan2D(spin1[0], spin1[1]) / 90.0;
    idk.par[DISK1][1] = racosD(spin1[2]) / 90.0;
  }
  if (spin2 != NULL) {
    idk.par[DISK2][0] = ratan2D(spin2[0], spin2[1]) / 90.0;
    idk.par[DISK2][1] = racosD(spin2[2]) / 90.0;
  }
  if (xyoff != NULL) {
    idk.par[XYOFF][0] = xyoff[0];
    idk.par[XYOFF][1] = - xyoff[1];
    parset[XYOFF] = TRUE;
  }
  if (XYbox != NULL) {
    idk.par[BOXXY][0] = XYbox[0];
    idk.par[BOXXY][1] = - XYbox[1];
    parset[BOXXY] = TRUE;
  }
  if (Vbox != NULL) {
    idk.par[BOXVR][0] = *Vbox;
    parset[BOXVR] = TRUE;
  }
  if (Rbox != NULL) {
    idk.par[BOXVR][1] = -0.4 * rlog10(*Rbox / 0.02);
    parset[BOXVR] = TRUE;
  }
}

//  _______________________________________________
//  showtext: display string at specified position.

void showtext(string str, int x, int y)
{
  int i;

  i = 0;
  while (str[i] != (long) NULL)
    if (str[i] != '#') {
      glRasterPos2f((GLfloat) x, (GLfloat) y);
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13, str[i]);
      x += 8;
      i++;
    } else if (str[i+1] == '#') {
      glRasterPos2f((GLfloat) x, (GLfloat) y);
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '#');
      x += 8;
      i += 2;
    } else if (str[i+1] == 'w') {
      glColor3f(1.0, 1.0, 1.0);
      i += 2;
    } else if (str[i+1] == 'y') {
      glColor3f(1.0, 1.0, 0.3);
      i += 2;
    } else if (str[i+1] == 'm') {
      glColor3f(1.0, 0.0, 1.0);
      i += 2;
    } else if (str[i+1] == 'c') {
      glColor3f(0.0, 1.0, 1.0);
      i += 2;
    } else if (str[i+1] == 'r') {
      glColor3f(1.0, 0.2, 0.2);
      i += 2;
    } else if (str[i+1] == 'g') {
      glColor3f(0.2, 1.0, 0.2);
      i += 2;
    } else if (str[i+1] == 'b') {
      glColor3f(0.2, 0.2, 1.0);
      i += 2;
    } else
      error("%s.showtext: undefined escape sequence \"#%c\"\n",
	    getargv0(), str[i+1]);
}

//  _______________________________________________
//  showcross: display cross at specified position.

void showcross(real xmid, real ymid, real size)
{
  glBegin(GL_LINES);
  glVertex2f(xmid - size/2, ymid);
  glVertex2f(xmid + size/2, ymid);
  glVertex2f(xmid, ymid - size/2);
  glVertex2f(xmid, ymid + size/2);
  glEnd(/* GL_LINES */);
}

//  ________________________________________________
//  showbox: display open box at specified position.

void showbox(real xmid, real ymid, real size)
{
  glBegin(GL_LINE_LOOP);
  glVertex2f(xmid - size/2, ymid - size/2);
  glVertex2f(xmid + size/2, ymid - size/2);
  glVertex2f(xmid + size/2, ymid + size/2);
  glVertex2f(xmid - size/2, ymid + size/2);
  glEnd(/* GL_LINE_LOOP */);
}

//  ___________________________________
//  rotmatrix: compute rotation matrix.

void rotmatrix(matrix rmat, real xrot, real yrot, real zrot)
{
  real sx = rsinD(xrot), sy = rsinD(yrot), sz = rsinD(zrot);
  real cx = rcosD(xrot), cy = rcosD(yrot), cz = rcosD(zrot);
  matrix xmat, ymat, zmat, tmp1;

  xmat[0][0] = 1.0;    xmat[0][1] = 0.0;    xmat[0][2] = 0.0;
  xmat[1][0] = 0.0;    xmat[1][1] =  cx;    xmat[1][2] =  sx;
  xmat[2][0] = 0.0;    xmat[2][1] = -sx;    xmat[2][2] =  cx;
  ymat[0][0] =  cy;    ymat[0][1] = 0.0;    ymat[0][2] = -sy;
  ymat[1][0] = 0.0;    ymat[1][1] = 1.0;    ymat[1][2] = 0.0;
  ymat[2][0] =  sy;    ymat[2][1] = 0.0;    ymat[2][2] =  cy;
  zmat[0][0] =  cz;    zmat[0][1] =  sz;    zmat[0][2] = 0.0;
  zmat[1][0] = -sz;    zmat[1][1] =  cz;    zmat[1][2] = 0.0;
  zmat[2][0] = 0.0;    zmat[2][1] = 0.0;    zmat[2][2] = 1.0;
  MULM(tmp1, xmat, ymat);
  MULM(rmat, zmat, tmp1);
}

//  _____________________________________________
//  padlist: pad string list to specified length.

string *padlist(string *old, int len, string name)
{
  int oldlen = xstrlen(old, sizeof(string)) - 1, i;
  string *new;
  
  if (oldlen == 0)
    error("%s: no value given for %s\n", getargv0(), name);
  if (oldlen >= len)
    return (old);
  new = (string *) allocate((len + 1) * sizeof(string));
  for (i = 0; i <= len; i++)
    new[i] = (i < len ? old[MIN(i, oldlen - 1)] : NULL);
  return (new);
}

//  ________________________________________________________________
//  findcenters: extract center positions and velocities from images.

void findcenters(void)
{
  real xmin[2], xmax[2], ymin[2], ymax[2];

  if (xyimage != NULL && centbbox(xmin, xmax, ymin, ymax, xyimage, imgsize)) {
    cpos1[0] = (xmin[0] + xmax[0]) / 2;
    cpos1[1] = (ymin[0] + ymax[0]) / 2;
    cpos2[0] = (xmin[1] + xmax[1]) / 2;
    cpos2[1] = (ymin[1] + ymax[1]) / 2;
    cpos1[2] = cpos2[2] = 0;
    cposflag = TRUE;
  }
  if (xvimage != NULL && centbbox(xmin, xmax, ymin, ymax, xvimage, imgsize)) {
    cvel1[0] = ymin[0];
    cvel1[1] = ymax[0];
    cvel2[0] = ymin[1];
    cvel2[1] = ymax[1];
    cvelflag = TRUE;
  }
  if (vyimage != NULL && centbbox(xmin, xmax, ymin, ymax, vyimage, imgsize)) {
    cvel1[0] = (cvelflag ? (xmin[0] + cvel1[0]) / 2 : xmin[0]);
    cvel1[1] = (cvelflag ? (xmax[0] + cvel1[1]) / 2 : xmax[0]);
    cvel2[0] = (cvelflag ? (xmin[1] + cvel2[0]) / 2 : xmin[1]);
    cvel2[1] = (cvelflag ? (xmax[1] + cvel2[1]) / 2 : xmax[1]);
    cvelflag = TRUE;
  }
}    

//  ___________________________________________________________________
//  centbbox: determine bounding boxes around cyan and magenta markers.

bool centbbox(real xmin[2], real xmax[2], real ymin[2], real ymax[2],
	      GLubyte *img, int size)
{
  int n1 = 0, n2 = 0, ix, iy, red, grn, blu;
  real x, y;

  for (ix = 0; ix < size; ix++)
    for (iy = 0; iy < size; iy++) {
      red = img[3 * (ix + size * iy) + 0];
      grn = img[3 * (ix + size * iy) + 1];
      blu = img[3 * (ix + size * iy) + 2];
      x = ix / ((real) (size - 1)) - 0.5;
      y = 0.5 - iy  / ((real) (size - 1));
      if (red < 1 && grn > 254 && blu > 254) {
	xmin[0] = (n1 == 0 ? x : MIN(xmin[0], x));
	xmax[0] = (n1 == 0 ? x : MAX(xmax[0], x));
	ymin[0] = (n1 == 0 ? y : MIN(ymin[0], y));
	ymax[0] = (n1 == 0 ? y : MAX(ymax[0], y));
	n1++;
      }
      if (red > 254 && grn < 1 && blu > 254) {
	xmin[1] = (n2 == 0 ? x : MIN(xmin[1], x));
	xmax[1] = (n2 == 0 ? x : MAX(xmax[1], x));
	ymin[1] = (n2 == 0 ? y : MIN(ymin[1], y));
	ymax[1] = (n2 == 0 ? y : MAX(ymax[1], y));
	n2++;
      }
    }
  return (n1 > 0 && n2 > 0);
}

//  ________________________________________________________________
//  buildsphere: build approximation to sphere based on icosahedron.

#define SQRT15  0.447213595
#define SQRT45  0.894427191
#define SIN072  0.951056516
#define COS072  0.309016994
#define SIN144  0.587785252
#define COS144 -0.809016994
#define SIN216 -0.587785252
#define COS216 -0.809016994
#define SIN288 -0.951056516
#define COS288  0.309016994

void buildsphere(int lev)
{
  vector vert_a = {                0,                0,       1 };
  vector vert_b = {                0,           SQRT45,  SQRT15 };
  vector vert_c = {  SIN072 * SQRT45,  COS072 * SQRT45,  SQRT15 };
  vector vert_d = {  SIN144 * SQRT45,  COS144 * SQRT45,  SQRT15 };
  vector vert_e = {  SIN216 * SQRT45,  COS216 * SQRT45,  SQRT15 };
  vector vert_f = {  SIN288 * SQRT45,  COS288 * SQRT45,  SQRT15 };
  vector vert_g = {  SIN144 * SQRT45, -COS144 * SQRT45, -SQRT15 };
  vector vert_h = {  SIN072 * SQRT45, -COS072 * SQRT45, -SQRT15 };
  vector vert_i = {                0,          -SQRT45, -SQRT15 };
  vector vert_j = {  SIN288 * SQRT45, -COS288 * SQRT45, -SQRT15 };
  vector vert_k = {  SIN216 * SQRT45, -COS216 * SQRT45, -SQRT15 };
  vector vert_l = {                0,                0,      -1 };

  ball = (facet *) allocate(NUMICOS * (1 << (2 * lev)) * sizeof(facet));
  nfacet = 0;
  buildtriangle(vert_a, vert_b, vert_c, lev);
  buildtriangle(vert_a, vert_c, vert_d, lev);
  buildtriangle(vert_a, vert_d, vert_e, lev);
  buildtriangle(vert_a, vert_e, vert_f, lev);
  buildtriangle(vert_a, vert_f, vert_b, lev);
  buildtriangle(vert_b, vert_c, vert_g, lev);
  buildtriangle(vert_c, vert_d, vert_h, lev);
  buildtriangle(vert_d, vert_e, vert_i, lev);
  buildtriangle(vert_e, vert_f, vert_j, lev);
  buildtriangle(vert_f, vert_b, vert_k, lev);
  buildtriangle(vert_g, vert_h, vert_c, lev);
  buildtriangle(vert_h, vert_i, vert_d, lev);
  buildtriangle(vert_i, vert_j, vert_e, lev);
  buildtriangle(vert_j, vert_k, vert_f, lev);
  buildtriangle(vert_k, vert_g, vert_b, lev);
  buildtriangle(vert_l, vert_g, vert_h, lev);
  buildtriangle(vert_l, vert_h, vert_i, lev);
  buildtriangle(vert_l, vert_i, vert_j, lev);
  buildtriangle(vert_l, vert_j, vert_k, lev);
  buildtriangle(vert_l, vert_k, vert_g, lev);
  if (nfacet != NUMICOS * (1 << (2 * lev)))
    error("%s.buildsphere: facet miscount\n");
}

//  ___________________________________________________
//  buildtriangle: recursively build a triangular face.

void buildtriangle(vector v1, vector v2, vector v3, int lev)
{
  vector vmid, v4, v5, v6;
  real vnorm, a, b, c, x, y;

  if (lev == 0) {				// recursion terminates?
    SETV(ball[nfacet].vert1, v1);
    SETV(ball[nfacet].vert2, v2);
    SETV(ball[nfacet].vert3, v3);
    ADDV(vmid, v1, v2);
    ADDV(vmid, vmid, v3);
    vnorm = absv(vmid);
    DIVVS(vmid, vmid, vnorm);
    SETV(ball[nfacet].midp, vmid);
    nfacet++;					// count another facet
#if defined(LISTFACET)
    a = distv(v1, v2);
    b = distv(v2, v3);
    c = distv(v3, v1);
    x = (a*a + c*c - b*b) / (2 * a);
    y = rsqrt(c*c - x*x);
    printf("%10.8f  %10.8f  %10.8f  %10.8f  %5d\n", a*y/2, a, b, c, nfacet);
#endif
  } else {					// recursion continues?
    ADDV(v4, v1, v2);
    vnorm = absv(v4);
    DIVVS(v4, v4, vnorm);			// find segment midpoints
    ADDV(v5, v2, v3);
    vnorm = absv(v5);
    DIVVS(v5, v5, vnorm);
    ADDV(v6, v3, v1);
    vnorm = absv(v6);
    DIVVS(v6, v6, vnorm);
    buildtriangle(v4, v5, v6, lev - 1);		// middle triangle fist
    buildtriangle(v1, v4, v6, lev - 1);		// then three outer ones
    buildtriangle(v2, v4, v5, lev - 1);
    buildtriangle(v3, v5, v6, lev - 1);
  }
}
