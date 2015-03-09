/*
 * SetPolyhedron.cpp
 *
 *  Created on: Mar 4, 2015
 *      Author: afonsoal
 */
#include "SetPolyhedron.h"
#include <fstream>
#include <iostream>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <math.h> // log(), pow
#include <deal.II/base/point.h>
#include <cassert>
#include <string>
#include <math.h> // log(), pow

// Copied from Mirtich; must change appropriately
//#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
//#define MAX_FACES 100     /* maximum number of polyhedral faces */
//#define MAX_POLYGON_SZ 10 /* maximum number of verts per polygonal face */
//#define MAX_LINES 5 /* maximum number of LINES per polygonal face Added Afonso*/

//#define X 0 // DEFINED AS ENUM IN .h file
//#define Y 1
//#define Z 2
#define MAX_POLYGON_SZ 	10 /* maximum number of verts per polygonal face */
#define MAX_VERTS 100     /* maximum number of polyhedral vertices */

using namespace dealii;

SetPolyhedron::SetPolyhedron(POLYHEDRON *p_)
		:	p(p_) { }

SetPolyhedron::~SetPolyhedron(){ /*delete p;*/ }

void SetPolyhedron::ReadPolyhedron(/*char *name,*/ /*POLYHEDRON *p*/)
{
  double dx1, dy1, dz1, dx2, dy2, dz2, nx, ny, nz, len;

  FACE * f = new FACE ();
  LINE * l = new LINE ();

  p->numVerts = 8;
  p->verts.resize(p->numVerts); // Have to change accordingly later;
  std::ifstream read_file("cube_afonso.txt");
  assert(read_file.is_open());
  for (unsigned int i=0; i<8; i++)
  {
	  read_file >> p->verts[i][X] >> p->verts[i][Y] >> p->verts[i][Z];
  }
  read_file.close();


  p->numFaces = 6;
  // Read file of vertices indices; the orientation must match that of the points coordinates.
  std::ifstream read_file2("cube_afonso2.txt");
  assert(read_file2.is_open());
  for (int j = 0; j < p->numFaces; ++j)
  {
	  f = &p->faces[j];
	  // The size has to be changed accordingly later, with the number of vertices of the face.
	  f->verts_idx.resize(MAX_POLYGON_SZ);
//	  	  f->poly = p; // ERROR; cannot convert...
//	   p = ‘SetPolyhedron::POLYHEDRON* {aka SetPolyhedron::polyhedron*}’
//	   f->poly = ‘polyhedron*’  ? Why is it different?
	  read_file2 >> f->verts_idx[0] >> f->verts_idx[1] >> f->verts_idx[2] >> f->verts_idx[3];
  }
  read_file2.close();

  double X_centroid = 0;
  double Y_centroid = 0;
  double Z_centroid = 0;

  for (int i = 0; i < p->numFaces; ++i) {
	  f = &p->faces[i];
	  f->numVerts = 4; // Change accordingly!
	  f->verts.resize(/*MAX_VERTS*/4); // Change accordingly!
	  // Input the coordinates of all vertices for each face. Note that the vertex numbering
	  // "maps" the vertex that are on each face.
	  for (int j = 0; j < f->numVerts; ++j) {
		  f->verts[j] = p->verts[f->verts_idx[j]];
	  }
	  // Compute the centroid of the polyhedron. May be important to define the direction of
	  // the normal vectors.
	  for (int j = 0; j<p->numVerts; ++j) {
		  X_centroid += p->verts[j][X];
		  Y_centroid += p->verts[j][Y];
		  Z_centroid += p->verts[j][Z];
	  }
  }
  Point<3> centroid (X_centroid/f->numVerts,Y_centroid/f->numVerts,Z_centroid/f->numVerts);
  p->centroid = centroid;
  std::cout << "Test p->numFaces: " << p->numFaces  << "\n";
  for (int i = 0; i < p->numFaces; ++i) {
	  f = &p->faces[i];
	  /* compute face normal and offset w from first 3 vertices */
	  dx1 = p->verts[f->verts_idx[1]][X] - p->verts[f->verts_idx[0]][X];
	  dy1 = p->verts[f->verts_idx[1]][Y] - p->verts[f->verts_idx[0]][Y];
	  dz1 = p->verts[f->verts_idx[1]][Z] - p->verts[f->verts_idx[0]][Z];
	  dx2 = p->verts[f->verts_idx[2]][X] - p->verts[f->verts_idx[1]][X];
	  dy2 = p->verts[f->verts_idx[2]][Y] - p->verts[f->verts_idx[1]][Y];
	  dz2 = p->verts[f->verts_idx[2]][Z] - p->verts[f->verts_idx[1]][Z];
	  nx = dy1 * dz2 - dy2 * dz1;
	  ny = dz1 * dx2 - dz2 * dx1;
	  nz = dx1 * dy2 - dx2 * dy1;
	  len = sqrt(nx * nx + ny * ny + nz * nz);
	  f->norm[X] = nx / len;
	  f->norm[Y] = ny / len;
	  f->norm[Z] = nz / len;
	  f->w = - f->norm[X] * p->verts[f->verts_idx[0]][X]
	         - f->norm[Y] * p->verts[f->verts_idx[0]][Y]
	         - f->norm[Z] * p->verts[f->verts_idx[0]][Z];

	  // Compute the centroid of the face. May be important to choose the orientation of the normal vector.
	  double X_centroid = 0;
	  double Y_centroid = 0;
	  double Z_centroid = 0;
	  for (int j = 0; j<f->numVerts; ++j) {
		  X_centroid += f->verts[j][X];
		  Y_centroid += f->verts[j][Y];
		  Z_centroid += f->verts[j][Z];
	  }
	  Point<3> centroid (X_centroid/f->numVerts,Y_centroid/f->numVerts,Z_centroid/f->numVerts);
	  f->centroid = centroid;

	  //Compute the normal vector of each line.
	  // THE normal for 3D line is not defined! It makes sense only for the projected line (in 2D plane)
	  // Get the appropriate mapping of this face to the alfa-beta (A B) plane.
	  nx = fabs(f->norm[X]);
	  ny = fabs(f->norm[Y]);
	  nz = fabs(f->norm[Z]);			// A => alfa, B=> beta, C=>gamma
	  if (nx > ny && nx > nz) f->C = X;  // X = 0, Y = 1, Z = 3 (defined in the beginning)
	  else f->C = (ny > nz) ? Y : Z;
	  f->A = (f->C + 1) % 3;
	  f->B = (f->A + 1) % 3;

	  // Setting up the lines of each face and inputting their vertices.
	  f->numLines = 4;
	  for (int k = 0; k < f->numLines; ++k) {
		  l = &f->lines[k];
		  l->X_real.resize(2); // For this one I am sure the size is 2!
		  l->X_projection.resize(2);
		  // Organize vertices of each line: As there are 4 lines in one face, 2 vertices each,
		  // and the face has 4 vertices, one must duplicate the vertices that "join" lines together.
		  // The order is given by the filAe cube_afonso, and passed to f->verts_idx above.
		  l->X_real[0] =  f->verts[k];
		  l->X_real[1] =  (k != (f->numLines-1)) ? f->verts[k+1] : f->verts[0];

		  // Get the coordinates of the projected line.
		  // Added A,B,C in the left
//		  l->X_projection[0][0] = l->X_real[0][f->A];
//		  l->X_projection[0][1] = l->X_real[0][f->B];
//		  l->X_projection[1][0] = l->X_real[1][f->A];
//		  l->X_projection[1][1] = l->X_real[1][f->B];

		  for (int q=0; q < 2; ++q)
		  {
			  l->X_projection[q][f->A] = l->X_real[q][f->A];
			  l->X_projection[q][f->B] = l->X_real[q][f->B];
			  l->X_projection[q][f->C] = (l->X_real[q][f->A]*f->norm[f->A]+
					  l->X_real[q][f->B]*f->norm[f->B]+f->w )/(-f->norm[f->C]) ;
		  }

		  // Compute the normal face of the projected line and reorders the points accordingly.
		  CompLineNormal(f,l);
	  }
  }

  // Plot each plane in a .txt file for debugging purposes. Visualize in Gnuplot.
  for (int i = 0; i < p->numFaces; ++i) {

	  std::string filename = "PLANE_from_code-";
	  filename += ('0' + i);
	  filename += ".txt";
	  std::ofstream PLANE;
	  PLANE.open(filename.c_str());
	  assert(PLANE.is_open());
	  f = &p->faces[i];
	  for (int k = 0; k < f->numLines; ++k) {
		  l = &f->lines[k];
		  for(int i = 0; i < 3; ++i) {
			  PLANE << l->X_real[0][i] << " ";
		  }
		  PLANE << std::endl;
		  for(int i = 0; i < 3; ++i) {
			  PLANE << l->X_real[1][i] << " ";
		  }
		  PLANE << std::endl;
		  PLANE << std::endl;
	  }
	  PLANE.close();
  }

// Plot the planes projected onto A,B surface. Useful for debugging purposes.
/*  for (int i = 0; i < p->numFaces; ++i) {

	  std::string filename = "Projected_PLANE_from_code-";
	  filename += ('0' + i);
	  filename += ".txt";
	  std::ofstream PLANE;
	  PLANE.open(filename.c_str());
	  assert(PLANE.is_open());
	  f = &p->faces[i];
	  for (int k = 0; k < f->numLines; ++k) {
		  l = &f->lines[k];
		  for(int i = 0; i < 2; ++i) {
			  PLANE << l->X_projection[0][i] << " ";
		  }
		  PLANE << std::endl;
		  for(int i = 0; i < 2; ++i) {
			  PLANE << l->X_projection[1][i] << " ";
		  }
		  PLANE << std::endl;
		  PLANE << std::endl;
	  }
	  PLANE.close();
  }*/

/* Gnuplot code to plot plane by plane.

  set parametric ???
  set autoscale
  set style data points
  set xlabel "x"
  set ylabel "y"
  set zlabel "z"
  */

// splot for [i=0:5] 'PLANE_from_code-'.i.'.txt' w lp title 'Plane '.i


//  delete p;
//  delete f;
//  delete l;
}

void SetPolyhedron::CompLineNormal(FACE *f, LINE *l)
{
	// Computes the normal face of the projected line and reorders the points accordingly.
	// Also computes the length of the projected line.

	Point<2> centroid_projected;
	centroid_projected[0] = f->centroid[f->A];
	centroid_projected[1] = f->centroid[f->B];

//	X_projection[initial point,final point][X,Y coordinate]
	double dy = l->X_projection[1][f->B]-l->X_projection[0][f->B];
	double dx = l->X_projection[1][f->A]-l->X_projection[0][f->A];
	double length = sqrt(dx*dx+dy*dy);

	dy =  dy / length;
	dx =  dx / length;

	l->length_projection = length;

	// If this is the right normal vector, it means that the line is going from X0 to X1!
	Point<2> n1(dy,-dx);
	// If this is the right normal vector, it means that the line is going from X1 to X0!
	// We can rearrange this Vector so that it becomes oriented from X0 to X1. This
	// will fix the parametric equation in return_face_integration.
	Point<2> n2(-dy,dx);


	// Exctract only the relevant (alfa, beta) components of the X_projection points.
	Point<2> aux_n00;
	Point<2> aux_n01;
	aux_n00[0] = l->X_projection[0][f->A];
	aux_n00[1] = l->X_projection[0][f->B];

	aux_n01[0] = l->X_projection[1][f->A];
	aux_n01[1] = l->X_projection[1][f->B];

	// aux_ni represents a point beginning on the middle of the line and being pointed by the normal
	// vector. If the resulting point points outwards the plane, this is the right normal vector.
	Point<2> aux_n1;
	Point<2> aux_n2;
	aux_n1 = (n1+(aux_n00+aux_n01)/2);
//	aux_n1 = (n1+(l->X_projection[0]+l->X_projection[1])/2);
//	Point<2> aux_n2 = (n2+(l->X_projection[0]+l->X_projection[1])/2);
	aux_n2 = (n2+(aux_n00+aux_n01)/2);

	Point<3> Xtemp;

	if ( aux_n1.distance(centroid_projected) > aux_n2.distance(centroid_projected) )
		l->norm_projection = n1;

	else
	{
		l->norm_projection = n2;
////		 Reorder X0, X1
		Xtemp = l->X_projection[0];
		l->X_projection[0] = l->X_projection[1];
		l->X_projection[1] = Xtemp;
	}
}

