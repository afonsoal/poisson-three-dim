/*
 * SetPolyhedron.h
 *
 *  Created on: Mar 4, 2015
 *      Author: afonsoal
 */

#ifndef SETPOLYHEDRON_H_
#define SETPOLYHEDRON_H_


#include <deal.II/base/point.h>
using namespace dealii;

//template<int dim>
class SetPolyhedron
{

#define MAX_VERTS 100     /* maximum number of polyhedral vertices */
#define MAX_FACES 	 	100     /* maximum number of polyhedral faces */
#define MAX_POLYGON_SZ 	10 /* maximum number of verts per polygonal face */
#define MAX_LINES 5 /* maximum number of LINES per polygonal face Added Afonso*/

private:
public:


	typedef struct line {
		std::vector<Point<3>> X_real; // Vector of points of this line. Size always 2.
		 // Every line will be projected onto a alfa, beta plane, and this is the vector of
		// points for the projected line. The coordinates will be the same as of the real line,
		// excluding one (it can be X,Y or Z)
		std::vector<Point<3>> X_projection;
		// normal vector of the projected line.
		Point<2> norm_projection;
//		double length;
		double length_projection;
		struct FACE *face_of_line;
//		struct polyhedron *poly;
	} LINE;

	typedef struct face {
	  int numVerts, numLines;
	//  double norm[3];
	  Point<3> norm;
	  double w;
	  std::vector<int > verts_idx;
	  std::vector<Point<3> > verts;
	  struct polyhedron *poly;
	  LINE lines[MAX_LINES];
	  Point<3> centroid;
	  // Vector of values A,B,C, corresponding to the value of X,Y,Z of the mapped face.
//	  std::vector<int> map_projection_plane;
	  int A,B,C;
	} FACE;
//	public:
	typedef struct polyhedron {
	  int numVerts, numFaces;
//	  double verts[MAX_VERTS][3];
	  std::vector<Point<3> > verts; // p->verts.resize(MAX_VERTS);
//	  FACE faces[MAX_FACES];
	  struct face faces[MAX_FACES];
	  Point<3> centroid;
	} POLYHEDRON;
	//symbol POLYHEDRON defined as an alias to 'struct polyhedron'

//	POLYHEDRON *p;

//public:
	SetPolyhedron(POLYHEDRON *p);
	~SetPolyhedron();
	void ReadPolyhedron(/*char *name,*/ /*POLYHEDRON *p*/);
	polyhedron *p;

	enum { X = 0, Y = 1, Z = 2};
//	enum { A	, B,	 C};
	void CompLineNormal(FACE *f, LINE *l);
};


#endif /* SETPOLYHEDRON_H_ */
