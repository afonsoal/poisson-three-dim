/*
 * cut_cell_integration.h
 *
 *  Created on: Nov 5, 2014
 *      Author: afonsoal*/


#ifndef CUT_CELL_3D_H_
#define CUT_CELL_3D_H_
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
//#include "SetPolyhedron.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/LaplaceBeltramiProblem_general/Problems_3D/PoissonStationary/SetPolyhedron.h"
/*If you want to declare Obj1 and Obj2 in your .h file, add extern in the .h file like so:

extern SA Obj1, Obj2;
but you should declare the objects in a .cpp file in your project:

main.cpp

SA Obj1, Obj2;

1. Add a forward declaration of library_class_1 in header.h before using it.

class library_class_1;
2. Make sure to #include the .h file that defines library_class_1 in header.cc.*/


using namespace dealii;

class cut_cell_3D
{
private:

	FE_Q<3> fe1;
	QGauss<3> quadrature_formula1;
	QGauss<1> face_quadrature_formula1;
	int   dofs_per_cell;
	int   n_q_points;
	int 	n_face_q_points;
	FullMatrix<double> jacobian_inverse;
	FullMatrix<double> coefficients;
	std::vector<double > cell_quadrature_weights;
	std::vector<double > face_quadrature_weights;
	std::vector<Point<1> > face_quadrature_points;
	double jacobian_determinant;
	int dim;

//	FE_Q<2>& fe1;
//	const unsigned int   n_q_points    = quadrature_formula.size();
//	const unsigned int n_face_q_points = face_quadrature_formula.size();

public:

	cut_cell_3D(FEValues<3> const & fe_values,
			FE_Q<3> const& fe,
			QGauss<3> const &quadrature_formula,
			QGauss<1> const &face_quadrature_formula, const int dim);

	double return_face_integration (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const int dof_j, const double face_length);

	double return_rhs_face_integration (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i,
			const double face_length);

	double getTermC (const Point<2> &X0, const Point<2> &X1,
			const double face_length, const Point<2> &face_normal_vector,
			const int dof_i, const int dof_j);

	double getTermD (const Point<2> &X0, const Point<2> &X1,
				const double face_length, const int dof_i,
				const int dof_j, const double alfa);

	double getTermJ (FEFaceValues<2> const & fe_face_values,
			FEFaceValues<2> const & fe_face_values_neighborCell,
			const int dof_i, const int dof_j,
			const std::vector<int> & local_dof_K,
					const std::vector<int> & local_dof_K_neighbor,
					const FEValuesExtractors::Scalar uvariable );

	double getTermD2(const Point<2> &X0, const Point<2> &X1,
			const Point<2> &face_normal_vector, const int dof_i, const double alfa,
			const double g_D, const double face_length);
	double mass_matrix (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,
			const int dof_i, const int dof_j, const double face_length);
	double getTermConstraintBoundary (const Point<2> &X0, const Point<2> &X1,
			const double face_length, const int dof_i, const int dof_j);
	double getTermBeltramiBoundary (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
			const double face_length);

	double getTermBeltramiBoundaryRHS (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i,
			const double face_length, const double fs);
	double constraintVector (const Point<2> &X0,
			const Point<2> &X1,	const int dof_i,
			const double face_length);

	double getTermCoupling (const Point<2> &X0,
			const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
			const double face_length, const int corrector_u_i, const int corrector_u_j,
			const int corrector_p_i, const int corrector_p_j, const double b_B, const double b_S);

	double getTermJ_mixed (FEFaceValues<2> const & /*NULL_*/fe_face_values,
				FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
				/*const*/ int dof_i, /*const */int dof_j,
				const std::vector<int> & local_dof_K,
				const std::vector<int> & local_dof_K_neighbor);

	double getTermJ_OneVar (FEFaceValues<2> const & /*NULL_*/fe_face_values,
			FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
			/*const*/ int dof_i, /*const */int dof_j,
			const std::vector<int> & local_dof_K,
			const std::vector<int> & local_dof_K_neighbor);

	int CompACoeff (const int index, const int dof, const int var);
	int GetDelta (const int var1, const int var2);
	SetPolyhedron::POLYHEDRON * p;
	void InitializePolyhedron(SetPolyhedron::POLYHEDRON *p_);
	double CompVolumeIntegral ();
	Point<2> CompHVector(const double xt, const double yt, const int m, const int n);
	Point<2> CompTxyVector ( SetPolyhedron::FACE *face,
			SetPolyhedron::LINE *line);
	Point<2> CompTxxVector (SetPolyhedron::FACE *face,
			SetPolyhedron::LINE *line);
	Point<2> CompTyyVector (SetPolyhedron::FACE *face,
				SetPolyhedron::LINE *line);
	Point<2> CompTxxVectorTraditional (const double xt, const double yt, SetPolyhedron::FACE *face,
			SetPolyhedron::LINE *line);

	enum { X = 0, Y = 1, Z = 2};
	/* face integrals */


	void compVolumeIntegrals(/*SetPolyhedron::POLYHEDRON *p*/);
	double T0, T1[3], T2[3], TP[3];

	void compFaceIntegrals(SetPolyhedron::FACE *f);
	void compProjectionIntegrals(SetPolyhedron::FACE *f);

	double Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
	double P1 , Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;
		/* volume integrals */




};

#endif


