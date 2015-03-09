/* * cut_cell_integration.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: afonsoal*/



#include "cut_cell_3D.h"
#include <fstream>
#include <iostream>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <iostream>
#include <math.h> // log(), pow
#include <cassert> // log(), pow
//#include "SetPolyhedron.h"
#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/LaplaceBeltramiProblem_general/Problems_3D/PoissonStationary/SetPolyhedron.h"
using namespace dealii;

//#define X 0
//#define Y 1
//#define Z 2

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

// This function will pass all the relevant objects (fe, fe_values, quadrature)
// to private, in order to be accessible to all the methods
// For some strange reason, I could not pass FEValues the same way I passed
// quadrature_formula & fe (to private) . WHY????
// However, I can pass individually for each member function.
//template <int dim>
cut_cell_3D::cut_cell_3D(FEValues<3> const & fe_values,
		FE_Q<3> const& fe,
		QGauss<3> const &quadrature_formula,
		QGauss<1> const &face_quadrature_formula
		,int dim)
:   fe1(fe), quadrature_formula1(quadrature_formula),
		face_quadrature_formula1(face_quadrature_formula),dim(3)

{
	face_quadrature_weights = face_quadrature_formula1.get_weights();
	face_quadrature_points = face_quadrature_formula1.get_points();

//	jacobian_inverse.reinit(2,2);
//	coefficients.reinit(4,4);
//	dofs_per_cell = fe1.dofs_per_cell;
//	n_q_points    = quadrature_formula1.size();
//	for (unsigned int i=0; i<2; ++i)
//		for (unsigned int j=0; j<2; ++j)
//		{
//			jacobian_inverse(i,j) = fe_values.inverse_jacobian(0)[i][j];
////			std::cout << "jacobian_inverse(i,j): "
////			<< fe_values.inverse_jacobian(0)[i][j] << "\n";
//		}
//
//
//	for (unsigned int i = 0;i<4;++i)
//		for (unsigned int j = 0;j<4;++j)
//		{
//			if (i==j || (i==0 && j==3))
//				coefficients(i,j) = 1;
//			else if((i==0 && (j==1 || j==2)) ||
//					(j==3 && (i==1 || i==2)) )
//				coefficients(i,j) = -1;
//			else coefficients(i,j) = 0;
//		}
//	cell_quadrature_weights = quadrature_formula1.get_weights();
//
//	face_quadrature_weights = face_quadrature_formula1.get_weights();
//	face_quadrature_points = face_quadrature_formula1.get_points();
//	n_face_q_points = face_quadrature_formula1.size();
//
//	jacobian_determinant = fe_values.JxW(0)/
//			(cell_quadrature_weights[0] /**cell_quadrature_weights[1]*/ );
////	 cell_quadrature_weights is already wi*wi (= 0.25*0.25 = 0.625)
////	std::cout << "cell_quadrature_weights[0]: " << cell_quadrature_weights[0] << "\n";
////				std::cout << "face_quadrature_weights[0]: " << face_quadrature_weights[0] << "\n";
//	//			std::cout << "fe_values.JxW (0)1: " << fe_values.JxW(0) << "\n";
//	//			std::cout << "fe_values.shape_value (0,0) " << fe_values.shape_value (0,0) << "\n";
//


}
void cut_cell_3D::InitializePolyhedron(SetPolyhedron::POLYHEDRON *p_) { p = p_; }

double cut_cell_3D::return_rhs_face_integration (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i,
		const double face_length)
{
//	Vector<double> exponents_matrix_x(4);
//	Vector<double> exponents_matrix_y(4);
//	exponents_matrix_x(1) = 1;
//	exponents_matrix_x(3) = 1;
//	exponents_matrix_y(2) = 1;
//	exponents_matrix_y(3) = 1;
//
//	double xt;
//	double yt;
//	double final_face_integration = 0;
//	for (unsigned int q_point = 0;q_point<face_quadrature_points.size();++q_point)
//	{
//		xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Point<2> row_sum;
//		double face_integration_scalar = 0;
//		for (unsigned int i = 0;i<4;++i) {
//			Point<2> new_Xt;
//			// First row, representing d phi_i/dx * d phi_j/dx
//			new_Xt[0] = pow(xt,(exponents_matrix_x(i)+1))/(2*(1+exponents_matrix_x(i)))
//																																																													* pow(yt,exponents_matrix_y(i));
//			new_Xt[1] = pow(yt,(exponents_matrix_y(i)+1))/(2*(1+exponents_matrix_y(i)))
//																																																													* pow(xt,exponents_matrix_x(i));
//			new_Xt *= coefficients(dof_i,i);
//			row_sum += new_Xt;
//		}
//
//		face_integration_scalar = row_sum*normal; // point*point = scalar (Is this the correct form?)
//		face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];
//		face_integration_scalar *= face_length;
//		final_face_integration += face_integration_scalar;
//	} // end quadrature sum
//	return final_face_integration;
	return 0;
}

int cut_cell_3D::GetDelta (const int var1, const int var2)
{
	// Return the value of the "delta" multiplier, which acts to select the appropriate multiplication
	// of phi_i, phi_j.
	int delta;

	if (var1 == var2)	delta = 1;
	else 				delta = 0;

	return delta;
}




int cut_cell_3D::CompACoeff (const int index, const int dof, const int var)
{
/*	FullMatrix c_coeff(8,8);
	Vector<doubleint> a_coeff(7);
	a_coeff(0) = delta(0,var)*c_coeff(1,dof)+delta(1,var)*c_coeff(2,dof)+delta(2,var)*c_coeff(3,dof);
	a_coeff(1) = delta(1,var)*c_cpeff(4,dof)+delta(2,var)*c_coeff(6,dof);
	a_coeff(2) = delta(0,var)*c_cpeff(4,dof)+delta(2,var)*c_coeff(5,dof);
	a_coeff(3) = delta(1,var)*c_cpeff(5,dof)+delta(1,var)*c_coeff(6,dof);
	a_coeff(4) = delta(0,var)*c_cpeff(7,dof);
	a_coeff(5) = delta(1,var)*c_cpeff(7,dof);
	a_coeff(6) = delta(2,var)*c_cpeff(7,dof);
	return a_coeff(index); */
	return 0;
}
// Returns face integration for the stiffness (A) matrix.

double cut_cell_3D::return_face_integration (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,
		const int dof_i, const int dof_j, const double face_length)
{
	return 0;
}
Point<2> cut_cell_3D::CompHVector(const double xt, const double yt, const int m, const int n)
{
	/*

	Vector<double> k_p(23);

	k_p(0) = CompACoeff(0,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(1) = CompACoeff(0,dof_i,m)*CompACoeff(1,dof_j,n)+CompACoeff(1,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(2) = CompACoeff(0,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(1,dof_i,m)*CompACoeff(2,dof_j,n)
	+CompACoeff(2,dof_i,m)*CompACoeff(1,dof_j,n)+CompACoeff(6,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(3) = CompACoeff(1,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(3,dof_i,m)*CompACoeff(6,dof_j,n)+
			CompACoeff(4,dof_i,m)*CompACoeff(1,dof_j,n)+CompACoeff(5,dof_i,m)*CompACoeff(2,dof_j,n)+
			CompACoeff(6,dof_i,m)*CompACoeff(3,dof_j,n);
	k_p(4) = CompACoeff(0,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(1,dof_i,m)*CompACoeff(3,dof_j,n)+
			CompACoeff(3,dof_i,m)*CompACoeff(1,dof_j,n)+CompACoeff(5,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(5) = CompACoeff(0,dof_i,m)*CompACoeff(2,dof_j,n)+CompACoeff(2,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(6) = CompACoeff(6,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(2,dof_i,m)*CompACoeff(3,dof_j,n)+
				CompACoeff(3,dof_i,m)*CompACoeff(2,dof_j,n)+CompACoeff(4,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(7) = CompACoeff(0,dof_i,m)*CompACoeff(3,dof_j,n)+CompACoeff(3,dof_i,m)*CompACoeff(0,dof_j,n);
	k_p(8) = CompACoeff(2,dof_i,m)*CompACoeff(2,dof_j,n);
	k_p(9) = CompACoeff(6,dof_i,m)*CompACoeff(1,dof_j,n)+CompACoeff(1,dof_i,m)*CompACoeff(7,dof_j,n);
	k_p(10) = CompACoeff(5,dof_i,m)*CompACoeff(6,dof_j,n)+CompACoeff(6,dof_i,m)*CompACoeff(5,dof_j,n);
	k_p(11) = CompACoeff(6,dof_i,m)*CompACoeff(6,dof_j,n);
	k_p(12) = CompACoeff(1,dof_i,m)*CompACoeff(5,dof_j,n)+CompACoeff(5,dof_i,m)*CompACoeff(1,dof_j,n);
	k_p(13) = CompACoeff(5,dof_i,m)*CompACoeff(5,dof_j,n);
	k_p(14) = CompACoeff(2,dof_i,m)*CompACoeff(2,dof_j,n);
	k_p(15) = CompACoeff(2,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(4,dof_i,m)*CompACoeff(2,dof_j,n);
	k_p(16) = CompACoeff(2,dof_i,m)*CompACoeff(6,dof_j,n)+CompACoeff(6,dof_i,m)*CompACoeff(2,dof_j,n);
	k_p(17) = CompACoeff(3,dof_i,m)*CompACoeff(3,dof_j,n);
	k_p(18) = CompACoeff(3,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(4,dof_i,m)*CompACoeff(3,dof_j,n);
	k_p(19) = CompACoeff(3,dof_i,m)*CompACoeff(5,dof_j,n)+CompACoeff(5,dof_i,m)*CompACoeff(3,dof_j,n);
	k_p(20) = CompACoeff(4,dof_i,m)*CompACoeff(4,dof_j,n);
	k_p(21) = CompACoeff(5,dof_i,m)*CompACoeff(4,dof_j,n)+CompACoeff(4,dof_i,m)*CompACoeff(5,dof_j,n);
	k_p(22) = CompACoeff(4,dof_i,m)*CompACoeff(6,dof_j,n)+CompACoeff(6,dof_i,m)*CompACoeff(4,dof_j,n);

	std::vector<std::vector<int> > alfa(23,3);
//	FullMatrix<double> alfa(23,3);
	int alfa[][3] = {{0,0,0},{1,0,0},{1,1,0},
					  {1,1,1},{1,0,1},{0,1,0},
					  {0,1,1},{0,0,1},{2,0,0},
					  {2,1,0},{2,1,1},{2,2,0},
					  {2,0,1},{2,0,2},{0,2,0},
					  {0,2,1},{1,2,0},{0,0,2},
					  {0,1,2},{1,0,2},{0,2,2},
					  {1,1,2},{1,2,1}};

	Point<2> H;
	Point<2> sum_H;
	// normal vector of the plane F (before projection) = (n_A,n_B,n_C)
	for (unsigned int p = 0; i<22; ++i) {
		beta_1_x= alfa[p][0]+1;
		beta_1_y= alfa[p][0];
		beta_2_x= alfa[p][1];
		beta_2_y= alfa[p][1]+1;
		beta_1_z= alfa[p][0];
		beta_2_z= alfa[p][1];
		if (alfa[p][2] == 0) {
			H(0) = pow(A,beta_1_x+1)*pow(B,beta_2_x)*n_A/(beta_1_x+1)
					+pow(A,beta_1_y+1)*pow(B,beta_2_y)*n_B/(beta_1_y+1)
					+pow(A,beta_1_z+2)*pow(B,beta_2_z)*n_A*n_C/(beta_1_z+2)
					+pow(A,beta_1_z+1)*pow(B,beta_2_z+1)*n_B*n_C/(beta_1_z+1)
					+pow(A,beta_1_z+1)*pow(B,beta_2_z)*n_C/(beta_1_z+1);
			H(1) = H(2) = 0;
		}
		else if (alfa[p][2] == 1) {

			H(0) = pow(A,beta_1_x+2)*pow(B,beta_2_x)*n_A*n_A/(beta_1_x+2)
				+pow(A,beta_1_x+1)*pow(B,beta_2_x+1)*n_A*n_B/(beta_1_x+1)
				+pow(A,beta_1_x+1)*pow(B,beta_2_x)*n_A*w/(beta_1_x+1)
				+pow(A,beta_1_y+2)*pow(B,beta_2_y)*n_A*n_B/(beta_1_y+2)
				+pow(A,beta_1_y+1)*pow(B,beta_2_y+1)*n_B*n_B/(beta_1_y+1)
				+pow(A,beta_1_y+1)*pow(B,beta_2_y)*n_B*w/(beta_1_y+1)
				+pow(A,beta_1_z+3)*pow(B,beta_2_z)*n_A*n_A*n_C/(beta_1_z+3)
				+pow(A,beta_1_z+2)*( pow(B,beta_2_z+1)*2*n_A*n_B*n_C
						+ pow(B,beta_2_z)*2*n_A*n_C*w ) /(beta_1_z+2)
				+pow(A,beta_1_z+1)*( pow(B,beta_2_z+2)*n_B*n_B*n_C
						+ pow(B,beta_2_z)*w*w*n_C+ 2*pow(B,beta_2_z+1)*2*n_B*n_C*w ) /(beta_1_z+1);
			H(1) = H(2) = 0;
		}
		else if (alfa[p][2] == 2) {
			H(0) = 	n_A*pow(A,beta_1_x+3)*pow(B,beta_2_x)*n_A*n_A/(beta_1_x+3)+
					n_A*pow(A,beta_1_x+2)*( pow(B,beta_2_x+1)*2*n_A*n_B
							+ pow(B,beta_2_x)*2*n_A*w )/(beta_1_x+2)
					+n_A*pow(A,beta_1_x+1)*( pow(B,beta_2_x+2)*n_B*n_B +
									pow(B,beta_2_x)*w*w + pow(B,beta_2_x+1)*2*n_B*w )/(beta_1_x+1)
					+// y(beta)-part
					n_B*pow(A,beta_1_y+3)*pow(B,beta_2_y)*n_A*n_A/(beta_1_y+3)+
					n_B*pow(A,beta_1_y+2)*( pow(B,beta_2_y+1)*2*n_A*n_B
							+ pow(B,beta_2_y)*2*n_A*w )/(beta_1_y+2)
					+n_B*pow(A,beta_1_y+1)*( pow(B,beta_2_y+2)*n_B*n_B +
									pow(B,beta_2_y)*w*w + pow(B,beta_2_y+1)*2*n_B*w )/(beta_1_y+1)
					+ // z(gamma)-part
					n_C*pow(A,beta_1_z+4)*pow(B,beta_2_z)*n_A*n_A*n_A/(beta_1_y+4)
					+n_C*pow(A,beta_1_z+3)*( 3*n_A*n_B*pow(B,beta_2_z+1)
							+ pow(B,beta_2_z)*n_A*n_A*w )/(beta_1_y+3)
					+n_C*pow(A,beta_1_z+2)*(
							3*w*w*n_A*pow(B,beta_2_z) +
							+3*n_A*pow(B,beta_2_z+1) +
							+6*n_A*n_B*w*pow(B,beta_2_z+1) ) /(beta_1_z+2)
					+n_C*pow(A,beta_1_z+1)*(
							n_A*n_A*n_A*pow(B,beta_2_z+3) +
							+w*w*w*pow(B,beta_2_z) +
							+3*n_B*n_B*pow(B,beta_2_z+2)
							+3*n_B*w*w*pow(B,beta_2_z+1) ) /(beta_1_z+1) ;
		} // end if alfa 3 == ...
//		H*=k_p[p];
		sum_H+=H*k_p[p];
	} // end loop over p=0,22;

	return sum_H; */
	Point<2> sum_H;
	return sum_H;
}
/* compute various integrations over projection of face */
//Point<2> cut_cell_3D::CompTxxVector ( SetPolyhedron::FACE *f)
//


void cut_cell_3D::compVolumeIntegrals(/*SetPolyhedron::POLYHEDRON *p*/)
{
	// Output checked: A,B,C; f->norm; f->w; a0,b0, a1, b1;
	// This means that my functions to compute the variables above yielded the same result as the
	// computation by Mirtich.

	SetPolyhedron::FACE *f;
	int i, j, n;

	T0 = T1[X] = T1[Y] = T1[Z]
	                        = T2[X] = T2[Y] = T2[Z]
	                                             = TP[X] = TP[Y] = TP[Z] = 0;

	for (i = 0; i < p->numFaces; ++i) {

		f = &p->faces[i];

		int A = f->A;
		int B = f->B;
		int C = f->C;

		std::cout << "\n";
		std::cout << "FACE: " << i << "\n";

		compFaceIntegrals(f);

		T0 += f->norm[X] * ((A == X) ? Fa : ((B == X) ? Fb : Fc));
		//T1[x] = T1[0]

		T1[A] += f->norm[A] * Faa;
		T1[B] += f->norm[B] * Fbb;
		T1[C] += f->norm[C] * Fcc;
		// T2[0] = Txx
		T2[A] += f->norm[A] * Faaa;
		T2[B] += f->norm[B] * Fbbb;
		T2[C] += f->norm[C] * Fccc;

		TP[A] += f->norm[A] * Faab;
		TP[B] += f->norm[B] * Fbbc;
		TP[C] += f->norm[C] * Fcca;
	}

	T1[X] /= 2; T1[Y] /= 2; T1[Z] /= 2;
	T2[X] /= 3; T2[Y] /= 3; T2[Z] /= 3;
	TP[X] /= 2; TP[Y] /= 2; TP[Z] /= 2;

	// Checked output: All results are the same!
	std::cout << "Tx = "  << T1[X] << "\n";
	std::cout << "Ty = "  << T1[Y] << "\n";
	std::cout << "Tz = "  << T1[Z] << "\n";
	std::cout << "\n";
	std::cout << "Txx: " << T2[X] << "\n";
	std::cout << "Tyy: " << T2[Y] << "\n";
	std::cout << "Tzz: " << T2[Z] << "\n";
	std::cout << "\n";
	std::cout << "Txy: " << TP[X] << "\n";
	std::cout << "Tyz: " << TP[Y] << "\n";
	std::cout << "Tzx: " << TP[Z] << "\n";

}

void cut_cell_3D::compFaceIntegrals(SetPolyhedron::FACE *f)
{
  double w;
  double k1, k2, k3, k4;

  std::cout << "Call to compFaceIntegrals \n";
  compProjectionIntegrals(f);

  w = f->w;
  Point<3> n;
  n = f-> norm;

  int A = f->A;
  int B = f->B;
  int C = f->C;

  k1 = 1 / n[C]; k2 = k1 * k1; k3 = k2 * k1; k4 = k3 * k1;

  Fa = k1 * Pa;
  Fb = k1 * Pb;
  Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

  Faa = k1 * Paa;
  Fbb = k1 * Pbb;
  Fcc = k3 * (SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb
	 + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faaa = k1 * Paaa;
  Fbbb = k1 * Pbbb;
  Fccc = -k4 * (CUBE(n[A])*Paaa + 3*SQR(n[A])*n[B]*Paab
	   + 3*n[A]*SQR(n[B])*Pabb + CUBE(n[B])*Pbbb
	   + 3*w*(SQR(n[A])*Paa + 2*n[A]*n[B]*Pab + SQR(n[B])*Pbb)
	   + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

  Faab = k1 * Paab;
  Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
  Fcca = k3 * (SQR(n[A])*Paaa + 2*n[A]*n[B]*Paab + SQR(n[B])*Pabb
	 + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
}

void cut_cell_3D::compProjectionIntegrals(SetPolyhedron::FACE *f)
{
	/* projection integrals */
	std::cout << "Call to compProjectionIntegrals \n";

	double a0, a1, da;
	double b0, b1, db;
	double a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
	double a1_2, a1_3, b1_2, b1_3;
	double C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
	double Cab, Kab, Caab, Kaab, Cabb, Kabb;

	P1 = 0.0;
	Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

	SetPolyhedron::LINE *line;

	std::cout << "f->numVerts: " <<  f->numVerts << "\n";
	for (int i = 0; i < f->numVerts; ++i){
		line = &f->lines[i];

		a0 = f->verts[i][f->A];
		b0 = f->verts[i][f->B];

		a1 = f->verts[(i+1) % f->numVerts][f->A];
		b1 = f->verts[(i+1) % f->numVerts][f->B];


			da = a1 - a0;
			db = b1 - b0;
			a0_2 = a0 * a0; a0_3 = a0_2 * a0; a0_4 = a0_3 * a0;
			b0_2 = b0 * b0; b0_3 = b0_2 * b0; b0_4 = b0_3 * b0;
			a1_2 = a1 * a1; a1_3 = a1_2 * a1;
			b1_2 = b1 * b1; b1_3 = b1_2 * b1;

			C1 = a1 + a0;
			Ca = a1*C1 + a0_2; Caa = a1*Ca + a0_3; Caaa = a1*Caa + a0_4;
			Cb = b1*(b1 + b0) + b0_2; Cbb = b1*Cb + b0_3; Cbbb = b1*Cbb + b0_4;
			Cab = 3*a1_2 + 2*a1*a0 + a0_2; Kab = a1_2 + 2*a1*a0 + 3*a0_2;
			Caab = a0*Cab + 4*a1_3; Kaab = a1*Kab + 4*a0_3;
			Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
			Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

			P1 += db*C1;
			Pa += db*Ca;
			Paa += db*Caa;
			Paaa += db*Caaa;
			Pb += da*Cb;
			Pbb += da*Cbb;
			Pbbb += da*Cbbb;
			Pab += db*(b1*Cab + b0*Kab);
			Paab += db*(b1*Caab + b0*Kaab);
			Pabb += da*(a1*Cabb + a0*Kabb);
	}
		P1 = P1/2.0;
		Pa = Pa/6.0;
		Paa =Paa/ 12.0;
		Paaa = Paaa/ 20.0;
		Pb = Pb/-6.0;
		Pbb = Pbb/-12.0;
		Pbbb = Pbbb/ -20.0;
		Pab = Pab/ 24.0;
		Paab = Paab/ 60.0;
		Pabb = Pabb/ -60.0;
}


Point<2> cut_cell_3D::CompTxxVectorTraditional ( const double alfa_t, const double beta_t,
		SetPolyhedron::FACE *face,
		SetPolyhedron::LINE *line)
{
	// For testing purposes only. Try to calculate the integral Txy over the cube volume to
	// compare with the results obtained by Mirtich (1996).
	// m and n don't matter here.
	Point<2> H_t;
	if (face->C == Y || face->C == Z) {
		if (face->A == X) {
			H_t[0] = pow(alfa_t,4)/4;
//			assert(0);
		}
		if (face->B == X) {
			H_t[0] = pow(beta_t,4)/4;
//			assert(0);
		}
	}
	else if (face->C == X) {
		if (face->A == Z) assert(0);
		H_t[0] =
				pow(line->norm_projection[0],3)*pow(alfa_t,4)/4
				+
				(
					3*pow(line->norm_projection[1]*beta_t,1)
					+
					3*face->w
				)*pow(alfa_t,3)*pow(line->norm_projection[0],2)/3
				+
				(
					3*pow(line->norm_projection[1]*beta_t,2)
					+
					6*line->norm_projection[1]*beta_t*face->w
					+
					3*pow(face->w,2)
				)*pow(alfa_t,2)*pow(line->norm_projection[0],1)/2

				+
				(
					  pow(line->norm_projection[1]*beta_t,3)*pow(face->w,0)
					+
					3*pow(line->norm_projection[1]*beta_t,2)*pow(face->w,1)
					+
					3*pow(line->norm_projection[1]*beta_t,1)*pow(face->w,2)
					+
					  pow(face->w,3)
				)*alfa_t;
		H_t[0] *= pow(-1/face->norm[face->C],3); // * 1/-n_gamma (from h(alfa,beta) = ...
	}
//	else if (face->C == X) {
//		if (face->A == Z) assert(0);
//		H_t[0] =
//				pow(face->norm[face->A],3)*pow(alfa_t,4)/4
//				+
//				(
//					3*pow(face->norm[face->B]*beta_t,1)
//					+
//					3*face->w
//				)*pow(alfa_t,3)*pow(face->norm[face->A],2)/3
//				+
//				(
//					3*pow(face->norm[face->B]*beta_t,2)
//					+
//					6*face->norm[face->B]*beta_t*face->w
//					+
//					3*pow(face->w,2)
//				)*pow(alfa_t,2)*pow(face->norm[face->A],1)/2
//
//				+
//				(
//					  pow(face->norm[face->B]*beta_t,3)*pow(face->w,0)
//					+
//					3*pow(face->norm[face->B]*beta_t,2)*pow(face->w,1)
//					+
//					3*pow(face->norm[face->B]*beta_t,1)*pow(face->w,2)
//					+
//					  pow(face->w,3)
//				)*alfa_t;
//		H_t[0] *= pow(-1/face->norm[face->C],3); // * 1/-n_gamma (from h(alfa,beta) = ...
//	}

	H_t[1] = 0;
	std::cout << "H_t0: "<< H_t << "\n";
	H_t[0] *= face->norm[X]; // * n_x*1/3 *PLANE)
	std::cout << "H_t1: "<< H_t << "\n";
	H_t[0] = H_t[0]*1/3; // * n_x*1/3 *PLANE)
	H_t[0] *= face->norm[face->A]*line->length_projection; // * n_x*1/3
	std::cout << "H_t2: "<< H_t << "\n";
	H_t[0] *= (1/fabs(face->norm[face->C]));

	return H_t;
}

Point<2> cut_cell_3D::CompTyyVector ( SetPolyhedron::FACE *face,
		SetPolyhedron::LINE *line)
{
	// For testing purposes only. Try to calculate the integral Txy over the cube volume to
	// compare with the results obtained by Mirtich (1996).
	// m and n don't matter here.
	Point<2> H_t;

	H_t[0] = pow(line->X_projection[1][Y],4)+

			 pow(line->X_projection[1][Y],3)
			*pow(line->X_projection[0][Y],1)
				+
			 pow(line->X_projection[1][Y],2)
			*pow(line->X_projection[0][Y],2)
				+
			 pow(line->X_projection[1][Y],1)
			*pow(line->X_projection[0][Y],3)
				+
			 pow(line->X_projection[0][Y],4);
	H_t[1] = 0;
	H_t[0] *= face->norm[Y]/3; // * n_x*1/3
	H_t[0] = -H_t[0]/20; // * n_x*1/3
//	H_t[0] *= line->norm_projection[Y]*line->length_projection; // * n_x*1/3

	return H_t;
}

Point<2> cut_cell_3D::CompTxyVector ( SetPolyhedron::FACE *face,
		SetPolyhedron::LINE *line)
{
	// For testing purposes only. Try to calculate the integral Txy over the cube volume to
	// compare with the results obtained by Mirtich (1996).
	// m and n don't matter here.
	Point<2> H_t;

	H_t[0] = pow(line->X_projection[1][X],4)+

			 pow(line->X_projection[1][X],3)
			*pow(line->X_projection[0][X],1)
				+
			 pow(line->X_projection[1][X],2)
			*pow(line->X_projection[0][X],2)
				+
			 pow(line->X_projection[1][X],1)
			*pow(line->X_projection[0][X],3)
				+
			 pow(line->X_projection[0][X],4);
	H_t[1] = 0;
	H_t[0] *= face->norm[1]/2; // * n_x*1/3
	H_t[0] *= face->norm[1]/20; // * n_x*1/3



	return H_t;
}



double cut_cell_3D::CompVolumeIntegral ()
{
	SetPolyhedron::FACE *face;
	SetPolyhedron::LINE *line;
	double alfa_t, beta_t;
	n_q_points = 2;
	std::vector<Point<2> > row_sum(4);
//	std::vector <double > face_quadrature_points_hack;
//	face_quadrature_points_hack.push_back(-pow(1/3,1/2));
//	face_quadrature_points_hack.push_back(+pow(1/3,1/2));
//	int jacobian_determinant = 1;
	double final_integral = 0;
//	for (unsigned int l = 0;l<3;++l) {	// Refers to the element of the InvJac matrix that mult.d phi_i /d m
//		for (unsigned int m = 0;m<3;++m) {   	// Refers to the derivative (X,Y,Z) : d phi_i /d m
//			for (unsigned int n = 0;n<3;++n) {	// Refers to the derivative (X,Y,Z) : d phi_i /d n

//				for (unsigned int p = 0;p<22;++p) { // Already into integral_H... important because
				// the integral formula changes depending on the p! (alfa[][2], to be exact)
				// Get inside a face of the polyhedron!
				double sum_faces = 0;
					for (int it_face = 0;it_face < p->numFaces; ++it_face) {
						// Get information about this face belonging to polyhedron p.
						std::cout << "NEW FACE.......... " << it_face << "\n";
						face = &p->faces[it_face];
						std::cout << "face->verts[ii]: \n";
						for (int ii = 0; ii<face->numVerts; ++ii)
							std::cout << face->verts[ii] << "\n";

//						std::cout << "Test 1 \n";
					    // Get inside a LINE of the FACE!
					    double sum_lines = 0;
						for (int it_line = 0; it_line<face->numLines;++it_line) {
//							&line = plane(it_line);
							line = &face->lines[it_line];
							std::cout << "NEW LINE.......... " << it_line << "\n";
							// l has two points: X_real[0], X_real[1], each with coordinates [0;1;2],
							// corresponding to X,Y,Z. Here, I use the "mappped" A,B,C indices:
							// They will assign the two right projected coordinates to get from this line.
							// Relevant info: line_normal, X0 and X1 (unit coordinates of the line
							// on plane)
//							std::cout << "Test 3 \n";

							std::cout << "line->X_real[0]: "<< line->X_real[0] << "\n";
							std::cout << "line->X_real[1]: "<< line->X_real[1] << "\n";

							std::cout << "line->X_projection[0]: "<< line->X_projection[0] << "\n";
							std::cout << "line->X_projection[1]: "<< line->X_projection[1] << "\n";

							std::cout << "line->X_projection[0]A: "<< line->X_projection[0][face->A] << "\n";
							std::cout << "line->X_projection[0]B: "<< line->X_projection[0][face->B] << "\n";

							std::cout << "line->X_projection[1]A: "<< line->X_projection[1][face->A] << "\n";
							std::cout << "line->X_projection[1]B: "<< line->X_projection[1][face->B] << "\n";

							double sum_quadrature = 0;
							for (int q_point = 0;q_point</*1*/n_q_points /*line quad.,2*/ ;++q_point) {

//								std::cout << "Test 4 \n";

//								alfa_t = line->X_projection[0][X]+(line->X_projection[1][X]
//										- line->X_projection[0][X]) *face_quadrature_points[q_point](0);
//								beta_t = line->X_projection[0][Y]+(line->X_projection[1][Y]
//								        - line->X_projection[0][Y]) *face_quadrature_points[q_point](0);
								// 			 (point order)  v  v (refers to alfa coord)
								alfa_t = line->X_projection[0][face->A]+(line->X_projection[1][face->A]
								       - line->X_projection[0][face->A] )*face_quadrature_points[q_point](0);
								beta_t = line->X_projection[0][face->B]+(line->X_projection[1][face->B]
								       - line->X_projection[0][face->B] )*face_quadrature_points[q_point](0);
								Point<2> H_t;

//								std::cout << "A: " << face->A << "\n";
//								std::cout << "B: " << face->B << "\n";
//								std::cout << "C: " << face->C << "\n";
								std::cout << "w: " << face->w << "\n";

								std::cout << "alfa_t: "<< alfa_t << "\n";
								std::cout << "beta_t: "<< beta_t << "\n";

//								H_t = CompHVector(alfa_t,beta_t,m,n);
//								H_t = CompTxxVector(face,line);
//								H_t = CompTyyVector(face,line);


								H_t = CompTxxVectorTraditional(alfa_t,beta_t,face,line);

								H_t[0] *= face_quadrature_weights[q_point];

								std::cout << "face->norm: "<< face->norm << "\n";
								std::cout << "line->norm_projection: "<< line->norm_projection << "\n";
//								std::cout << "line->length_projection: "<< line->length_projection << "\n";
//								std::cout << "face_quadrature_weights[q_point]: "<< face_quadrature_weights[q_point] << "\n";
								/// DON'T FORGET TO CHANGE ACCORDINGLY |BELOW...
								// X because the formula is Txx = f(x²)
								sum_quadrature += H_t[X];
								std::cout << "sum_quadrature: "<< sum_quadrature << "\n";
							} // End loop quadrature points

							sum_lines += sum_quadrature;
							std::cout << "AFTER sum_lines: "<< sum_lines << "\n";
							std::cout << "\n";
						} // End loop lines
						std::cout << "Final sum of this face: " << sum_lines <<  "\n";
						// Multiply the sum over faces by the negative of the abs. of the Gamma (C)
						// term of the normal vector. This is equivalent to the det(J) of this projection.

						sum_faces +=sum_lines;
//						sum_faces +=sum_lines/(+(face->norm[face->C]));

						std::cout << "fabs(face->norm[face->C]): " << fabs(face->norm[face->C]) << "\n";
						std::cout << "Partial sum of face: " << sum_faces <<  "\n";
						std::cout << "\n";
						std::cout << "\n";
						std::cout << "...................................... " << "\n";
					} // End loop faces
					std::cout << "Final sum of all faces: " << sum_faces <<  "\n";
					final_integral+=sum_faces; // *J(l,n)*J(l,m);
//					total_sum +=integral;
//			} // end loop m
//		} // end loop n
//	} // end loop l (Jacobian)
	// Final volume integral of the polyhedron.
//	final_integral *= jacobian_determinant;

	return final_integral;
}

double cut_cell_3D::getTermC (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const Point<2> &face_normal_vector,
		const int dof_i, const int dof_j)
// Checked if input are really from the new face integration; OK
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//		//		std::cout << "new_Xt: " << new_Xt << "\n";
//		// This will be J⁻¹*fe.shape_grad(dof_i j ,new_Xt)
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		Vector<double> J_fe_shape_grad_j(2);
//		tmp = 0;
//		fe1.shape_grad(dof_j,new_Xt).unroll(tmp);
//		jacobian_inverse.vmult(J_fe_shape_grad_j,tmp);
//		// Checked, multiplication is ok
//		face_integration +=
//				-((J_fe_shape_grad_i[0]*face_normal_vector[0]
//				                                           + J_fe_shape_grad_i[1]*face_normal_vector[1])
//						*fe1.shape_value(dof_j,new_Xt)*face_quadrature_weights[q_point]
//						                                                       *face_length
//			          +
//			          (J_fe_shape_grad_j[0]*face_normal_vector[0]
//			                         + J_fe_shape_grad_j[1]*face_normal_vector[1])
//
//			     *fe1.shape_value(dof_i,new_Xt)*face_quadrature_weights[q_point]*face_length
//				);
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::getTermD (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const int dof_i, const int dof_j, const double alfa)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//				(  fe1.shape_value(dof_i,new_Xt)
//						*fe1.shape_value(dof_j,new_Xt)
//						*face_quadrature_weights[q_point]  )
//						*alfa*face_length ;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

// Before, I used a getTermJ where every local_dof variable was multiplied by two, and it could be used
// only for the ubulk assembly. This one can assemble both usurface and ubulk, with the difference that
// the local_dof_K and local_dof_K_neighbor must be set as if they were a 8 dof cell, ie,
// EXTENDED_local_dof_K   = {0,2,4,6,-1,-1} for usurface
// EXTENDED_local_dof_K_p = {1,3,5,7,-1,-1} for ubulk
// same for neighbor.

double cut_cell_3D::getTermJ (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor,const FEValuesExtractors::Scalar uvariable)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell[uvariable].gradient /*.shape_grad*/(/*2**/local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell[uvariable].gradient/* .shape_grad*/ (/*2**/local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF i");
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell[uvariable].gradient(/*2**/local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell[uvariable].gradient(/*2**/local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values[uvariable].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF j");
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}

//	return j;
	return 0;

}
// Calculate J term when the neighbor cell is inside
double cut_cell_3D::getTermJ_mixed (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	const FEValuesExtractors::Scalar pressure /*usurface */(1);
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF i");
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(/*2**/local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values[pressure].gradient/* .shape_grad*/(/*2**/local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else assert(false && "No common DOF j");
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}
//
//	return j;
	return 0;

}

double cut_cell_3D::getTermJ_OneVar (FEFaceValues<2> const & /*NULL_*/fe_face_values,
		FEFaceValues<2> const & /*NULL_*/fe_face_values_neighborCell,
		/*const*/ int dof_i, /*const */int dof_j,
		const std::vector<int> & local_dof_K,
		const std::vector<int> & local_dof_K_neighbor)
{
//	double j = 0;
//	double dof_i_jump = 0;
//	double dof_j_jump = 0;
//
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> normal = fe_face_values.normal_vector(q_point);
//		// Commom global DOF's for K and K'
//		if ( local_dof_K[dof_i] != -1 && local_dof_K_neighbor[dof_i] != -1) {
//			dof_i_jump =  normal *fe_face_values.shape_grad(local_dof_K[dof_i] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_i] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_i] == -1 && local_dof_K_neighbor[dof_i] != -1){
//			dof_i_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_i] ,q_point)) ;
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_i] == -1 && local_dof_K[dof_i] != -1){
//			dof_i_jump = ( normal *fe_face_values.shape_grad(local_dof_K[dof_i] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else (assert(0));
//
//		// Common global DOF's for K and K'
//		if ( local_dof_K[dof_j] != -1 && local_dof_K_neighbor[dof_j] != -1) {
//			dof_j_jump =  normal *fe_face_values.shape_grad(local_dof_K[dof_j] ,q_point) -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_j] ,q_point);
//		}
//		// Global DOF's are not in the cell K
//		else if ( local_dof_K[dof_j] == -1 && local_dof_K_neighbor[dof_j] != -1){
//
//			dof_j_jump = ( 0 -
//					normal * fe_face_values_neighborCell.shape_grad(local_dof_K_neighbor[dof_j] ,q_point) );
//		}
//		// Global DOF's are not in the cell K'
//		else if ( local_dof_K_neighbor[dof_j] == -1 && local_dof_K[dof_j] != -1){
//
//			dof_j_jump = ( normal *fe_face_values.shape_grad(local_dof_K[dof_j] ,q_point) - 0) ;
//		}
//		// No common global DOF - should not happen
//		else (assert(0));
//
//		j+= /*gamma_1*h**/dof_i_jump*dof_j_jump*fe_face_values.JxW(q_point);
//	}
//	return j;
	return 0;

}

double cut_cell_3D::getTermD2 (const Point<2> &X0, const Point<2> &X1,
		const Point<2> &face_normal_vector,	const int dof_i, const double alfa,
		const double g_D, const double face_length)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//		//		std::cout << "new_Xt: " << new_Xt << "\n";
//		// This will be J⁻¹*fe.shape_grad(dof_i j ,new_Xt)
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		tmp = 0;
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad_i = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		// Checked, multiplication is ok
//
//		face_integration +=
//				(	alfa*fe1.shape_value(dof_i,new_Xt)
//						-
//						(J_fe_shape_grad_i[0]*face_normal_vector[0]
//					   + J_fe_shape_grad_i[1]*face_normal_vector[1])	)
//						*g_D*face_quadrature_weights[q_point]*face_length;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::mass_matrix (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,
		const int dof_i, const int dof_j, const double face_length)
{
//	std::vector<int> exponents_matrix_x =  {0,1,2,0,0,1,2,1,2 };
//	std::vector<int> exponents_matrix_y =  {0,0,0,1,2,1,1,2,2 };
//
//	std::vector<double> new_coefficients =
//	{
//			coefficients(dof_i,0)*coefficients(dof_j,0),
//			coefficients(dof_i,0)*coefficients(dof_j,1)+coefficients(dof_i,1)*coefficients(dof_j,0),
//			coefficients(dof_i,1)*coefficients(dof_j,1),
//			coefficients(dof_i,0)*coefficients(dof_j,2)+coefficients(dof_i,2)*coefficients(dof_j,0),
//			coefficients(dof_i,2)*coefficients(dof_j,2),
//
//			coefficients(dof_i,0)*coefficients(dof_j,3)+coefficients(dof_i,1)*coefficients(dof_j,2) +
//			coefficients(dof_i,2)*coefficients(dof_j,1)+coefficients(dof_i,3)*coefficients(dof_j,0),
//
//			coefficients(dof_i,1)*coefficients(dof_j,3)+coefficients(dof_i,3)*coefficients(dof_j,1),
//			coefficients(dof_i,2)*coefficients(dof_j,3)+coefficients(dof_i,3)*coefficients(dof_j,2),
//			coefficients(dof_i,3)*coefficients(dof_j,3) };
//
//	double xt;
//	double yt;
//	double final_face_integration = 0;
//	for (unsigned int q_point = 0;q_point<face_quadrature_points.size();++q_point)
//	{
//		xt = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		yt = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Point<2> row_sum;
//		double face_integration_scalar = 0;
//		for (unsigned int i = 0;i<9;++i) {
//			Point<2> new_Xt;
//			// First row, representing d phi_i/dx * d phi_j/dx
//			new_Xt[0] = pow(xt,(exponents_matrix_x[i]+1))/(2*(1+exponents_matrix_x[i]))
//									* pow(yt,exponents_matrix_y[i]);
//			new_Xt[1] = pow(yt,(exponents_matrix_y[i]+1))/(2*(1+exponents_matrix_y[i]))
//									* pow(xt,exponents_matrix_x[i]);
//
//			new_Xt *= new_coefficients[i];
//			row_sum += new_Xt;
//		}
//
//		face_integration_scalar = row_sum*normal; // point*point = scalar (Is this the correct form?)
//		face_integration_scalar *= jacobian_determinant*face_quadrature_weights[q_point];
//		face_integration_scalar *= face_length;
//		final_face_integration += face_integration_scalar;
//	} // end quadrature sum
//	return final_face_integration;
	return 0;
}


double cut_cell_3D::getTermConstraintBoundary (const Point<2> &X0, const Point<2> &X1,
		const double face_length, const int dof_i, const int dof_j)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//				(  fe1.shape_value(dof_i,new_Xt)
//						*face_quadrature_weights[q_point]  )
//						*face_length ;
//	} // end quadrature sum
//	return face_integration;
	return 0;
}


double cut_cell_3D::getTermBeltramiBoundary (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
		const double face_length)
{
//	FullMatrix<double> identity_matrix (IdentityMatrix(2));
//	FullMatrix<double> beltrami_operator(2,2);
//	Point<2> beltrami_0;
//	Point<2> beltrami_1;
//	beltrami_operator(0,0) = -normal[0] * normal[0];
//	beltrami_operator(0,1) = -normal[1] * normal[0];
//	beltrami_operator(1,0) = -normal[0] * normal[1];
//	beltrami_operator(1,1) = -normal[1] * normal[1];
//	beltrami_operator.add(1,identity_matrix);
//
//	beltrami_0(0) = beltrami_operator(0,0);
//	beltrami_0(1) = beltrami_operator(0,1);
//	beltrami_1(0) = beltrami_operator(1,0);
//	beltrami_1(1) = beltrami_operator(1,1);
//
//
//
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		Vector<double> J_fe_shape_grad_i(2);
//		// Need to convert tensor fe.shape_grad to vector.
//		Vector<double> tmp(2);
//		fe1.shape_grad(dof_i,new_Xt).unroll(tmp);
//		//J_fe_shape_grad = J⁻¹*fe.shape_grad(dof_i,new_Xt)
//		jacobian_inverse.vmult(J_fe_shape_grad_i,tmp);
//		Vector<double> J_fe_shape_grad_j(2);
//		tmp = 0;
//		fe1.shape_grad(dof_j,new_Xt).unroll(tmp);
//		jacobian_inverse.vmult(J_fe_shape_grad_j,tmp);
//
//
////		 Ignoring Beltrami operator ("=0 for 2D")
//
////		face_integration +=
////						(       (  J_fe_shape_grad_i[0]
////								 * J_fe_shape_grad_j[0] )
////								 +
////								 ( J_fe_shape_grad_i[1]
////		                         * J_fe_shape_grad_j[1] ) )
////				*face_quadrature_weights[q_point]*face_length;
//
//
////		Considering v2 of Beltrami Operator. Here n(x)n is not equal to zero. The results obviously differ
////		from the above.
//		tmp = 0;
//		beltrami_operator.vmult(tmp,J_fe_shape_grad_i);
//		J_fe_shape_grad_i = tmp;
//		tmp = 0;
//		beltrami_operator.vmult(tmp,J_fe_shape_grad_j);
//		J_fe_shape_grad_j = tmp;
//
//		face_integration +=
//				(       (  J_fe_shape_grad_i[0]
//				                             * J_fe_shape_grad_j[0] )
//						+
//						( J_fe_shape_grad_i[1]
//						                    * J_fe_shape_grad_j[1] ) )
//						                    *face_quadrature_weights[q_point]*face_length;
//
//
//	} // end quadrature sum
//	return face_integration;
	return 0;
}

double cut_cell_3D::getTermBeltramiBoundaryRHS (const Point<2> &X0,
		const Point<2> &X1,	const int dof_i,
		const double face_length, const double fs)
{
//	double face_integration = 0;
//		for (int q_point = 0;q_point<n_face_q_points;++q_point)
//		{
//			Point<2> new_Xt;
//			new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//			new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//			face_integration +=
//					( fs*fe1.shape_value(dof_i,new_Xt) )
//					*face_quadrature_weights[q_point]*face_length;
//
//		} // end quadrature sum
//		return face_integration;
	return 0;
}

double cut_cell_3D::constraintVector (const Point<2> &X0,
		const Point<2> &X1,	const int dof_i,
		const double face_length)
{
//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=  ( fe1.shape_value(dof_i,new_Xt) )
//								*face_quadrature_weights[q_point]*face_length;
//
//	} // end quadrature sum
//	return face_integration;
	return 0;

}

double cut_cell_3D::getTermCoupling (const Point<2> &X0,
		const Point<2> &X1,	const Point<2> &normal,	const int dof_i, const int dof_j,
		const double face_length,
		const int multiplier_u_surface_i,
		const int multiplier_u_surface_j,
		const int multiplier_u_bulk_i,
		const int multiplier_u_bulk_j, const double b_B, const double b_S)
{

//	double face_integration = 0;
//	for (int q_point = 0;q_point<n_face_q_points;++q_point)
//	{
//		Point<2> new_Xt;
//		new_Xt[0] = X0(0)+(X1(0)-X0(0))*face_quadrature_points[q_point](0);
//		new_Xt[1] = X0(1)+(X1(1)-X0(1))*face_quadrature_points[q_point](0);
//
//		face_integration +=
//		 (  	b_B*b_B*
//				multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
//				 	 * multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
//				-
//				b_B*b_S*
//				multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
//					 * multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt)
//				-
//				b_B*b_S
//					 * multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt)*multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
//				+
//				b_S*b_S*
//				multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt)
//					 * multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt) )
//					*face_quadrature_weights[q_point]*face_length;
//
//		// This is equivalent to the above
////				( (multiplier_u_bulk_i*fe1.shape_value(dof_i,new_Xt)
////						-  multiplier_u_surface_i*fe1.shape_value(dof_i,new_Xt) )
////						*
////						(multiplier_u_bulk_j*fe1.shape_value(dof_j,new_Xt)
////								-  multiplier_u_surface_j*fe1.shape_value(dof_j,new_Xt)) )
////						*face_quadrature_weights[q_point]*face_length;
//
//	} // end quadrature sum
//	return face_integration;
	return 0;

}

