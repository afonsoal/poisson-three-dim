/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 */


#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <math.h> // log()

#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/base/derivative_form.h>
// If one wants to use ILU preconditioner
#include <deal.II/lac/sparse_ilu.h>

//Added 14/10
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>

#include <cmath>
#include "SetPolyhedron.h"
//#include "/home/afonsoal/Documents/dealii_8_2_1_unzipped/dealii-8.2.1/my_programs/my_includes/cut_cell_3D.h"
#include "cut_cell_3D.h"


namespace cut_cell_method
{
using namespace dealii;

template <int dim>
class PoissonProblem_3D
{
public:
	PoissonProblem_3D ();
	void run ();
private:
	void make_grid();
	void assemble_system();
	int test;
	Triangulation<dim>     		triangulation;
	FE_Q<dim>           		fe;
	DoFHandler<dim>     	   	dof_handler;
};

template <int dim>
PoissonProblem_3D<dim>::PoissonProblem_3D ()
: fe(1), dof_handler()
{  }


template <int dim>
void PoissonProblem_3D<dim>::make_grid ()
{
	std::cout << "Call to make_grid \n";

	GridGenerator::hyper_cube (triangulation, -2, 2);
//	triangulation.refine_global(1);
	dof_handler.initialize (triangulation,fe);

	dof_handler.distribute_dofs (fe);

	std::string filename_new = "triangulation";
	filename_new += ".eps";

	std::ofstream out (filename_new.c_str());
	GridOut grid_out;
	grid_out.write_eps (triangulation, out);
	std::cout << "New Triangulation created \n";
}

template <int dim>
void PoissonProblem_3D<dim>::assemble_system ()
{
	std::cout << "Call to assemble_system \n";
	// 1st step: Test the formulation developed for integration of a cut cell volume. The
	// first test case will be the simplest possible: integration of phi_i*phi_j in a cube [-2,-2],
	// i = 1, j = 1.

	QGauss<dim>  quadrature_formula(2);
	QGauss<dim-2>  face_quadrature_formula(2); // I am interested in the line integrals

	FEValues<dim> fe_values (fe, quadrature_formula,
			update_values | update_gradients | update_JxW_values
			| update_quadrature_points | update_jacobians |
			update_support_jacobians | update_inverse_jacobians);
	// Just 1 cell in this example...
	int   dofs_per_cell = fe.dofs_per_cell; // 8, ok
	std::cout << "dofs_per_cell: " << dofs_per_cell << "\n";
	for (typename DoFHandler<dim>::active_cell_iterator
			cell = dof_handler.begin_active();
			cell != dof_handler.end(); ++cell)
	{
		std::cout << "cell: " << cell << "\n";
	}

	  SetPolyhedron::POLYHEDRON * p = new  SetPolyhedron::POLYHEDRON();
	  SetPolyhedron Obj_SetPolyhedron(p);
	  Obj_SetPolyhedron.ReadPolyhedron();
	  // Ok; p is ready to be used here.
	  cut_cell_3D Obj_cut_cell_3D
	  (fe_values,FE_Q<dim>(1),quadrature_formula,face_quadrature_formula, dim);

	  Obj_cut_cell_3D.InitializePolyhedron(p);
//	  double volumeIntegral = Obj_cut_cell_3D.CompVolumeIntegral();
//	  double volumeIntegral;
	  std::cout << "volumeIntegral (by Mirtich) [-1,1]^3 = " << +0.333 << std::endl;
//	  std::cout << "volumeIntegral = " << volumeIntegral << std::endl;
	  Obj_cut_cell_3D.compVolumeIntegrals(/*p*/);
//	  Obj_SetPolyhedron.~SetPolyhedron();
	  // BAD; p is still accessible here.
//	  std::cout << "test3: " << p->numVerts << std::endl;
//	  delete p;
	  // BAD2; p is still accessible here, BUT points somewhere (not original p)
//	  std::cout << "test4: " << p->numVerts << std::endl;
}

template <int dim>
void PoissonProblem_3D<dim>::run ()
{
	std::cout<< "Call to run \n";
	make_grid();
	assemble_system();

//	std::cout << "0%3: " << 0%3 << std::endl;
//	std::cout << "1%3: " << 1%3 << std::endl;
//	std::cout << "2%3: " << 2%3 << std::endl;
}
} // End namespace cut_cell_method

int main ()
{
	using namespace dealii;
	using namespace cut_cell_method;
	const unsigned int dim = 3;
	PoissonProblem_3D<dim> Obj_PoissonProblem;
	Obj_PoissonProblem.run ();




}
