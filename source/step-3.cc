/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *
 * Authors: Wolfgang Bangerth, 1999,
 *          Guido Kanschat, 2011
 *          Luca Heltai, 2021
 */
#include "step-3.h"

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

Step3::Step3()
  : ParameterAcceptor("Step3")
  , dof_handler(triangulation)
{
  add_parameter("Number of global refinements", n_refinements);
  add_parameter("Output filename", output_name);
  add_parameter("Finite element degree", fe_degree);
  add_parameter("Forcing term expression", forcing_term_expression);
  add_parameter("Boundary condition expression", boundary_contition_expression);
  add_parameter("Problem constants", function_constants);
  add_parameter("Grid generator function", grid_generator_function);
  add_parameter("Grid generator arguments", grid_generator_arguments);
}



void
Step3::make_grid(const std::string &params)
{
  ParameterAcceptor::initialize(params);
  // GridGenerator::hyper_cube(triangulation, -1, 1);
  // GridGenerator::hyper_L(triangulation);
  GridGenerator::generate_from_name_and_arguments(triangulation,
                                                  grid_generator_function,
                                                  grid_generator_arguments);
  // triangulation.begin_active()->face(0)->set_boundary_id(1);
  // // for (auto &face : triangulation.active_face_iterators())
  // //   if (((std::fabs(face->center()[0]) < 1e-12 &&
  // //         std::fabs(face->center()[1] - 0.5) < 1e-12) ||
  // //        (std::fabs(face->center()[1]) < 1e-12 &&
  // //         std::fabs(face->center()[0] - 0.5) < 1e-12)))
  //     face->set_boundary_id(1);
  triangulation.refine_global(n_refinements);
  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}



void
Step3::setup_system()
{
  FE_Q<2> fe{fe_degree};
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);
  system_matrix.reinit(sparsity_pattern);
  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());

  forcing_term.initialize("x,y", forcing_term_expression, function_constants);
  boundary_condition.initialize("x,y",
                                boundary_contition_expression,
                                function_constants);
}



void
Step3::assemble_system()
{
  QGauss<2>          quadrature_formula(fe_degree + 1);
  FEValues<2>        fe_values(dof_handler.get_fe(),
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values |
                          update_quadrature_points);
  const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell_matrix = 0;
      cell_rhs    = 0;
      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx
          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            forcing_term.value(
                              fe_values.quadrature_point(q_index)) * // f(x_q)
                            fe_values.JxW(q_index));                 // dx
        }
      cell->get_dof_indices(local_dof_indices);
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));
      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }
  std::map<types::global_dof_index, double> boundary_values;
  // VectorTools::interpolate_boundary_values(dof_handler,
  //                                          0,
  //                                          Functions::ZeroFunction<2>(),
  //                                          boundary_values);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           boundary_condition,
                                           boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}



void
Step3::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());
}



void
Step3::output_results() const
{
  DataOut<2> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();
  std::ofstream output(output_name);
  data_out.write_vtk(output);
}



void
Step3::run(const std::string &params)
{
  make_grid(params);
  setup_system();
  assemble_system();
  solve();
  output_results();
}
