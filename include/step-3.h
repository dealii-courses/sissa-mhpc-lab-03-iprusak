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

// Make sure we don't redefine things
#ifndef step3_include_file
#define step3_include_file



#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_acceptor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>



// Forward declare the tester class
class Step3Tester;

using namespace dealii;
class Step3 : ParameterAcceptor
{
public:
  Step3();
  void
  run(const std::string &params);
  void
  initialize();

protected:
  void
  make_grid(const std::string &params);
  void
  setup_system();
  void
  assemble_system();
  void
  solve();
  void
  output_results() const;

  Triangulation<2>     triangulation;
  DoFHandler<2>        dof_handler;
  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;
  Vector<double>       solution;
  Vector<double>       system_rhs;

  FunctionParser<2> forcing_term;
  FunctionParser<2> boundary_condition;

  friend class Step3Tester;

  unsigned int                  n_refinements                 = 3;
  std::string                   output_name                   = "solution.vtk";
  unsigned int                  fe_degree                     = 1;
  std::string                   forcing_term_expression       = "1";
  std::string                   boundary_contition_expression = "0";
  std::map<std::string, double> function_constants;
  std::string                   grid_generator_function  = "hyper_L";
  std::string                   grid_generator_arguments = "-1: 1: false";
};

#endif