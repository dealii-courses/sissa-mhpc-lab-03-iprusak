#include "step-3.h"

int
main()
{
  deallog.depth_console(2);
  Step3 laplace_problem;
  laplace_problem.run("parameters.prm");
  Step3 laplace_problem_ball;
  laplace_problem_ball.run("parameters/hyper_ball.prm");
  Step3 laplace_problem_cube;
  laplace_problem_cube.run("parameters/hyper_cube.prm");
  Step3 laplace_problem_shell;
  laplace_problem_shell.run("parameters/hyper_shell.prm");
  Step3 laplace_problem_sinbc;
  laplace_problem_sinbc.run("parameters/sin_bc.prm");
  Step3 laplace_problem_sinrhs;
  laplace_problem_sinrhs.run("parameters/sin_rhs.prm");
  return 0;
}