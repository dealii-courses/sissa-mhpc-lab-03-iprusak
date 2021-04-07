#include <gtest/gtest.h>

#include <fstream>

#include "step-3.h"

using namespace dealii;

// Test Fixture for step-3
class Step3Tester : public ::testing::Test, public Step3
{
public:
  Step3Tester() = default;
};


TEST_F(Step3Tester, HyperCube)
{
  make_grid("parameters/hyper_cube.prm");
}

TEST_F(Step3Tester, HyperBall)
{
  make_grid("parameters/hyper_ball.prm");
}


TEST_F(Step3Tester, HyperShell)
{
  make_grid("parameters/hyper_shell.prm");
}
TEST_F(Step3Tester, SinBC)
{
  make_grid("parameters/sin_bc.prm");
}
TEST_F(Step3Tester, SinRHS)
{
  make_grid("parameters/sin_rhs.prm");
}
