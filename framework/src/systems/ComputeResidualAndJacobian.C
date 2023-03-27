//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeResidualAndJacobian.h"
#include "FEProblemBase.h"

ComputeResidualAndJacobian::ComputeResidualAndJacobian(FEProblemBase & fe_problem)
  : _fe_problem(fe_problem)
{
}

void
ComputeResidualAndJacobian::residual_and_jacobian(const NumericVector<Number> & u,
                                                  NumericVector<Number> * R,
                                                  SparseMatrix<Number> * J,
                                                  NonlinearImplicitSystem & sys)
{
  mooseAssert(R, "This should be non-null");
  mooseAssert(J, "This should be non-null");

  // libMesh never threw any curve balls here before, and now half of
  // MOOSE seems to break if this identity doesn't hold, so let's have
  // an assertion here rather than weird segfaults or erroneous
  // evaluations later.
  mooseAssert(&u == sys.current_local_solution.get(),
              "MOOSE was asked for a residual not at current_local_solution");

  // If we allow a residual to depend on constrained DoFs, then
  // the corresponding true Jacobian doesn't have a proper zero
  // block for that dependency.
  std::unique_ptr<NumericVector<Number>> u_unconstrained = u.clone();
  sys.get_dof_map().enforce_constraints_exactly(
      sys, sys.current_local_solution.get(), /* homogeneous */ true);

  _fe_problem.computingNonlinearResid(true);
  _fe_problem.computeResidualAndJacobian(*sys.current_local_solution, *R, *J);
  _fe_problem.computingNonlinearResid(false);

  // At this point we've added constraint row terms to the jacobian
  // matrix, but because we've assembled the residual vector from
  // an a priori unspecified set of subresiduals per element we can't
  // safely add heterogeneous terms anywhere there.  So we'll apply
  // those terms once, at the global level.  We need the
  // *unconstrained* soln here because that's where the constraint
  // mismatches are.
  sys.get_dof_map().enforce_constraints_on_residual(sys, R, u_unconstrained.get());
}
