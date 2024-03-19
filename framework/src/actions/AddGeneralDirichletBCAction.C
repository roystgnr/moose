//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddGeneralDirichletBCAction.h"

// MOOSE includes
#include "DisplacedProblem.h"
#include "FEProblem.h"
#include "Function.h"
#include "InputParameters.h"
#include "MooseMesh.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/function_base.h"

namespace
{
// A shim to translate a MOOSE Function to a libMesh::FunctionBase
class FunctionToFunction : public FunctionBase<Number>
{
public:
  FunctionToFunction(const Function & f) : _f(&f) { this->_initialized = true; }

  // The DirichletBoundary usage only goes through the multi-component output API
  virtual Number operator()(const Point &, const Real = 0)
  {
    mooseError("Internal error - this code should be unreachable.");
  }

  virtual void operator()(const Point & p, const Real t, DenseVector<Number> & output)
  {
    // We can't clone the underlying function in our clone() method,
    // so we're left locking it to prevent threading bugs.
    Threads::spin_mutex::scoped_lock lock(_function_operator_mutex);

    output.resize(1);
    output.zero();
    output(0) = _f->value(t, p);
  }

  virtual std::unique_ptr<FunctionBase<Number>> clone() const
  {
    return std::make_unique<FunctionToFunction>(*_f);
  }

private:
  const Function * _f;

  Threads::spin_mutex _function_operator_mutex;
};
}

registerMooseAction("MooseApp", AddGeneralDirichletBCAction, "add_general_dirichlet_bc");

InputParameters
AddGeneralDirichletBCAction::validParams()
{
  InputParameters params = Action::validParams();

  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundaries", "Names/ids of boundaries on which to enforce the Dirichlet condition");
  params.addRequiredParam<std::vector<VariableName>>(
      "variables", "Variables for which to apply the Dirichlet condition");

  params.addRequiredParam<std::vector<FunctionName>>(
      "functions", "Functions that specify the Dirichlet values for each variable");

  params.addClassDescription(
      "Action that adds Dirichlet boundary conditions for general continuous variable types");
  return params;
}

AddGeneralDirichletBCAction::AddGeneralDirichletBCAction(const InputParameters & params)
  : Action(params)
{
}

void
AddGeneralDirichletBCAction::addDirichletBoundary(System & sys)
{
  const auto & var_names = getParam<std::vector<VariableName>>("variables");
  const auto & func_names = getParam<std::vector<FunctionName>>("functions");

  auto boundaries = getParam<std::vector<BoundaryName>>("boundaries");
  auto boundary_ids = MooseMeshUtils::getBoundaryIDs(sys.get_mesh(), boundaries, false);

  const std::set<boundary_id_type> bcid_set{boundary_ids.begin(), boundary_ids.end()};

  for (std::size_t i : index_range(var_names))
  {
    const VariableName & var = var_names[i];
    if (!sys.has_variable(var))
      continue;

    unsigned int var_num = sys.variable_number(var);

    // Get the Function we're asked for.  Non-constant, since we'll
    // need to do setup for it.
    Function & fi = _problem->getFunction(func_names[i]);

    // We don't normally setup MOOSE Function objects until *after*
    // init(), but libMesh System::init() is going to need initialized
    // functions to evaluate, so we'll set up the ones we need early.
    for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
      fi.initialSetup();

    // This goes out of scope but it's okay; libMesh takes its clone()
    FunctionToFunction f(fi);

    DirichletBoundary d{bcid_set, {var_num}, f};

    sys.get_dof_map().add_dirichlet_boundary(d);
  }
}

void
AddGeneralDirichletBCAction::act()
{
  auto displaced_problem = _problem->getDisplacedProblem();

  // We'll need scalars to be properly sized early in case any of our
  // functions are ParsedFunctions, since we're about to do function
  // initialSetup early too.
  //
  // We only need thread 0 data for libMesh init() right now, but
  // let's try to future-proof a little.
  for (THREAD_ID tid = 0; tid < libMesh::n_threads(); tid++)
  {
    _problem->reinitScalars(tid);
    if (displaced_problem)
      displaced_problem->reinitScalars(tid);
  }

  auto & eq = _problem->es();
  for (const auto i : make_range(eq.n_systems()))
    this->addDirichletBoundary(eq.get_system(i));
  if (displaced_problem)
  {
    auto & deq = displaced_problem->es();
    for (const auto i : make_range(deq.n_systems()))
      this->addDirichletBoundary(deq.get_system(i));
  }
}
