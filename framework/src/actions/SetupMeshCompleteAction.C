/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SetupMeshCompleteAction.h"
#include "MooseMesh.h"
#include "Moose.h"
#include "Adaptivity.h"
#include "MooseApp.h"

template<>
InputParameters validParams<SetupMeshCompleteAction>()
{
  InputParameters params = validParams<Action>();
  return params;
}

SetupMeshCompleteAction::SetupMeshCompleteAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

bool
SetupMeshCompleteAction::completeSetup(MooseMesh *mesh)
{
  bool prepared = mesh->prepared();

  if (!prepared)
  {
    Moose::setup_perf_log.push("Prepare Mesh","Setup");
    mesh->prepare();
    Moose::setup_perf_log.pop("Prepare Mesh","Setup");
  }

  return prepared;
}

void
SetupMeshCompleteAction::act()
{
  if (!_mesh)
    mooseError("No mesh file was supplied and no generation block was provided");

  /**
   * If possible we'd like to refine the mesh here before the equation systems
   * are setup to avoid doing expensive projections. If however we are doing a
   * file based restart and we need uniform refinements, we'll have to postpone
   * those refinements until after the solution has been read in.
   */
  if (_current_task == "uniform_refine_mesh" && _app.setFileRestart() == false && _app.isRecovering() == false)
  {
    Adaptivity::uniformRefine(_mesh.get());

    if (_displaced_mesh)
      Adaptivity::uniformRefine(_displaced_mesh.get());
  }
  else
  {
    completeSetup(_mesh.get());

    if (_displaced_mesh)
      completeSetup(_displaced_mesh.get());
  }
}
