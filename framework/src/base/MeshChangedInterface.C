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

#include "MeshChangedInterface.h"

#include "FEProblem.h"

template<>
InputParameters validParams<MeshChangedInterface>()
{
  InputParameters params = emptyInputParameters();
  return params;
}

MeshChangedInterface::MeshChangedInterface(InputParameters & params) :
    _mci_feproblem(*params.get<FEProblem *>("_fe_problem"))
{
  _mci_feproblem.notifyWhenMeshChanges(this);
}
