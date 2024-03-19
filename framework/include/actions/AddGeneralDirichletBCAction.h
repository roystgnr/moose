//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"

#include "libmesh/dirichlet_boundaries.h"

class MooseMesh;

/**
 * This Action adds a strongly-enforced Dirichlet boundary to the
 * problem.  Note that DirichletBoundary objects are not MooseObjects
 * so you need not specify a type for these boundaries.
 */
class AddGeneralDirichletBCAction : public Action
{
public:
  static InputParameters validParams();

  AddGeneralDirichletBCAction(const InputParameters & params);

  virtual void act() override;

protected:
  void addDirichletBoundary(System & s);
};
