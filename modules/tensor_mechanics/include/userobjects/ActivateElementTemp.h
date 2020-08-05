//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"
#include "Function.h"

class ActivateElementTemp;

template <>
InputParameters validParams<ActivateElementTemp>();

class ActivateElementTemp : public ElementUserObject
{
public:
  ActivateElementTemp(const InputParameters & parameters);

  const std::map<dof_id_type, Real> & getActivatedElementsMap() const
  {
    return _activated_elem_map;
  };

  const std::vector<dof_id_type> & getNewlyActivatedElements() const
  {
    return _newly_activated_elem;
  };

  void initialize() override{};
  void execute() override;
  void threadJoin(const UserObject & /*uo*/) override{};
  void finalize() override;

protected:
  std::map<dof_id_type, Real> _activated_elem_map;
  std::vector<dof_id_type> _newly_activated_elem;

  /// activate/inactive subdomain IDs
  const subdomain_id_type _active_subdomain_id;
  /// Spatial tolerance for checking if an element should be activated
  const Real _tol;
  /// path of the heat source, x, y, z components
  const Function & _function_x;
  const Function & _function_y;
  const Function & _function_z;
};
