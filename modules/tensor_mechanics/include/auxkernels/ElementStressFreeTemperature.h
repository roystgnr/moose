//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "ActivatedElementsMarkerUO.h"

class ElementStressFreeTemperature;

template <>
InputParameters validParams<ElementStressFreeTemperature>();

class ElementStressFreeTemperature : public AuxKernel
{
public:
  /**
   * ElementStressFreeTemperature
   */
  ElementStressFreeTemperature(const InputParameters & parameters);

  virtual ~ElementStressFreeTemperature() {}

protected:
  virtual Real computeValue();

  const VariableValue & _temperature;
  Real _melt_temperature;

  const ActivatedElementsMarkerUO * _marker_uo;
  const std::vector<dof_id_type> * _newly_activated_elem;
};
