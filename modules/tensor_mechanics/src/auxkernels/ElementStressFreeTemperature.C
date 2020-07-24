//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementStressFreeTemperature.h"

registerMooseObject("TensorMechanicsApp", ElementStressFreeTemperature);

template <>
InputParameters
validParams<ElementStressFreeTemperature>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Determine activated elements.");
  params.addRequiredCoupledVar("temp_aux",
                               "Temperature aux variable used to determine activated elements.");
  params.addRequiredParam<Real>("melt_temperature", "Melt temperature.");
  params.addRequiredParam<UserObjectName>("marker_uo", "Marker UserObject");
  return params;
}

ElementStressFreeTemperature::ElementStressFreeTemperature(const InputParameters & parameters)
  : AuxKernel(parameters),
    _temperature(coupledValue("temperature")),
    _melt_temperature(getParam<Real>("melt_temperature")),
    _marker_uo(isParamValid("marker_uo") ? &getUserObjectByName<ActivatedElementsMarkerUO>(
                                           getParam<UserObjectName>("marker_uo"))
                                           : nullptr)
{
  if (_marker_uo)
    _newly_activated_elem = &(_marker_uo->getNewlyActivatedElements());
  else
    _newly_activated_elem = nullptr;
}

Real
ElementStressFreeTemperature::computeValue()
{
  if (!isNodal()) // not sure if should be nodal or element variable
    mooseError("must run on a nodal variable");

  // check if element is newly activated
  auto iter = std::find(_newly_activated_elem->begin(), _newly_activated_elem->end(), _current_elem->id());
  // if element is just activated: stress_free_temp =  current temperature
  // else: stress_free_temp = stress free temperature at previous step
  if (iter != _newly_activated_elem->end())
    return _temperature[_qp];
  else
    return _u_old[_qp];
}
