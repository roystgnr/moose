#include "HeatSource.h"

template<>
InputParameters validParams<HeatSource>()
{
  InputParameters params = validParams<BodyForce>();

  // Override defaults and documentation, weak form is identical to BodyForce in MOOSE
  params.addParam<Real>("value", 1.0, "Value of heat source. Multiplied by function if present.");
  params.addParam<FunctionName>("function", "1", "Function describing the volumetric heat source");
  return params;
}

HeatSource::HeatSource(const std::string & name, InputParameters parameters)
  :BodyForce(name, parameters)
{
}
