/*****************************************/
/* Written by andrew.wilkins@csiro.au    */
/* Please contact me if you make changes */
/*****************************************/

//  "PowerGas" form of relative permeability
//
#include "RichardsRelPermPowerGas.h"

template<>
InputParameters validParams<RichardsRelPermPowerGas>()
{
  InputParameters params = validParams<RichardsRelPerm>();
  params.addRequiredRangeCheckedParam<Real>("simm", "simm >= 0 & simm < 1", "Immobile saturation.  Must be between 0 and 1.   Define s = (seff - simm)/(1 - simm).  Then relperm = 1 - (n+1)(1-s)^n + n(1-s)^(n+1)");
  params.addRequiredRangeCheckedParam<Real>("n", "n >= 2", "Exponent.  Must be >= 2.   Define s = (seff - simm)/(1 - simm).  Then relperm = 1 - (n+1)(1-s)^n + n(1-s)^(n+1)");
  params.addClassDescription("Power form of relative permeability that might be useful for gases.  Define s = (seff - simm)/(1 - simm).  Then relperm = 1 - (n+1)(1-s)^n + n(1-s)^(n+1) if s<simm, otherwise relperm=1");
  return params;
}

RichardsRelPermPowerGas::RichardsRelPermPowerGas(const std::string & name, InputParameters parameters) :
    RichardsRelPerm(name, parameters),
    _simm(getParam<Real>("simm")),
    _n(getParam<Real>("n"))
{
}


Real
RichardsRelPermPowerGas::relperm(Real seff) const
{
  if (seff >= 1.0)
    return 1.0;

  if (seff <= _simm)
    return 0.0;

  Real s_internal = (seff - _simm)/(1.0 - _simm);
  Real krel = 1 - (_n + 1)*std::pow(1 - s_internal, _n) + _n*std::pow(1 - s_internal, _n + 1);

  // bound, just in case
  if (krel < 0) { krel = 0;}
  if (krel > 1) { krel = 1;}
  return krel;
}


Real
RichardsRelPermPowerGas::drelperm(Real seff) const
{
  if (seff >= 1.0)
    return 0.0;

  if (seff <= _simm)
    return 0.0;

  Real s_internal = (seff - _simm)/(1.0 - _simm);
  Real krelp = (_n + 1)*_n*std::pow(1 - s_internal, _n - 1) - _n*(_n + 1)*std::pow(1- s_internal, _n);
  return krelp/(1.0 - _simm);
}


Real
RichardsRelPermPowerGas::d2relperm(Real seff) const
{
  if (seff >= 1.0)
    return 0.0;

  if (seff <= _simm)
    return 0.0;

  Real s_internal = (seff - _simm)/(1.0 - _simm);
  Real krelpp = -(_n + 1)*_n*(_n - 1)*std::pow(1 - s_internal, _n - 2) + _n*(_n + 1)*_n*std::pow(1 - s_internal, _n - 1);
  return krelpp/std::pow(1.0 - _simm, 2);
}

