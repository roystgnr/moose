/*****************************************/
/* Written by andrew.wilkins@csiro.au    */
/* Please contact me if you make changes */
/*****************************************/

//  This post processor returns the derivative of density wrt pressure
//
#include "RichardsDensityPrimeAux.h"

template<>
InputParameters validParams<RichardsDensityPrimeAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("pressure_var", "The variable that represents the pressure");
  params.addRequiredParam<UserObjectName>("density_UO", "Name of user object that defines the density.");
  params.addClassDescription("auxillary variable which is d(density)/dp");
  return params;
}

RichardsDensityPrimeAux::RichardsDensityPrimeAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    _pressure_var(coupledValue("pressure_var")),
    _density_UO(getUserObject<RichardsDensity>("density_UO"))
{}

Real
RichardsDensityPrimeAux::computeValue()
{
  return _density_UO.ddensity(_pressure_var[_qp]);
}
