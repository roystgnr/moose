#include "OrderParameterFunctionMaterial.h"

template<>
InputParameters validParams<OrderParameterFunctionMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addCoupledVar("eta", "Order parameter variable");
  params.addParam<std::string>("function_name", "f", "actual name for f(eta), i.e. 'h' or 'g'");
  return params;
}

OrderParameterFunctionMaterial::OrderParameterFunctionMaterial(const std::string & name,
                                                               InputParameters parameters) :
    DerivativeMaterialInterface<Material>(name, parameters),
    _eta(coupledValue("eta")),
    _eta_var(coupled("eta")),
    _eta_name(getVar("eta", 0)->name()),
    _function_name(getParam<std::string>("function_name")),
    _prop_f(declareProperty<Real>(_function_name)),
    _prop_df(declarePropertyDerivative<Real>(_function_name, _eta_name)),
    _prop_d2f(declarePropertyDerivative<Real>(_function_name, _eta_name, _eta_name))
{
}
