#ifndef TENSORMECHANICSHARDENINGCUTEXPONENTIAL_H
#define TENSORMECHANICSHARDENINGCUTEXPONENTIAL_H

#include "TensorMechanicsHardeningModel.h"

class TensorMechanicsHardeningCutExponential;


template<>
InputParameters validParams<TensorMechanicsHardeningCutExponential>();

/**
 * CutExponential hardening
 * The value = _val_res + (val_0 - val_res)*exp(-rate*(internal_parameter - _intnl_0)), for internal_parameter >= _intnl_0, otherwise value = _val_0
 * Note that while this is not smooth at internal_parameter = _intnl_0,
 * which can produce bad numerical problems.
 */
class TensorMechanicsHardeningCutExponential : public TensorMechanicsHardeningModel
{
 public:
  TensorMechanicsHardeningCutExponential(const std::string & name, InputParameters parameters);

  virtual Real value(const Real & intnl) const;

  virtual Real derivative(const Real & intnl) const;

 private:

  /// The value = _val_res + (val_0 - val_res)*exp(-rate*(internal_parameter - _intnl_0)), for internal_parameter >= _intnl_0, otherwise value = _val_0
  Real _val_0;

  /// The value = _val_res + (val_0 - val_res)*exp(-rate*(internal_parameter - _intnl_0)), for internal_parameter >= _intnl_0, otherwise value = _val_0
  Real _val_res;

  /// The value = _val_res + (val_0 - val_res)*exp(-rate*(internal_parameter - _intnl_0)), for internal_parameter >= _intnl_0, otherwise value = _val_0
  Real _intnl_0;

  /// The value = _val_res + (val_0 - val_res)*exp(-rate*(internal_parameter - _intnl_0)), for internal_parameter >= _intnl_0, otherwise value = _val_0
  Real _rate;

};

#endif // TENSORMECHANICSHARDENINGCUTEXPONENTIAL_H
