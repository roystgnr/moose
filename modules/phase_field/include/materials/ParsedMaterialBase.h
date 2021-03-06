#ifndef PARSEDMATERIALBASE_H
#define PARSEDMATERIALBASE_H

#include "InputParameters.h"

// Forward Declarations
class ParsedMaterialBase;

template<>
InputParameters validParams<ParsedMaterialBase>();

/**
 * Helper class for ParsedMaterial and DerivativeParsedMaterial
 * to declare and read the input parameters.
 */
class ParsedMaterialBase
{
public:
  ParsedMaterialBase(const std::string & name,
                     InputParameters parameters);

protected:
  /// function expression
  std::string _function;

  /// constant vectors
  std::vector<std::string> _constant_names;
  std::vector<std::string> _constant_expressions;

  /// tolerance vectors
  std::vector<std::string> _tol_names;
  std::vector<Real> _tol_values;

  /// material property names
  std::vector<std::string> _mat_prop_names;
};

#endif //PARSEDMATERIALBASE_H
