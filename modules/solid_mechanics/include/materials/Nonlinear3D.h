#ifndef NONLINEAR3D_H
#define NONLINEAR3D_H

#include "Nonlinear.h"

// Forward declarations
class MaterialModel;
class VolumetricModel;

namespace SolidMechanics
{

/**
 * Nonlinear3D is the base class for all 3D nonlinear solid mechanics material models.
 */
class Nonlinear3D : public Nonlinear
{
public:
  Nonlinear3D( SolidModel & solid_model,
               const std::string & name,
               InputParameters parameters );

  virtual ~Nonlinear3D();

protected:

  VariableGradient & _grad_disp_x;
  VariableGradient & _grad_disp_y;
  VariableGradient & _grad_disp_z;
  VariableGradient & _grad_disp_x_old;
  VariableGradient & _grad_disp_y_old;
  VariableGradient & _grad_disp_z_old;

  virtual void computeDeformationGradient( unsigned int qp, ColumnMajorMatrix & F);

  virtual Real volumeRatioOld(unsigned qp) const;

  virtual void computeIncrementalDeformationGradient( std::vector<ColumnMajorMatrix> & Fhat);

};

} // namespace solid_mechanics


#endif
