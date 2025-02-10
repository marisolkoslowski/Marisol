/// Computes pressure from a jones-wilkens-lee equation of state

#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

#include "ComputeStrainBase.h"
//#include "HeatSource.h"

#include "MathUtils.h"
#include "PiecewiseBilinear.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeJWLdependentT : public Material //DerivativeMaterialInterface<Material>
{
public:
  ComputeJWLdependentT(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  const   MaterialPropertyName _property_name;
          MaterialProperty<Real> & _property;

  /// Our Variables
  const   VariableValue & _temperature;

  const Real               _A;
  const Real               _B;
  const Real              _R1;
  const Real              _R2;
  const Real           _omega;
  const Real              _Cv;
  const Real _ref_temperature;

  const MaterialProperty<Real> &              _V;
  const MaterialProperty<Real> &        _density;
};


























//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
