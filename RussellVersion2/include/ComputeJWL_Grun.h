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

class ComputeJWL_Grun : public Material //DerivativeMaterialInterface<Material>
{
public:
  ComputeJWL_Grun(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  const MaterialProperty<RankTwoTensor>  & _mechanical_strain;

  const   MaterialPropertyName _property_name;
          MaterialProperty<Real> & _property;
  const   MaterialPropertyName _property_dT_name;
          MaterialProperty<Real> & _property_dT;

  /// Our Variables
  const std::string _specific_heat_name;
  const MaterialProperty<Real> &  _Cv;
  const   VariableValue & _temperature;

  const Real               _A;
  const Real               _B;
  const Real              _R1;
  const Real              _R2;
  const Real           _omega;
  const Real _ref_temperature;

  const MaterialProperty<Real> &        _density;
};


























//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
