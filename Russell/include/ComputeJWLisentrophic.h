/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

#include "ComputeStressBase.h"
//#include "HeatSource.h"

#include "MathUtils.h"
#include "PiecewiseBilinear.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeJWLisentrophic : public Material //public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  ComputeJWLisentrophic(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override; 
  virtual void initQpStatefulProperties() override;

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  const MaterialPropertyName _property_name;
  MaterialProperty<Real> & _property;

  /// Our Variables
  const VariableValue & _temperature;

  const Real     _A;
  const Real     _B;
  const Real     _C;
  const Real    _R1;
  const Real    _R2;
  const Real _omega;

  const MaterialProperty<Real> &   _bulk_modulus;
  const MaterialProperty<Real> &              _V;
};


























//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
