/// Calculates stress for anisortopic crack propagation
/// Includes artificial viscosity and Mie Gruneisen Equation of State

#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

//#include "ComputeStressBase.h"
//#include "HeatSource.h"
#include "MathUtils.h"
//#include "RankTwoTensor.h"
#include "PiecewiseBilinear.h"
//#include "RankFourTensor.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class ComputeMieGrun : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  ComputeMieGrun(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override; 
  virtual void initQpStatefulProperties() override;

  /// Base name prepended to all material property names to allow for
  /// multi-material systems
  const std::string _base_name;

  const MaterialPropertyName _property_name;
  MaterialProperty<Real> & _property;

  const std::string _elasticity_tensor_name;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Our Variables
  const VariableValue & _temperature;

  const std::string       _density_name;
  const std::string _specific_heat_name;

  const MaterialProperty<Real> &        _density;
  const MaterialProperty<Real> &  _specific_heat;

  const MaterialProperty<Real> &              _V;

  const Real                 _Gamma;
  const Real            _slope_UsUp;
  const Real _reference_temperature;

  const MaterialProperty<RankTwoTensor> & _mechanical_strain;

  usingTensorIndices(i_, j_, k_, l_);
};












//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
