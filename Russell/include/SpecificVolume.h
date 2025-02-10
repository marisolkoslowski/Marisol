#pragma once

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"

#include "MathUtils.h"
#include "PiecewiseBilinear.h"
#include "ExternalPetscSolverApp.h"
#include "petscblaslapack.h"

class SpecificVolume : public Material //DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  SpecificVolume(const InputParameters & parameters);

protected:
  virtual void      computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

        MaterialProperty<Real> &  _bulk_modulus;
        MaterialProperty<Real> &       _ss_prop;
        MaterialProperty<Real> &             _V;
        MaterialProperty<Real> &  _density_corr;
        MaterialProperty<Real> &            _mu;
  const MaterialProperty<Real> &       _density;
  const MaterialProperty<RankTwoTensor>  & _mechanical_strain;
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;
};












//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
