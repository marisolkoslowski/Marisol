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

class SimonGlatzelMelt : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();
  SimonGlatzelMelt(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override; 
  virtual void initQpStatefulProperties() override;

  const Real                 _temp_melt_ref;
  const Real                 _pres_melt_ref;
  const Real                        _a_melt;
  const Real                        _c_melt;
  const MaterialProperty<Real> & _pressure_eos;
        MaterialProperty<Real> & _temp_melt;

  usingTensorIndices(i_, j_, k_, l_);
};
















//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
