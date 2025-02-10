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

#include <MaterialPropertyInterface.h>

class PressureAddUp : public Material
{
public:
  static InputParameters validParams();
  PressureAddUp(const InputParameters & parameters);

protected:
  virtual void      computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  const MaterialPropertyName _property_name;
  MaterialProperty<Real> & _property;

  std::vector<std::string> _sum_materials;
  unsigned int _num_materials;

  std::vector<const VariableValue *> _values;
  const unsigned int _n_values;

  std::vector<const MaterialProperty <Real> *> _summand_F;


};













//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
