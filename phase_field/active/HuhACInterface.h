//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * Compute the Allen-Cahn interface term with constant Mobility and Interfacial parameter
 */
class HuhACInterface : public Kernel
{
public:
  static InputParameters validParams();

  HuhACInterface(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);


  /// Step Function
  const MaterialProperty<Real> & _si;

  std::vector<MaterialPropertyName> _sk_names;
  std::vector<const MaterialProperty<Real> *> _prop_sk;
  /// Mobilities
  std::vector<MaterialPropertyName> _Mik_names;
  std::vector<const MaterialProperty<Real> *> _prop_Mik;

  unsigned int _num_k;
  
  /// Interfacial parameter
  std::vector<MaterialPropertyName> _Epsik_names;
  std::vector<const MaterialProperty<Real> *> _prop_Epsik;

  /// Gradient of the coupled variable
  // std::vector<VariableName> _etak_names;
  std::vector<const VariableGradient *> _grad_etak;
  std::vector<unsigned int> _etak_var;
  /// Index of the coupled variable

  const MaterialProperty<Real> & _num_phases;
};
