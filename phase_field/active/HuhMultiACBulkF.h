//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HuhMultiACBulkBase.h"

// Forward Declarations

/**
 * KKSMultiACBulkBase child class for the free energy term
 * \f$ \sum_j \frac{\partial h_j}{\partial \eta_i} F_j + w_i \frac{dg}{d\eta_i} \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter \f$ eta_i \f$.
 */
class HuhMultiACBulkF : public HuhMultiACBulkBase
{
public:
  static InputParameters validParams();

  HuhMultiACBulkF(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);


  /// double well height parameter
  std::vector<MaterialPropertyName> _wik_names;
  std::vector<const MaterialProperty<Real> *> _prop_wik;


};
