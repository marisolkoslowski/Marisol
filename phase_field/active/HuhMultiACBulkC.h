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
 * KKSACBulkBase child class for the phase concentration term
 * \f$ - \sum_j \frac{dF_1}{dc_1} \frac{dh_j}{d\eta_i} (c_j) \f$
 * in the the Allen-Cahn bulk residual.
 *
 * The non-linear variable for this Kernel is the order parameter 'eta_i'.
 */
class HuhMultiACBulkC : public HuhMultiACBulkBase
{
public:
  static InputParameters validParams();

  HuhMultiACBulkC(const InputParameters & parameters);

protected:
  virtual Real computeDFDOP(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// Names of phase concentration variables
  MaterialPropertyName _c1_name;

  const VariableValue & _ci;
  unsigned int _ci_var;

  const std::vector<const VariableValue *> _cks;
  const std::vector<unsigned int> _cks_var;

  /// Derivative of the free energy function \f$ \frac d{dc_1} F_1 \f$
  const MaterialProperty<Real> & _prop_dF1dc1;
  /// Second derivative of the free energy function \f$ \frac {d^2}{dc_1^2} F_1 \f$
  //const MaterialProperty<Real> & _prop_d2F1dc12;
  /// Mixed partial derivatives of the free energy function wrt c1 and
  /// any other coupled variables \f$ \frac {d^2}{dc_1 dv} F_1 \f$
  
  
  // std::vector<const MaterialProperty<Real> *> _prop_d2F1dc1darg;
};
