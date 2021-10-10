//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HuhACBulk.h"

// Forward Declarations

/**
 * ACBulk child class that sets up necessary variables and materials for
 * calculation of residual contribution \f$ \frac{\partial f}{\partial \eta_i} \f$
 * by child classes KKSMultiACBulkF and KKSMultiACBulkC.
 *
 * The non-linear variable for this Kernel is the order parameter \f$ \eta_i \f$.
 */
class HuhMultiACBulkBase : public HuhACBulk<Real>
{
public:
  static InputParameters validParams();

  HuhMultiACBulkBase(const InputParameters & parameters);

  virtual void initialSetup();

protected:
  /// name of order parameter that derivatives are taken wrt (needed to retrieve the derivative material properties)
  // VariableName _etai_name;

  /// index of order parameter that derivatives are taken wrt
  // unsigned int _etai_var;

  /// Value of the variable free energy function
  const MaterialProperty<Real> & _prop_Fi;

  std::vector<const MaterialProperty<Real> *> _prop_dFidarg;

  const MaterialProperty<Real> & _prop_si;


  /// Names of free energy functions for each phase \f$ F_j \f$
  std::vector<MaterialPropertyName> _Fk_names;
  unsigned int _num_k;

  /// Values of the free energy functions for each phase \f$ F_j \f$
  std::vector<const MaterialProperty<Real> *> _prop_Fk;

  /// Derivatives of the free energy functions (needed for off-diagonal Jacobians)
  std::vector<std::vector<const MaterialProperty<Real> *>> _prop_dFkdarg;

  /// step function names
  std::vector<MaterialPropertyName> _sk_names;

  /// values of the step function
  std::vector<const MaterialProperty<Real> *> _prop_sk;

  /// switching function names
  std::vector<MaterialPropertyName> _hk_names;

  /// Values of the switching functions for each phase \f$ h_k \f$
  std::vector<const MaterialProperty<Real> *> _prop_hk;

  // /// Derivatives of the switching functions wrt the order parameter for this kernel
  // std::vector<const MaterialProperty<Real> *> _prop_dhkdetai;

  /// Names of the mobilities used in this thing
  std::vector<MaterialPropertyName> _Mik_names;

  /// Values of the mobilities used in this thing
  std::vector<const MaterialProperty<Real> *> _prop_Mik;

  const MaterialProperty<Real> & _num_phases;


};
