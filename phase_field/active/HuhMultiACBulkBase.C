//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HuhMultiACBulkBase.h"

InputParameters
HuhMultiACBulkBase::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Multi-order parameter KKS model kernel for the Bulk Allen-Cahn. This "
                             "operates on one of the order parameters 'eta_i' as the non-linear "
                             "variable");
  params.addRequiredParam<MaterialPropertyName>("Fi_name", "Name of free energy corresponding to the given variable.");
  params.addRequiredParam<MaterialPropertyName>("si_name", "Name of the step function corresponding to the given variable.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "Fk_names", "List of free energies for each phase. Place in same order of least to greatest!");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "sk_names", "List of step functions for each phase. Place in order of least to greatest!");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "hk_names", "Switching Function Materials that provide h. Place in order of least to greatest!, where k != i.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>(
      "Mik_names", "Mobility pairs associated with each coupled phase. Place in order of least to greatest k.");
  params.addRequiredParam<MaterialPropertyName>(
      "num_phases", "The number of phases located at the current point");
  // params.addRequiredCoupledVar("eta_i",
  //                              "Order parameter that derivatives are taken with respect to");
  return params;
}

HuhMultiACBulkBase::HuhMultiACBulkBase(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    // _etai_name(getVar("eta_i", 0)->name()),
    // _etai_var(coupled("eta_i", 0)),
    _prop_Fi(getMaterialProperty<Real>("Fi_name")),
    _prop_si(getMaterialProperty<Real>("si_name")),
    _Fk_names(getParam<std::vector<MaterialPropertyName>>("Fk_names")),
    _num_k(_Fk_names.size()),
    _prop_Fk(_num_k),
    _prop_dFkdarg(_num_k),
    _sk_names(getParam<std::vector<MaterialPropertyName>>("sk_names")),
    _prop_sk(num_k),
    _hk_names(getParam<std::vector<MaterialPropertyName>>("hk_names")),
    _prop_hk(_num_k),
    _prop_dhkdetai(_num_k),
    _Mik_names(getParam<std::vector<MaterialPropertyName>>("Mik_names")),
    _num_phases(getMaterialProperty<Real>("num_phases")),
    _prop_Mik(_num_k)
    
{
  // check passed in parameter vectors
  if (_num_k != _hk_names.size())
    paramError("hk_names", "Need to pass in as many hk_names as Fk_names-1");

  }
  // reserve space and set phase material properties
  for (unsigned int n = 0; n < _num_k; ++n)
  {
    // get phase free energy
    _prop_Fk[n] = &getMaterialPropertyByName<Real>(_Fk_names[n]);
    _prop_dFkdarg[n].resize(_n_args);

    // get step function value
    _prop_sk[n] = &getMaterialPropertyByName<Real>(_sk_names[n]);

    _prop_Mik[n] = &getMaterialPropertyByName<Real>(_Mik_names[n]);

    for (unsigned int i = 0; i < _n_args; ++i)
    {
      // Get derivatives of all Fj wrt all coupled variables
      _prop_dFkdarg[n][i] = &getMaterialPropertyDerivative<Real>(_Fk_names[n], i);
    }
  }

  for (unsigned int i = 0; i < _n_args)
  {
    _prop_dFidarg = &getMaterialPropertyDerivative<Real>("Fi_name", i);
  }

  for(unsigned int n = 0; n< _num_k; ++n)
  {
    // get switching function and derivatives wrt eta_i, the nonlinear variable
    _prop_hk[n] = &getMaterialPropertyByName<Real>(_hk_names[n]);
    // _prop_dhkdetai[n] = &getMaterialPropertyDerivative<Real>(_hk_names[n], _etai_name);

  }
}

void
HuhMultiACBulkBase::initialSetup()
{
  ACBulk<Real>::initialSetup();

  // Test. If this line breaks, delete it.
  validateNonlinearCoupling<Real>("Fi_name");

  for (unsigned int n = 0; n < _num_k; ++n)
  {
    validateNonlinearCoupling<Real>(_Fk_names[n]);
  }
}
