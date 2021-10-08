//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HuhACInterface.h"

registerMooseObject("PhaseFieldApp", HuhACInterface);

InputParameters
HuhACInterface::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Gradient energy for Allen-Cahn Kernel with constant Mobility and "
                             "Interfacial parameter for a coupled order parameter variable.");
  params.addRequiredCoupledVar("etak_names", "Coupled variable that the Laplacian is taken of");
  params.addMaterialProperyName("si_name", "The step function corresponding to the ith variable.")
  params.addParam<MaterialPropertyName>("Mik_names", "Mik", "The mobility pairs");
  params.addParam<MaterialPropertyName>("Epsik_names", "Eps_ik", "The interfacial energies epsilon used with kernel");
  params.addParam<MaterialPropertyName>("sk_names", "The step function corresponding to the kth eta.");
  params.addParam<MaterialPropertyName>("num_phases", "A material property specifying the number of phases");
  return params;
}

HuhACInterface::HuhACInterface(const InputParameters & parameters)
  : Kernel(parameters),
    _Mik_names(getParam<std::vector<MaterialPropertyName>>("Mik_names")),
    _num_k(_Mik_names.size())
    _Epsik_names(getParam<std::vector<MaterialPropertyName>>("Epsik_names")),
    _etak_var(coupledIndices("etak_names")),
    _grad_etak(coupledGradients("etak_names")),
    _sk_names(getParam<std::vector<MaterialPropertyName>>("sk_names"),)
    _prop_si(getMaterialProperty<Real>("si_name")),
    _num_phases(getMaterialProperty<Real>("num_phases"))
{
  if (_num_k != _Epsik_names.size()){
    paramError("Epsik_names", "Need to pass in as many Epsik names as Mik names");
  }
  if (_num_k != _sk_names.size()){
    paramError("sk_names", "Need to pass in as many sk names as Mik names");
  }
  for (unsigned int n=0; n<_num_k; ++n){
    _prop_Mik[n] = &getMaterialPropertyByName<Real>(_Mik_names[n]);
    _prop_sk[n] = &getMaterialProperyByName<Real>(_sk_names[n]);
    _prop_Epsik[n] = &getMaterialPropertyByName<Real>(_Epsik_names[n]);
  }
}

Real
HuhACInterface::computeQpResidual()
{
  //return _grad_v[_qp] * _kappa[_qp] * _L[_qp] * _grad_test[_i][_qp];
  Real sum = 0;
  for (unsigned int n = 0; n<_num_k; ++n){

    // NOTE: This is an addition!
    sum += (*_prop_sk[n])[_qp] * _prop_si[_qp] * (*_prop_Mik[n])[_qp]* (*_prop_Epsik[n])[_qp]
      *(*_grad_etak[n])[_qp] * _grad_test[i][_qp];

    // NOTE: This is a subtraction!
    sum -= (*_prop_sk[n])[_qp] * _prop_si[_qp] * (*_prop_Mik[n])[_qp]* (*_prop_Epsik[n])[_qp]
      *_grad_u[_qp] * _grad_test[i][_qp];
  }
  
  return sum/_num_phases[_qp];
}

Real
HuhACInterface::computeQpOffDiagJacobian(unsigned int jvar)
{
  // if (jvar == _v_var)
  //   return _grad_phi[_j][_qp] * _kappa[_qp] * _L[_qp] * _grad_test[_i][_qp];

  return 0.0;
}
