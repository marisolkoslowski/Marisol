//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HuhMultiJunction.h"

registerMooseObject("PhaseFieldApp", HuhMultiJunction);

InputParameters
HuhMultiJunction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Gradient energy for Allen-Cahn Kernel with constant Mobility and "
                             "Interfacial parameter for a coupled order parameter variable.");
  params.addRequiredCoupledVar("etak_names", "Coupled variable that the Laplacian is taken of");
  params.addParam<MaterialPropertyName>("si_name", "The step function corresponding to the ith variable.");
  params.addParam<MaterialPropertyName>("Mik_names", "Mik", "The mobility pairs");
  params.addParam<MaterialPropertyName>("Epsil_names", "Eps_il", "The interfacial energies epsilon used with kernel");
  params.addParam<MaterialPropertyName>("Epskl_names", "Eps_kl", "Interfacial energies");
  params.addParam<MaterialPropertyName>("wil_names", "wil", "Barrier Energies");
  params.addParam<MaterialPropertyName>("wkl_names", "wkl", "Barrier Energies");
  params.addParam<MaterialPropertyName>("sk_names", "The step function corresponding to the kth eta.");
  params.addParam<MaterialPropertyName>("num_phases", "A material property specifying the number of phases");
  return params;
}

HuhMultiJunction::HuhMultiJunction(const InputParameters & parameters)
  : Kernel(parameters),
    _prop_si(getMaterialProperty<Real>("si_name")),
    _sk_names(getParam<std::vector<MaterialPropertyName>>("sk_names")),
    _Mik_names(getParam<std::vector<MaterialPropertyName>>("Mik_names")),
    _num_k(_Mik_names.size()),
    _num_extra(_Epsil_names.size()),
    _Epsil_names(getParam<std::vector<MaterialPropertyName>>("Epsil_names")),
    _Epskl_names(getParam<std::vector<MaterialPropertyName>>("Epskl_names")),
    _wil_names(getParam<std::vector<MaterialPropertyName>>("wil_names")),
    _wkl_names(getParam<std::vector<MaterialPropertyName>>("wkl_names")),
    _grad_etak(coupledGradients("etak_names")),
    _prop_etak(coupledValues("etak_names")),
    _etak_var(coupledIndices("etak_names")),
    _num_phases(getMaterialProperty<Real>("num_phases"))
{
  if (_num_k != _etak_var.size()){
    paramError("etak_names", "Need to pass in as many etak names as Mik names");
  }

  if (_num_extra != _Epskl_names.size()){
    paramError("Epskl_names", "Need to pass in as many Epskl names as Epsil names");
  }

  if (_num_extra != _wil_names.size()){
    paramError("Epskl_names", "Need to pass in as many wil names as Epsil names");
  }

  if (_num_extra != _wkl_names.size()){
    paramError("Epskl_names", "Need to pass in as many wkl names as Epsil names");
  }

  for (unsigned int n=0; n<_num_k; ++n){
    _prop_Mik[n] = &getMaterialPropertyByName<Real>(_Mik_names[n]);
    _prop_sk[n] = &getMaterialPropertyByName<Real>(_sk_names[n]);
  }


  for(unsigned int n=0; n<_num_extra; ++n){
    _prop_Epsil[n] = &getMaterialPropertyByName<Real>(_Epsil_names[n]);
    _prop_Epskl[n] = &getMaterialPropertyByName<Real>(_Epskl_names[n]);

    _prop_wil[n] = &getMaterialPropertyByName<Real>(_wil_names[n]);
    _prop_wkl[n] = &getMaterialPropertyByName<Real>(_wkl_names[n]);

  }
}


Real
HuhMultiJunction::computeQpResidual()
{
  //return _grad_v[_qp] * _kappa[_qp] * _L[_qp] * _grad_test[_i][_qp];
  Real sum = 0;
  unsigned int _num_l = _num_k - 1;

  for (unsigned int m = 0; m<_num_k; ++m){
    for(unsigned int n = 0; n<_num_k; ++n){

      // m functions as the "k" index
      // n functions as the "l" index
      unsigned int mn = m*_num_l+n;

      if(m == n){
        sum+=0;
      }

      else{
        sum += (*_prop_Mik[m])[_qp] * _prop_si[_qp] * (*_prop_sk[m])[_qp] * (*_prop_sk[n])[_qp] 
          *( (*_prop_wil[mn])[_qp] - (*_prop_wkl[mn])[_qp]) * (*_prop_etak[n])[_qp];


        sum += (*_prop_Mik[m])[_qp] * _prop_si[_qp] * (*_prop_sk[m])[_qp] * (*_prop_sk[n])[_qp]
          * ( (*_prop_Epsil[mn])[_qp] - (*_prop_Epskl[mn])[_qp]) *(*_grad_etak[n])[_qp] * _grad_test[_i][_qp];
      }


    }

  }
  
  return sum/_num_phases[_qp];
}
Real
HuhMultiJunction::computeQPJacobian(){

  Real Jacobian = 0;

  // so far I don't think that the Jacobian contributes anything.
  return Jacobian;
}

Real
HuhMultiJunction::computeQpOffDiagJacobian(unsigned int jvar)
{
  // if (jvar == _v_var)
  //   return _grad_phi[_j][_qp] * _kappa[_qp] * _L[_qp] * _grad_test[_i][_qp];

  return 0.0;
}
