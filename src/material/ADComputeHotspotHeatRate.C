#include "ADComputeHotspotHeatRate.h"

registerMooseObject("catApp", ADComputeHotspotHeatRate);

InputParameters
ADComputeHotspotHeatRate::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute heat rate for hotspot assignment");
    params.addRequiredParam<MaterialPropertyName>("target_name", "target temperature distribution name");
    params.addRequiredParam<MaterialPropertyName>("tau_name", "induction time distribution name");
    params.addRequiredParam<bool>("use_lump", "whether to compute lumped terms or not");
    return params;
}

ADComputeHotspotHeatRate::ADComputeHotspotHeatRate(const InputParameters & parameters)
  : Material(parameters),
    
    _rho(getADMaterialProperty<Real>("density")),
    _cv(getADMaterialProperty<Real>("specific_heat")),
    _target_name(getParam<MaterialPropertyName>("target_name")),
    _target(getADMaterialPropertyByName<Real>(_target_name)),
    _tau_name(getParam<MaterialPropertyName>("tau_name")),
    _tau(getADMaterialProperty<Real>(_tau_name)),

    _use_PK2(getParam<bool>("use_PK2")),
    //declare properties
    _q_elastic(declareADProperty<Real>("q_elastic")),
    _use_lump(getParam<bool>("use_lump"))
{   
}

void
ADComputeHotspotHeatRate::computeQpProperties()
{
    RankTwoTensor I2(RankTwoTensor::initIdentity);

    //component contribution from volumetric compression
    ADReal q_pressure;
    RankTwoTensor Ce = _Fe[_qp].transpose() * _Fe[_qp];

    //if use PK2, use work conjugate C.inverse()
    if (_use_PK2){
        q_pressure = - std::max(_T[_qp] * _dP_dT[_qp] * (Ce.inverse().doubleContraction(_Ee_dot[_qp])), 0.);
    }else{
        q_pressure = - std::max(_T[_qp] * _dP_dT[_qp] * _Ee_dot[_qp].trace(), 0.0);
    }

    //add artificial viscosity factor

    ADReal q_av;
    q_av = _beta_av * _pressure_av[_qp] * _Ee_dot[_qp].trace();

    ADReal q_tot = q_pressure + q_av;

    //activation for MISTERnet simulations
    if(_dirac_switch_react[_qp] > _thr_activation){
        q_tot *= 1.; //keep while activated
    }else{
        q_tot *= 0.; //set to zero before activation
    }

    _q_elastic[_qp] = q_tot;
}