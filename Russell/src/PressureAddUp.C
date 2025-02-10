
#include "PressureAddUp.h"
#include "libmesh/quadrature.h"

registerMooseObject("MooseApp", PressureAddUp);

InputParameters
PressureAddUp::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("mat");
  params.addRequiredParam<std::string>("property_name",     "The property name to declare");
  params.addParam<std::vector<std::string>>("sum_materials", "Base name of the parsed sum material property");
  params.addRequiredCoupledVar("values", "Vector of values to sum");

  return params;
}

PressureAddUp::PressureAddUp(const InputParameters & parameters)
  : Material(parameters),
    _property_name(getParam<std::string>("property_name")),
    _property(declareProperty<Real>(_property_name)),
    _sum_materials(getParam<std::vector<std::string>>("sum_materials")),
    _num_materials(_sum_materials.size()),
    _n_values(coupledComponents("values"))
{
  // we need at least one material in the sum
  if (_num_materials == 0)
    mooseError("Please supply the pressures ", name());

  // reserve space for summand material properties
  _summand_F.resize(_num_materials);

  for (unsigned int i = 0; i < _n_values; i++)
    _values.push_back(&coupledValue("values", i));

  for (unsigned int n = 0; n < _num_materials; ++n)
    _summand_F[n] = &getMaterialPropertyByName<Real>(_sum_materials[n]);
}

void
PressureAddUp::initQpStatefulProperties()
{
}

void
PressureAddUp::computeQpProperties()
{
  //_property[_qp] = 0;
  _property[_qp] = (*(_values[0]))[_qp] * (*_summand_F[0])[_qp];
  for (unsigned int n = 1; n < _num_materials; ++n)
        _property[_qp] += (*(_values[n]))[_qp] * (*_summand_F[n])[_qp];
}
























///-------------------------------------------------------------------------------------------------------------------------------------------------------------------------