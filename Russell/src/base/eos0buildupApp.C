#include "eos0buildupApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
eos0buildupApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

eos0buildupApp::eos0buildupApp(InputParameters parameters) : MooseApp(parameters)
{
  eos0buildupApp::registerAll(_factory, _action_factory, _syntax);
}

eos0buildupApp::~eos0buildupApp() {}

void
eos0buildupApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<eos0buildupApp>(f, af, s);
  Registry::registerObjectsTo(f, {"eos0buildupApp"});
  Registry::registerActionsTo(af, {"eos0buildupApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
eos0buildupApp::registerApps()
{
  registerApp(eos0buildupApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
eos0buildupApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  eos0buildupApp::registerAll(f, af, s);
}
extern "C" void
eos0buildupApp__registerApps()
{
  eos0buildupApp::registerApps();
}
