// SPDX-License-Identifier: AGPL-3.0-only
#include "AppFactory.h"
#include "HeatTransferApp.h"
#include "MooseSolvedApp.h"
#include "MooseSyntax.h"
#include "RayTracingApp.h"
#include "ShiftedBoundaryMethodApp.h"
#include "SolidMechanicsApp.h"
registerKnownLabel("MooseSolvedApp");
InputParameters MooseSolvedApp::validParams() {
  auto p = MooseApp::validParams();
  p.set<bool>("automatic_automatic_scaling") = false;
  p.set<bool>("use_legacy_material_output") = false;
  p.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return p;
}
MooseSolvedApp::MooseSolvedApp(const InputParameters &p) : MooseApp(p) {
  registerAll(_factory, _action_factory, _syntax);
}
void MooseSolvedApp::registerApps() {
  registerApp(MooseSolvedApp);
  SolidMechanicsApp::registerApps();
  HeatTransferApp::registerApps();
  RayTracingApp::registerApps();
  ShiftedBoundaryMethodApp::registerApps();
}
void MooseSolvedApp::registerAll(Factory &f, ActionFactory &a, Syntax &s) {
  Registry::registerObjectsTo(f, {"MooseSolvedApp"});
  Registry::registerActionsTo(a, {"MooseSolvedApp"});
  SolidMechanicsApp::registerAll(f, a, s);
  HeatTransferApp::registerAll(f, a, s);
  RayTracingApp::registerAll(f, a, s);
  ShiftedBoundaryMethodApp::registerAll(f, a, s);
}
extern "C" void MooseSolvedApp__registerAll(Factory &f, ActionFactory &a, Syntax &s) {
  MooseSolvedApp::registerAll(f, a, s);
}
extern "C" void MooseSolvedApp__registerApps() { MooseSolvedApp::registerApps(); }
