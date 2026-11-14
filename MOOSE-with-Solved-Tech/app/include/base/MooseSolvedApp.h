// SPDX-License-Identifier: AGPL-3.0-only
#pragma once
#include "MooseApp.h"
class MooseSolvedApp : public MooseApp {
public:
  static InputParameters validParams();
  MooseSolvedApp(const InputParameters &);
  static void registerApps();
  static void registerAll(Factory &, ActionFactory &, Syntax &);
};
