// SPDX-License-Identifier: AGPL-3.0-only
#include "MooseMain.h"
#include "MooseSolvedApp.h"
int main(int argc, char *argv[]) { return Moose::main<MooseSolvedApp>(argc, argv); }
