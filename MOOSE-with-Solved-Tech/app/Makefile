# SPDX-License-Identifier: AGPL-3.0-only
MOOSE_DIR ?= ../../moose
FRAMEWORK_DIR := $(MOOSE_DIR)/framework
MODULE_DIR := $(MOOSE_DIR)/modules
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk
SOLID_MECHANICS := yes
HEAT_TRANSFER := yes
INCLUDE_COMBINED := no
SKIP_LOADER := yes
USE_TEST_LIBS := no
include $(MODULE_DIR)/modules.mk
APPLICATION_DIR := $(CURDIR)
APPLICATION_NAME := moosesolved
BUILD_EXEC := yes
GEN_REVISION := no
BUILD_TEST_OBJECTS_LIB := no
include $(FRAMEWORK_DIR)/app.mk
