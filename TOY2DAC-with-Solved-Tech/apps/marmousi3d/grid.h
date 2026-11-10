// SPDX-License-Identifier: AGPL-3.0-only
// Copyright 2026 Solved Tech
#pragma once

#ifndef FWI_SPACING
#error "Specify FWI_SPACING (85, 75 or 65)."
#endif
#if FWI_SPACING == 85
#define FWI_EX 99
#define FWI_EY 14
#define FWI_EZ 35
#elif FWI_SPACING == 75
#define FWI_EX 112
#define FWI_EY 16
#define FWI_EZ 40
#elif FWI_SPACING == 65
#define FWI_EX 129
#define FWI_EY 18
#define FWI_EZ 46
#else
#error "Unsupported FWI_SPACING."
#endif
