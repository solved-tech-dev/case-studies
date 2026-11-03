// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 解得科技
/* CalculiX solver-entry declarations used by the PanguLU adapter.
 * These are application ABI declarations, not SPOOLES or PanguLU SDK code.
 * Include after CalculiX.h. This integration uses the default 32-bit ITG ABI.
 */
#ifndef CCX_PANGULU_ENTRY_H
#define CCX_PANGULU_ENTRY_H
#ifdef INTSIZE64
#error "This adapter requires the default 32-bit CalculiX integer ABI"
#endif
void spooles_factor(double*,double*,double*,double*,double*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*);
void spooles_solve(double*,ITG*);
void spooles_cleanup(void);
void spooles(double*,double*,double*,double*,double*,double*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*);
void spooles_factor_rad(double*,double*,double*,double*,double*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*,ITG*);
void spooles_solve_rad(double*,ITG*);
void spooles_cleanup_rad(void);
#endif
