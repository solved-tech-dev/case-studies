/* SPDX-License-Identifier: LGPL-2.1-or-later
 * PHG-Solid's unmodified 2011 source uses the legacy three-argument normal API.
 * This build-only C++ overload delegates to PHG 0.9.7's four-argument API. */
#pragma once
#include <cstring>
#include "phg.h"
inline const FLOAT *phgGeomGetFaceOutNormal(GRID *g, ELEMENT *e, int face) {
    static thread_local FLOAT normal[3];
    phgGeomGetFaceOutNormal(g,e,face,normal);
    return normal;
}
