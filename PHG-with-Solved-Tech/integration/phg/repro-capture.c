/* SPDX-License-Identifier: LGPL-2.1-or-later
 * Optional full native coefficient capture, outside native solver timing. */
#include "phg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
void phgReproCapture(DOF*u){
    const char*dir=getenv("PHG_REPRO_FIELDS");static unsigned id=0;static int busy=0;
    char file[4096];size_t n,i;FILE*f;int finite=1;
    if(!dir || !*dir || busy || !u || !u->data || !u->name || !u->g || !phgGridIsValid(u->g))return;
    if(strcmp(u->name,"Displacement") && strcmp(u->name,"u1_h") && strcmp(u->name,"u2_h"))return;
    n=DofGetDataCount(u);if(!n)return;busy=1;++id;
    for(i=0;i<n;++i)if(!isfinite((double)u->data[i]))finite=0;
    if(snprintf(file,sizeof(file),"%s/field-%u.f64",dir,id)>=(int)sizeof(file))phgError(1,"Capture path too long.\n");
    f=fopen(file,"wb");if(!f || fwrite(u->data,sizeof(FLOAT),n,f)!=n || fclose(f))phgError(1,"Native field capture failed.\n");
    snprintf(file,sizeof(file),"%s/field-%u.meta",dir,id);f=fopen(file,"w");if(!f)phgError(1,"Native metadata capture failed.\n");
    fprintf(f,"name=%s\ncount=%llu\nscalar_bytes=%u\nfinite=%d\nmesh_leaves=%lld\ndim=%d\n",u->name,(unsigned long long)n,(unsigned)sizeof(FLOAT),finite,(long long)u->g->nleaf_global,u->dim);
    fclose(f);
    /* Scalar examples already export VTK natively. Exporting their temporary
     * u_h interpolation fields from inside the exporter re-enters PHG VTK. */
    if(getenv("PHG_REPRO_VTK") && !strcmp(u->name,"Displacement")){
        snprintf(file,sizeof(file),"%s/field-%u.vtk",dir,id);phgExportVTK(u->g,file,u,NULL);
    }
    busy=0;
}
