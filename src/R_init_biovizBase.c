#include <R_ext/Rdynload.h>
#include "bin_offsets.h"

static const R_CallMethodDef callMethods[] = {
  {"scan_bam_bin_offsets", (DL_FUNC) & scan_bam_bin_offsets, 1},
  {NULL, NULL, 0}
};

void R_init_biovizBase(DllInfo * info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
