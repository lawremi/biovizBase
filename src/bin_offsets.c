#include <stdint.h>

#include "bin_offsets.h"

static int _scan_bin_chunk_count(Rbyte **b) {
  *b += 4; // skip the bin id
  int32_t count = *(int32_t *)(*b);
  *b += 4 + count * 16; // skip the actual counts
  return count;
}

static void _scan_bin_offsets(Rbyte **b, double **m) {
  int32_t bin = *(int32_t *)(*b); *b += 4;
  int32_t nchunks = *(int32_t *)(*b); *b += 4;
  for (int i = 0; i < nchunks; i++, *m += 5) {
    int64_t start = *(int64_t *)(*b); *b += 8;
    int64_t end = *(int64_t *)(*b); *b += 8;
    (*m)[0] = bin;
    (*m)[1] = start >> 16;
    (*m)[2] = end >> 16;
    (*m)[3] = (uint16_t)start;
    (*m)[4] = (uint16_t)end;
  }
}

static SEXP _scan_bin_offsets_seq(Rbyte **b) {
  int32_t nbin = *(int32_t *)(*b); *b += 4;
  
  Rbyte *tmp_bytes = *b;
  int nchunks = 0;
  for (int i = 0; i < nbin; i++) {
    nchunks += _scan_bin_chunk_count(&tmp_bytes);
  }
  
  // need real values, because 32 bits is generally not enough
  SEXP ans;
  PROTECT(ans = allocMatrix(REALSXP, 5, nchunks));
  double *m = REAL(ans);
  for (int i = 0; i < nbin; i++) {
    _scan_bin_offsets(b, &m);
  }

  // skip past the interval index
  *b += 4 + *(int32_t *)(*b) * 8;

  UNPROTECT(1);
  return ans;
}

SEXP scan_bam_bin_offsets(SEXP bytes) {
  SEXP ans;
  
  if (!IS_RAW(bytes))
    Rf_error("'bytes' must be a raw vector");
  Rbyte *b = RAW(bytes);

  if (strncmp(b, "BAI\1", 4))
    Rf_error("wrong magic number");
  b += 4;
  
  int32_t nref = *(int32_t *)b; b += 4;
  
  PROTECT(ans = NEW_LIST(nref));
  for (int i = 0; i < nref; i++)
    SET_VECTOR_ELT(ans, i, _scan_bin_offsets_seq(&b));

  UNPROTECT(1);
  return ans;
}
