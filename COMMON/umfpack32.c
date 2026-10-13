#include "umfpack.h"
#include <stdint.h>
#include <stdlib.h>

/* Experimental real, 32-bit-index path; Numeric remains an opaque pointer. */
void *neo2_umf32_factor(int32_t n, int32_t nz, const int32_t *pcol,
                       const int32_t *irow, const double *values,
                       int32_t *status)
{
    void *symbolic = NULL, *numeric = NULL;
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    *status = UMFPACK_ERROR_invalid_matrix;
    if (n <= 0 || n == INT32_MAX || nz <= 0 || nz == INT32_MAX) return NULL;
    int32_t *ap = malloc(((size_t)n + 1) * sizeof(*ap));
    int32_t *ai = malloc((size_t)nz * sizeof(*ai));
    if (ap == NULL || ai == NULL) {
        free(ap);
        free(ai);
        *status = UMFPACK_ERROR_out_of_memory;
        return NULL;
    }
    for (int32_t j = 0; j <= n; ++j) ap[j] = pcol[j] - 1;
    for (int32_t j = 0; j < nz; ++j) ai[j] = irow[j] - 1;
    umfpack_di_defaults(control);
    *status = umfpack_di_symbolic(n, n, ap, ai, values, &symbolic,
                                 control, info);
    if (*status >= 0) {
        *status = umfpack_di_numeric(ap, ai, values, symbolic, &numeric,
                                    control, info);
    }
    umfpack_di_free_symbolic(&symbolic);
    free(ap);
    free(ai);
    if (*status < 0) umfpack_di_free_numeric(&numeric);
    return numeric;
}

int32_t neo2_umf32_solve(int32_t n, void *numeric, double *b)
{
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    double *x = malloc((size_t)n * sizeof(*x));
    if (x == NULL) return UMFPACK_ERROR_out_of_memory;
    umfpack_di_defaults(control);
    control[UMFPACK_IRSTEP] = 0;
    int status = umfpack_di_solve(UMFPACK_A, NULL, NULL, NULL, x, b,
                                 numeric, control, info);
    if (status >= 0) {
        for (int32_t j = 0; j < n; ++j) b[j] = x[j];
    }
    free(x);
    return status;
}

void neo2_umf32_free(void *numeric)
{
    umfpack_di_free_numeric(&numeric);
}
