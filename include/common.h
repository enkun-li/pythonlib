/* ===============================================
 * File Name: common.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-10-21 17:09:56
 * =============================================== 
 */

#ifndef __COMMON__
#define __COMMON__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>

#define _TRUE_ 1 /* integer associated to true statement */
#define _FALSE_ 0 /* integer associated to false statement */

#define _SUCCESS_ 0 /* integer returned after successful call of a function */
#define _FAILURE_ 1 /* integer returned after failure call of a function */

#define _PI_ 3.1415926535897932384626433832795e0 /* the number of pi*/

#define _C_ 299792458.0e0 /* speed of light in m/s */

typedef struct {
    char* info;
    int n;
    double *x;
    double *y;
    double *z;
} Obs_Data;

#endif /* __COMMON__ */

