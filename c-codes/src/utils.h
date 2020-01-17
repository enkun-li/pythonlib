/* ===============================================
 * File Name: utils.h
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 10:35:23
 * =============================================== 
 */

#include <stdio.h>

#ifndef __UTILS__
#define __UTILS__

/* gaus leg integrate */
int gauleg(double x1, double x2, double x[], double w[]);
double qgausleg(double (*func)(double, double *), double a, double b, double *args );
double qgausleg2d_fix(double (*func)(double, double, double *), double a, double b, double c, double d, double *args );

/* */
double qgaus(double (*func)(double, double *), double a, double b, double *args);
double qgaus2d_fix(double (*func)(double, double, double *), double a, double b, double c, double d, double *args);

/* inverse matrix */
int gsl_matrix_inv(gsl_matrix *M);

#endif /* __UTILS__ */
