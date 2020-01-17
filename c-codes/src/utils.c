/* ===============================================
 * File Name: unitils.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-05-16 16:31:13
 * =============================================== 
 */

#include "common.h"

#define KMAX 50
#define EPS 1.0e-11

/**************************************************/
/* Using Gauss method to integrate function */
int NPOINT=5;

int gauleg(double x1, double x2, double x[], double w[])
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(NPOINT+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(_PI_*(i-0.25)/(NPOINT+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=NPOINT;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=NPOINT*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[NPOINT+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[NPOINT+1-i]=w[i];
	}
    return _SUCCESS_;
}

double qgausleg(double (*func)(double, double *),
        double a, double b,
        double *args )
{
    double x[NPOINT+1];
    double w[NPOINT+1];
    double s;
    int status, i;

    status = gauleg(a, b, x, w);
    if(status == _FAILURE_)
        printf("Failed initial gauleg.\n");

    s = 0;
    for (i=1; i<=NPOINT; i++){
        s += w[i] * (*func)(x[i], &args[0]);
    }
    return s;
}

double qgausleg2d_fix(double (*func)(double, double, double *),
        double a, double b, double c, double d,
        double *args )
{
    double xa[NPOINT+1], xc[NPOINT+1];
    double wa[NPOINT+1], wc[NPOINT+1];
    double s;
    int status, i, j;

    status = gauleg(a, b, xa, wa);
    if(status == _FAILURE_)
        printf("Failed initial gauleg.\n");
    status = gauleg(c, d, xc, wc);
    if(status == _FAILURE_)
        printf("Failed initial gauleg.\n");

    s = 0;
    for (i=1; i<=NPOINT; i++){
        for (j=1; j<=NPOINT; j++){
            s += wa[i] * wc[j] * (*func)(xa[i], xc[j], &args[0]);
        }
    }
    return s;
}

double qgaus(double (*func)(double, double *), 
        double a, double b, 
        double *args)
{
	int j;
	double xr,xm,dx,s;
    static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j] * ( (*func)(xm+dx, &args[0]) 
                + (*func)(xm-dx, &args[0]) );
	}
	return s *= xr;
}

double qgaus2d_fix(double (*func)(double, double, double *), 
        double a, double b, double c, double d, 
        double *args)
{
	int i,j;
	double xr,xm,dx,s;
    double yr,ym,dy;
    static double x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static double w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
    ym=0.5*(d+c);
    yr=0.5*(d-c);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
        for(i=1;i<=5;i++){
            dy = yr*x[i];
            s+= w[j]* w[i] * ( (*func)(xm+dx, ym+dy, &args[0]) 
                    + (*func)(xm-dx, ym-dy, &args[0]) );
        }
	}
    s = s*xr*yr;
	return s;
}

// ================================
// double qgaus(double (*func)(double, double *), double a, double b, double *args)
// {
// 	int i, j, k;
// 	double xr,xm,dx,s, s0;
//     double ai, bi;
// 	static double x[]={0.0,0.1488743389,0.4333953941,
// 		0.6794095682,0.8650633666,0.9739065285};
// 	static double w[]={0.0,0.2955242247,0.2692667193,
// 		0.2190863625,0.1494513491,0.0666713443};
// 
//     s0 = 1e30;
//     for (k=1; k< KMAX; k++){
//         s = 0;
//         for (i=0; i<k; i++){
//             ai = a + (b-a)/pow(2, k)*i;
//             bi = ai + (b-a)/k;
//             xm=0.5*(bi+ai);
//             xr=0.5*(bi-ai);
//             for (j=1;j<=5;j++) {
//                 dx=xr*x[j];
//                 s += w[j]*((*func)(xm+dx, &args[0])+(*func)(xm-dx, &args[0]));
//             }
//         }
//         s *= xr;
//         if(fabs((s - s0)/s0) <= EPS || (s == 0.0 || s0 == 0.0) )
//             return s;
//         else
//             s0 = s;
//     }
// 	// return s *= xr;
//     printf("Failed to get the precision.\n");
//     return 0;
// }
// 
// double qgaus2d_fix(double (*func)(double, double, double *), double a, double b, double c, double d, double *args)
// {
// 	int i,j, k, l;
// 	double xr,xm,dx,s;
//     double ai, bi, ci, di, s0;
//     double yr,ym,dy;
// 	static double x[]={0.0,0.1488743389,0.4333953941,
// 		0.6794095682,0.8650633666,0.9739065285};
// 	static double w[]={0.0,0.2955242247,0.2692667193,
// 		0.2190863625,0.1494513491,0.0666713443};
// 
//     s0 = 1e30;
//     for(k=1; k<KMAX; k++) {
//         s = 0;
//         for (l=0; l<k; l++){
//             ai = a + (b-a)/pow(2, k)*l;
//             bi = ai + (b-a)/k;
//             ci = c + (d-c)/k*l;
//             di = ci + (d-c)/k;
//         
//             xm=0.5*(bi+ai);
//             xr=0.5*(bi-ai);
//             ym=0.5*(di+ci);
//             yr=0.5*(di-ci);
//             for (j=1;j<=5;j++) {
//                 dx=xr*x[j];
//                 for(i=1;i<=5;i++){
//                     dy = yr*x[i];
//                     s+= w[j]* w[i] * ((*func)(xm+dx, ym+dy, &args[0]) 
//                             +(*func)(xm-dx, ym-dy, &args[0]));
//                 }
//             }
//         }
//         s = s*xr*yr;
//         if((s - s0)/s0 <= EPS || (s == 0 && s0 == 0) )
//             return s;
//         else
//             s0 = s;
//     }
//     // s = s*xr*yr;
// 	// return s;
//     printf("Failed to get the precision.\n");
//     return 0;
// }

/* ============================== */
// double func(double x, double * args)
// {
//     return *(args) * sin(*(args+1)*x);
// }
// 
// int test_qgaus(void)
// {
//     int i;
//     double x, ff, ff2, rr;
//     double pars[] = {1.0, 1.0};
// 
//     for(i=0; i<51; i++){
//         x = PI*2.0/50*i;
//         ff = qgaus(&func, 0.0, x, &pars[0]);
//         ff2 = qgausleg(&func, 0.0, x, &pars[0]); 
//         rr = 1-cos(x);
//         printf("int sinx dx = %12.9f  real = %12.9f  Delta = %12.9e  ", ff, rr, ff-rr );
//         printf("leg = %12.9f  Delta = %12.9e\n", ff2, ff2-rr);
//     }
//     return 0;
// }
// 
// /* ============================== */
// double func_2d(double x, double y, double * args)
// {
//     return pow(x, *(args)) + pow(y, *(args+1)) +x*y;
// }
// 
// int test_qgaus2d_fix(void)
// {
//     int i;
//     double x, ff, ff2, rr;
//     double pars[] = {2.0, 2.0};
// 
//     for(i=0; i<51; i++){
//         x = 2*PI/50*i;
//         ff = qgaus2d_fix(&func_2d, 0.0, x, 0.0, x, &pars[0]);
//         rr = 11.0/12*pow(x, 4.0);
//         // rr = pow(x, 3); // 2.0/3.0*pow(x,4);
//         printf("int sinx dx = %15.9f  real = %15.9f  Delta = %15.9e  ", ff, rr, ff-rr );
//         ff2 = qgausleg2d_fix(&func_2d, 0.0, x, 0.0, x, &pars[0]);
//         printf("leg = %15.9f  Delta = %15.9e\n", ff2, ff2-rr);
//     }
//     return 0;
// }
// 
// /* ============================== */
// int main(void)
// {
//     int status;
// 
//     status = test_qgaus();
//     status = test_qgaus2d_fix();
// 
//     if(status != 0)
//         printf("Failed");
// 
//     return 0;
// }

/* matrix inverse */
int gsl_matrix_inv(gsl_matrix *M)
{
    size_t n = M->size1;
    // size_t m = M->size2;
    
    gsl_matrix *temp1 = gsl_matrix_calloc(n,n);    
    gsl_matrix_memcpy(temp1, M);

    gsl_permutation *p = gsl_permutation_calloc(n);

    int sign=0;

    gsl_linalg_LU_decomp(temp1, p, &sign);
    gsl_matrix * inverse = gsl_matrix_calloc(n,n);

    gsl_linalg_LU_invert(temp1, p, inverse);
    gsl_matrix_memcpy(M, inverse);

    gsl_permutation_free(p);
    gsl_matrix_free(temp1);
    gsl_matrix_free(inverse);
    
    return _SUCCESS_;
}
