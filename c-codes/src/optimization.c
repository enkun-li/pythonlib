/* ===============================================
 * File Name: optimization.c
 * Author: ekli
 * Mail: ekli_091@mail.dlut.edu.cn  
 * Created Time: 2019-06-20 20:18:20
 * ===============================================
 */

#include "common.h"
#include "gapp.h"

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* The aim function or fitness function */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// Name  : aimfunc
// Input : *pars
// Output: fitness
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double aimfunc(const gsl_vector *v, void *pars)
{
    double sigf, lenf;
    int N = 31;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC_R18.txt";
    double chisq;

    sigf = gsl_vector_get(v, 0);
    lenf = gsl_vector_get(v, 1);
    
    initial_gapp(N, filename);
    setup_gapp(sigf, lenf);
    chisq = -loglikelihood();
    free_gapp();
    return chisq;
}

double aimfunc_SN(const gsl_vector *v, void *pars)
{
    double sigf, lenf;
    int N = 1048;
    //char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    //char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    double chisq;
    
    sigf = gsl_vector_get(v, 0);
    lenf = gsl_vector_get(v, 1);

    initial_gapp(N, filename);
    //initial_gapp_cov(filecov);
    setup_gapp(sigf, lenf);
    chisq = -loglikelihood();
    free_gapp();
    return chisq;
}


int multi_dim_mini(double(*func)(const gsl_vector *, void *))
{
    double par[] = {0.0};
    const gsl_multimin_fminimizer_type *T =
        gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    
    size_t iter = 0;
    int status;
    double size;
    
    /* Starting point */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, 130.0);
    gsl_vector_set (x, 1, 2.0);

    /* Set initial step sizes to 1 */
    ss = gsl_vector_alloc (2);
    gsl_vector_set_all (ss, 1.0);
    
    /* Initialize method and iterate */
    minex_func.n = 2;
    minex_func.f = (*func);
    minex_func.params = par;
  
    s = gsl_multimin_fminimizer_alloc (T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if (status)
            break;
        
        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);
        
        if (status == GSL_SUCCESS)  {
            printf ("converged to minimum at\n");
            printf ("%5ld %12.6f %12.6f f() = %20.8f size = %.8f\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->fval, size);
        }
        
        printf ("%5ld %10.3e %10.3e f() = %7.3f size = %.3f\n",
                iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->fval, size);
    } while (status == GSL_CONTINUE && iter < 100);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);
    
    return status;
}

int main(void)
{
    multi_dim_mini(&aimfunc);
    multi_dim_mini(&aimfunc_SN);
    return 0;
}
