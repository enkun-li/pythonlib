/* ===============================================
 * File Name: pso.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-14 20:57:23
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"

/* ================================= */
/*             PSO pars              */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// #define w 1.0
#define c1 1.49445
#define c2 1.49445
#define maxgen 1000 // iter times
#define sizepop 20 // size of pop
#define dim 2 // pars number

// iter pars
double Xmax[dim], Xmin[dim]; // max/min position
double Vmax[dim], Vmin[dim]; // max/min velocity
double X[sizepop][dim];      // position
double V[sizepop][dim];      // velocity
double pbest[sizepop][dim];  // personal best
double gbest[dim];           // global best
double p_fit[sizepop];       // personal fit
double g_fit = 1.0e30;                // global fit
double genbest[maxgen][dim]; // best pars per gen
double result[maxgen];       // best fit per gen


/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/* The aim function or fitness function */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// Name  : aimfunc
// Input : *pars
// Output: fitness
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
double aimfunc_OHD(double *pars)
{
    //double x=*pars;
    //double y=*(pars+1);

    //return x*x + y*y -10*(cos(2*_PI_*x) + cos(2*_PI_*y));
    double sigf = *pars, lenf=*(pars+1);
    int N = 31;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC_R18.txt";
    //double sigf=170.75260325, lenf=2.61947963;
    double chisq;

    initial_gapp(N, filename);
    setup_gapp(sigf, lenf);
    chisq = -loglikelihood();
    free_gapp();
    return chisq;
}

double aimfunc_Pantheon(double *pars)
{
    //double x=*pars;
    //double y=*(pars+1);

    //return x*x + y*y -10*(cos(2*_PI_*x) + cos(2*_PI_*y));
    double sigf = *pars, lenf=*(pars+1);
    //int N = 1048;
    //char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    //char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    int N = 1049;
    char * filename = "/home/ekli/myworks/cosmodata/sn_dat_to_dz.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_cov_to_dz.txt";

    double chisq;

    initial_gapp(N, filename);
    initial_gapp_cov(N, filecov);
    setup_gapp(sigf, lenf);
    chisq = -loglikelihood();
    free_gapp();
    return chisq;
}

//double aimfunc(double * pars)
//{
//    return aimfunc_SN(&pars[0]);
//}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*      Initial the population          */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
// Name  : initial_pop
// Input : *maxpars, *minpars
// Output: 
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int initial_pop(double *maxpars, double *minpars)
{
    int i, j, k;
    double r;

    for (i=0; i<dim; i++) {
        Xmax[i] = *(maxpars+i);
        Xmin[i] = *(minpars+i);
        Vmax[i] = 0.5 * (Xmax[i] - Xmin[i]);
        Vmin[i] = - Vmax[i];
    }

    for (j=0; j<sizepop; j++) {
        for (i=0; i<dim; i++) {
            r = (double)rand()/RAND_MAX;
            X[j][i] = Xmin[i] + (Xmax[i] - Xmin[i]) * r;
            r = (double)rand()/RAND_MAX;
            V[j][i] = Vmin[i] + (Vmax[i] - Vmin[i]) * r;
            pbest[j][i] = X[j][i]; //pbest position
        }
        p_fit[j] = aimfunc(&X[j][0]);
        if(p_fit[j] < g_fit) {
            g_fit = p_fit[j];
            for (k=0; k<dim; k++) {
                gbest[k] = X[j][k];
            }
        }
    }
    return _SUCCESS_;
}


/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*   pso update position and velocity   */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Name : pso_update                    */
/* Input: w                             */
/* Output:                              */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int pso_update(double w)
{
    int i, j, k;
    double r1, r2;
    double temp;
    
    for (j=0; j<sizepop; j++) {
        for (i=0; i<dim; i++) {
            // update velocity
            r1 = (double)rand()/RAND_MAX;
            r2 = (double)rand()/RAND_MAX;

            V[j][i] = w*V[j][i] +c1*r1*(pbest[j][i] - X[j][i]) 
                + c2*r2*(gbest[i] -X[j][i]);
            // whether big/less than max/min
            if(V[j][i] > Vmax[i])
                V[j][i] = Vmax[i];
            else if(V[j][i] < Vmin[i])
                V[j][i] = Vmin[i];

            // update position
            X[j][i] = X[j][i] + V[j][i];
            // whether 
            if(X[j][i] > Xmax[i])
                X[j][i] = Xmax[i];
            else if(X[j][i] < Xmin[i])
                X[j][i] = Xmin[i];
        }
        temp = aimfunc(X[j]);
        // update pbest
        if(temp < p_fit[j]) {
            p_fit[j] = temp;
            for (k=0; k<dim; k++) {
                pbest[j][k] = X[j][k];
            }
        }
        // update gbest
        if(temp < g_fit) {
            g_fit = temp;
            for (k=0; k<dim; k++) {
                gbest[k] = X[j][k];
            }
        }
    }
    return 0;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/*               pso iter               */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* Name:  pso_iter                      */
/* Input: w                             */
/* Output:                              */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int pso_iter(void)
{
    int i, k; //, j;
    double w;

    w = 0.72;
    //for (j=0; j< 10; j++) {
    //    w = 1.5 - (1.5-0.3)/10*j;
    for (k=0; k<maxgen; k++) {
        w = 0.9 - (0.9-0.2)*(k/maxgen)*(k/maxgen);
        pso_update(w);
        printf("%5d x = (", k);
        for (i=0; i<dim; i++) {
            genbest[k][i] = gbest[i];
            printf("%20.8f", gbest[i]);
        }
        printf("),   ");
        result[k] = g_fit;
        printf("best = %20.8f\n", g_fit);
    }
    //}
    return 0;
}

/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// main func
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(void)
{
    clock_t start,finish; //程序开始和结束时间
    int i;
    double maxpars[] = {300.0, 10.0};
    double minpars[] = {0.0, 0.0};

    start = clock(); //开始计时
    srand((unsigned)time(NULL)); // 初始化随机数种子

    initial_pop(maxpars, minpars);
    pso_iter();

    printf("迭代了%d次，最优值为:%lf.\n",maxgen,g_fit);
    printf("取到最优值的位置(");
    for (i=0; i<dim; i++) {
        printf("%lf,", gbest[i]);
    }
    printf(")\n");
    finish = clock(); //结束时间
    double duration = (double)(finish - start)/CLOCKS_PER_SEC; // 程序运行时间
    printf("程序运行耗时:%lf\n",duration);
    return 0;
}

