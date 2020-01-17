/* ===============================================
 * File Name: test_gapp.c
 * Author: ekli
 * mail: ekli_091@mail.dlut.edu.cn
 * Created Time: 2019-06-10 13:11:37
 * =============================================== 
 */

#include "common.h"
#include "gapp.h"
#include "gapp_rec.h"
#include "gapp_rec_fs8.h"
#include "cosmo_bk.h"
#include <time.h>

/* =============== main func ==================== */
int test_gapp_main(void)
{
    clock_t start, finish;
    int status;
    int i;
    double x, omm, h0, sig8, gam, ommh2;   
    // int N = 31;
    // char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC.txt";
    // double sigf=133.80863893, lenf=1.93610328;
    // int N = 23;
    // char * filename = "/home/ekli/myworks/cosmodata/OHD_23_BAO.txt";
    // double sigf=220.30220737, lenf=3.6377493;    
    int N = 51;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_51.txt";
    double sigf=170.75260325, lenf=2.61947963;

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();

    /* some test */
    for (i=0; i< 41; i++){
        x = 2.0*i/40.0;
        if(i == 0){
            printf("%5s %8s %9s ", "z", "H(z)", "sig");
            printf("%8s ", "lc_H");
            printf("%9s %9s %9s ", "dc", "sig", "lc_dc");
            printf("%9s %9s ", "Ez", "sig");
            printf("%9s %9s %9s ", "HD", "sig", "lc_HD");
            printf("%9s %9s %9s ", "fs8", "sig", "lc_fs8");
            printf("%9s %9s ", "fs8", "sig");
            printf("\n");
        }

        // rec mu
        printf("%5.2f %8.3f %9.6f ", x, rec_mu(x), sqrt( rec_covariance(x,x) ) );
        
        h0 = rec_mu(0.0);
        omm = 0.3;
        // initial lcdm
        initial_lcdm(omm, h0);
        // lcdm hz
        printf("%8.3f ", lcdm_hz(x));
        
        printf("%9.6f ", integrate_one_over_mu(0.0, x) );
        printf("%9.6f ", sqrt( integrate_cov_one_over_mu(0.0,x,0.0, x) ) );

        // lcdm dc
        printf("%9.6f ", lcdm_comoving_distance(x) );

        // rec mu/mu0
        printf("%9.6f %9.6f ", rec_mu_over_mu0(x),
                sqrt( rec_cov_mu_over_mu0(x,x) ) );

        // rec mu int 1/mu dx
        printf("%9.6f ", rec_mu_int_one_over_mu(x) );
        printf("%9.6f ", sqrt( rec_cov_mu_int_one_over_mu(x,x) ) );
        // lcdm H dc
        printf("%9.6f ", lcdm_comoving_distance(x)*lcdm_hz(x) );

        sig8 = 0.8;
        gam = 6.0/11.0;
        ommh2 = 0.3*pow(h0/100,2.0);
        // rec fs8
        printf("%9.6f ", rec_fs8_with_Hz(x, sig8, ommh2, gam));
        printf("%9.6f ", sqrt( rec_cov_fs8_with_Hz(x,x,sig8,ommh2,gam) ) );

        // lcdm fs8
        printf("%9.6f ", lcdm_growth(x, sig8, gam));

        // rec fs8 with Ez
        printf("%9.6f ", rec_fs8_with_Ez(x, sig8, omm, gam));
        printf("%9.6f ", rec_cov_fs8_with_Ez(x, x, sig8, omm, gam));

        printf("%9.6f ", rec_fs8_Mij_Hz(x,x,gam));

        printf("\n");
    }

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    
}

int test_gapp_fs8_loglike(void)
{
    int status;
    clock_t start, finish;
    // int i;
    // double x, omm, h0, sig8, gam, ommh2;   
    // int N = 31;
    // char * filename = "/home/ekli/myworks/cosmodata/OHD_31_CC.txt";
    // double sigf=133.80863893, lenf=1.93610328;
    // int N = 23;
    // char * filename = "/home/ekli/myworks/cosmodata/OHD_23_BAO.txt";
    // double sigf=220.30220737, lenf=3.6377493;    
    // int N = 51;
    // char * filename = "/home/ekli/myworks/cosmodata/OHD_51.txt";
    // double sigf=170.75260325, lenf=2.61947963;
    int N=54;
    char *filename = "/home/ekli/myworks/cosmodata/OHD_54_CC_BAO_P18.txt";
    double sigf=157.85759895, lenf=2.35798007;
    char * file_fs8 = "/home/ekli/myworks/fsigma8/code/src/fs8_ap_OHD_54_P18.txt";
    char * file_cov = "/home/ekli/myworks/fsigma8/code/src/fs8cov_ap_OHD_54_P18.txt";

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of gapp !\n");

    if((status = initial_fs8_obs(63, file_fs8, file_cov)) == _SUCCESS_ )
        printf("Successful initial fs8 data\n");
    else
        printf("Can not initial fs8 data\n");

    // loglike Hz
    //start = clock();
    //printf("%20.8e ", fs8_loglikelihood_Hz(0.8,0.147,6.0/11.0));
    //finish = clock();
    //printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    //// loglike Hz_t
    //start = clock();
    //printf("%20.8e ", fs8_loglikelihood_Hz_t(0.8,0.147,6.0/11.0));
    //finish = clock();
    //printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);
    //
    // loglike Hz_not
    start = clock();
    printf("%20.8e ", fs8_loglikelihood_Hz_not(0.8,0.147,6.0/11.0));
    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);
    
    // loglike Ez_not
    start = clock();
    printf("%20.8e ", fs8_loglikelihood_Ez_not(0.8,0.3,6.0/11.0));
    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    //printf("%20.8e ", fs8_loglikelihood_Ez(0.8,0.3,1.0/11.0));

    printf("\n");

    // free
    free_fs8obs();

    return _SUCCESS_;   
}

int test_gapp_loglike(void)
{
    int i;
    int N = 51;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_51.txt";
    //double sigf=170.75260325, lenf=2.61947963;
    double sigf, lenf;

    srand((unsigned)time(NULL));

    for (i=0; i<100; i++) {
        sigf = (double)rand()/RAND_MAX * 200.0;
        lenf = (double)rand()/RAND_MAX * 3.0;
        initial_gapp(N, filename);
        setup_gapp(sigf, lenf);
        printf("%9.3f %9.3f %20.9f\n", sigf, lenf, loglikelihood());
        free_gapp();
    }

    return _SUCCESS_;    

}

int test_gapp_panthon(void)
{
    clock_t start, finish;
    int status;
    int i;
    double x, omm, h0, sig8, gam, ommh2;   
    int N = 1048;
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    double sigf=24.46463208, lenf=0.14740619;

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();

    /* some test */
    for (i=0; i< 41; i++){
        x = 2.0*i/40.0;
        if(i == 0){
            printf("%5s %8s %9s ", "z", "H(z)", "sig");
            printf("%8s ", "lc_H");
            printf("%9s %9s %9s ", "dc", "sig", "lc_dc");
            printf("%9s %9s ", "Ez", "sig");
            printf("%9s %9s %9s ", "HD", "sig", "lc_HD");
            printf("%9s %9s %9s ", "fs8", "sig", "lc_fs8");
            printf("%9s %9s ", "fs8", "sig");
            printf("\n");
        }

        // rec mu
        printf("%5.2f %8.3f %9.6f ", x, rec_mu(x), sqrt( rec_covariance(x,x) ) );
        
        h0 = rec_mu(0.0);
        omm = 0.3;
        // initial lcdm
        initial_lcdm(omm, h0);
        // lcdm mu SN
        printf("%8.3f ", lcdm_mu_SN(x));
        //
        //printf("%9.6f ", integrate_one_over_mu(0.0, x) );
        //printf("%9.6f ", sqrt( integrate_cov_one_over_mu(0.0,x,0.0, x) ) );

        //// lcdm dc
        //printf("%9.6f ", lcdm_comoving_distance(x) );

        //// rec mu/mu0
        //printf("%9.6f %9.6f ", rec_mu_over_mu0(x),
        //        sqrt( rec_cov_mu_over_mu0(x,x) ) );

        //// rec mu int 1/mu dx
        //printf("%9.6f ", rec_mu_int_one_over_mu(x) );
        //printf("%9.6f ", sqrt( rec_cov_mu_int_one_over_mu(x,x) ) );
        //// lcdm H dc
        //printf("%9.6f ", lcdm_comoving_distance(x)*lcdm_hz(x) );

        //sig8 = 0.8;
        //gam = 6.0/11.0;
        //ommh2 = 0.3*pow(h0/100,2.0);
        //// rec fs8
        //printf("%9.6f ", rec_fs8_with_Hz(x, sig8, ommh2, gam));
        //printf("%9.6f ", sqrt( rec_cov_fs8_with_Hz(x,x,sig8,ommh2,gam) ) );

        //// lcdm fs8
        //printf("%9.6f ", lcdm_growth(x, sig8, gam));

        //// rec fs8 with Ez
        //printf("%9.6f ", rec_fs8_with_Ez(x, sig8, omm, gam));
        //printf("%9.6f ", rec_cov_fs8_with_Ez(x, x, sig8, omm, gam));

        //printf("%9.6f ", rec_fs8_Mij_Hz(x,x,gam));

        printf("\n");
    }

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    
}

int main(void)
{
    //test_gapp_main();
    //test_gapp_fs8_loglike();
    //test_gapp_loglike();
    test_gapp_panthon();
    return _SUCCESS_;
}
