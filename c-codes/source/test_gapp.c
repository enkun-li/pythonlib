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
#include "gapp_sn_lens.h"
#include "pantheon.h"
#include "cosmo_bk.h"
#include "stronglens.h"
#include "gapp_rec_SN.h"
#include "cosmo_gp.h"
#include <time.h>

/* =============== main func ==================== */
int test_gapp_pantheon(void)
{
    clock_t start, finish;
    int status;
    int i;
    double x;//, omm, h0, sig8, gam, ommh2;   
    int N = 1048;
    //char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    //double sigf=24.46463208, lenf=0.14740619;
    //
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();

    /* some test */
    for (i=0; i< 41; i++){
        x = 0.0001 + 2.0*i/40.0;
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
        x = log(x);
        printf("%5.2f %9.6f ", exp(x), rec_mu(x));
        printf("%9.6f ", sqrt( rec_covariance(x,x) ) );
        
        printf("\n");
    }

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    
}

int test_gapp_pantheon_like(void)
{
    clock_t start, finish;
    int i;
    int N = 1048;
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    double sigf, lenf;

    srand((unsigned)time(NULL));

    initial_gapp(N, filename);
    initial_gapp_cov(N, filecov);

    for (i=0; i<50; i++) {
        sigf = (double)rand()/RAND_MAX * 200.0;
        lenf = (double)rand()/RAND_MAX * 3.0;
        setup_gapp(sigf, lenf);
        printf("%9.3f %9.3f %20.9f\n", sigf, lenf, loglikelihood());
    }
    
    start = clock();

    /* some test */

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

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

    initial_gapp(N, filename);

    for (i=0; i<100; i++) {
        sigf = (double)rand()/RAND_MAX * 200.0;
        lenf = (double)rand()/RAND_MAX * 3.0;
        setup_gapp(sigf, lenf);
        printf("%9.3f %9.3f %20.9f\n", sigf, lenf, loglikelihood());
    }
    
    free_gapp();
    return _SUCCESS_;    

}

int test_gapp_sl_loglike(void)
{
    clock_t start, finish;
    int status;
    int i;
    double x, MB, omk, fE;//, omm, h0, sig8, gam, ommh2;   
    int N = 1048;
    //char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    //double sigf=24.46463208, lenf=0.14740619;
    //
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;

    char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs.txt";
    FILE *fp=NULL;

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();
    
    srand((unsigned)time(NULL));

    // initial strong lensing
    printf("Initial strong lensing data \n");
    //initial_SL(205, file_SL);
    printf("Done\n");

    fp = fopen(file_SL, "r");

    /* some test */
    for (i=0; i< 41; i++){
        x = 0.0001 + 2.0*i/40.0;
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
        x = log(x);
        printf("%5.2f %9.6f ", exp(x), rec_mu(x));
        printf("%9.6f ", sqrt( rec_covariance(x,x) ) );

        MB = ((double)rand()/RAND_MAX -0.5)*3 +23.83;
        omk = (double)rand()/RAND_MAX -0.5;
        fE = ((double)rand()/RAND_MAX -0.5)*0.1+1.0;
        printf("%12.6f %12.6f %12.6f ", MB, omk, fE);

        printf("%12.6f ", rec_distance_noM(exp(x))/exp(log(10)/5*23.83) );
        printf("%12.6f ", sqrt( rec_cov_distance_noM(exp(x), exp(x))
                    /pow(exp(log(10)/5*23.83), 2.0)) );
        //printf("%20.8e ", SL_loglikelihood(MB,omk,fE));

        printf("\n");
    }

    fclose(fp);

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    

}

int test_output_mean_cov_SL(void)
{
    clock_t start, finish;
    int status;
    //int i;
    //double x, MB, omk, fE;//, omm, h0, sig8, gam, ommh2;   
    int N = 1048;
    //char * filename = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;

    //int NN= 205;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs.txt";
    //int NN = 173;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real.txt";
    //int NN = 167;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_SN.txt";
    //int NN = 138;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real_SN.txt";
    //int NN=111;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real_1.4.txt";
    //int NN=120;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real_1.55.txt";
    //int NN=126;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real_1.8.txt";
    //int NN=129;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_real_2.0.txt";
    //char * out_mean = "SL_mean.txt";
    //char * out_cov = "SL_cov.txt";
    //int NN=54;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_lowmass.txt";
    //char * out_mean = "SL_mean_lowmass.txt";
    //char * out_cov = "SL_cov_lowmass.txt";
    //int NN=136;
    //char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_intmass.txt";
    //char * out_mean = "SL_mean_intmass.txt";
    //char * out_cov = "SL_cov_intmass.txt";
    int NN=15;
    char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs_higmass.txt";
    char * out_mean = "SL_mean_higmass.txt";
    char * out_cov = "SL_cov_higmass.txt";

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();
    
    srand((unsigned)time(NULL));

    // initial strong lensing
    printf("Initial strong lensing data \n");
    output_mean_cov_of_SL(NN, file_SL, out_mean, out_cov);
    //initial_SL(205, file_SL);
    printf("Done\n");

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    

}

int test_output_mean_cov_TDSL(void)
{
    clock_t start, finish;
    int status;
    int N = 1048;
    char * filecov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    char * filename = "/home/ekli/myworks/cosmodata/sn_full_lnz_dat.txt";
    double sigf=39.21533362, lenf=12.97267722;
    char * file_SL = "/home/ekli/myworks/cosmodata/time_delay_obs.txt";

    if((status = initial_gapp(N, filename)) == _SUCCESS_
            && (status = initial_gapp_cov(N, filecov)) == _SUCCESS_
            && (status = setup_gapp(sigf, lenf)) == _SUCCESS_ 
            )
        printf("Successful initial gapp\n");
    else
        printf("Can not setup covariance of !\n");
    
    start = clock();
    
    srand((unsigned)time(NULL));

    // initial strong lensing
    printf("Initial strong lensing data \n");
    output_mean_cov_of_TDSL(2, file_SL);
    //initial_SL(205, file_SL);
    printf("Done\n");

    finish = clock();
    printf("Cost time is %f s\n", (double)(finish - start)/CLOCKS_PER_SEC);

    free_gapp();

    return _SUCCESS_;    

}

int test_pantheon_like(void)
{
    int N = 1048;
    char * sn_dat = "/home/ekli/myworks/cosmodata/sn_full_dat.txt";
    char * sn_cov = "/home/ekli/myworks/cosmodata/sn_full_cov.txt";
    int i;
    double omm, omk, h0, chisq;

    initial_Pantheon(N, sn_dat, sn_cov);

    srand((unsigned)time(NULL));

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        omm = (double)rand()/RAND_MAX * 0.5;
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        printf("%12.6f %9.6f %9.6f ", h0, omm, omk);
        initial_lcdm(omm, omk, h0);
        chisq = sn_loglikelihood(&lcdm_mu_SN);
        printf("%20.8e ", chisq);
        printf("\n");
    }

    free_Pantheon();
    
    return _SUCCESS_;
}

int test_lcdm_stronglens_like(void)
{
    int i;
    double omm, omk, h0, fE, chisq;
    
    srand((unsigned)time(NULL));

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        omm = (double)rand()/RAND_MAX * 0.5;
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        fE = (double)rand()/RAND_MAX * 0.6 + 0.8;

        printf("%12.6f %9.6f %9.6f %9.6f", h0, omm, omk, fE);
        chisq = lcdm_stronglens(h0,omm,omk,fE);
        printf("%20.8e ", chisq);
        printf("\n");
    }
    
    free_stronglens();
    return _SUCCESS_;
}

int test_cosmography_stronglens_like(void)
{
    int i;
    double omk, a1, a2, fE, chisq;
    
    srand((unsigned)time(NULL));

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        a1 = (double)rand()/RAND_MAX * 2 - 1;
        a2 = (double)rand()/RAND_MAX * 5 -2.5;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        fE = (double)rand()/RAND_MAX * 0.6 + 0.8;

        printf("%12.6f %9.6f %9.6f %9.6f", a1, a2, omk, fE);
        chisq = cosmography_stronglens(omk,a1,a2,fE);
        printf("%20.8e ", chisq);
        printf("\n");
    }

    free_stronglens();
    return _SUCCESS_;
}

int test_return_SLlike(void)
{
    int i;
    double MB, omk, h0, fE, chisq;
    
    srand((unsigned)time(NULL));
    printf("test return SL loglike\n");

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        MB = (double)rand()/RAND_MAX * 3 - 21.0;
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        fE = (double)rand()/RAND_MAX * 0.6 + 0.8;

        printf("%12.6f %9.6f %9.6f %9.6f", h0, MB, omk, fE);
        chisq = return_SL_loglike(MB,omk,fE);
        printf("%20.8e ", chisq);
        printf("\n");
    }
    
    free_SL();
    return _SUCCESS_;
}

int test_margin_of_MB_SLlike(void)
{
    int i;
    double omk, h0, fE, chisq;
    
    srand((unsigned)time(NULL));
    printf("test return mag MB SL loglike\n");

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        fE = (double)rand()/RAND_MAX * 0.6 + 0.8;

        printf("%12.6f %9.6f %9.6f", h0, omk, fE);
        chisq = margin_of_MB_SL_loglike(omk,fE);
        printf("%20.8e ", chisq);
        if(chisq > 0)
            printf("************");
        if(chisq < -1e29)
            printf("+++++++++++++");
        printf("\n");
    }
    
    free_SL();
    return _SUCCESS_;
}

int test_return_TDSLlike(void)
{
    int i;
    double MB, omk, h0, fE, chisq;
    
    srand((unsigned)time(NULL));

    printf("this is test return tdsl loglike \n");

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        MB = (double)rand()/RAND_MAX * 3 - 21.0;
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;
        fE = (double)rand()/RAND_MAX * 0.6 + 0.8;

        printf("%12.6f %9.6f %9.6f %9.6f", h0, MB, omk, fE);
        chisq = return_TDSL_loglike(MB,omk,h0);
        printf("%20.8e ", chisq);
        if(chisq > 0)
            printf("**************");
        if(chisq < -1e29)
            printf("+++++++++++++");
        printf("\n");
    }
    
    free_TDSL();
    return _SUCCESS_;

}

int test_margin_of_MB_TDSLlike(void)
{
    int i;
    double omk, h0, chisq;
    
    srand((unsigned)time(NULL));
    printf("test return mag MB TDSL loglike\n");

    for (i=0; i<50; i++) {
        printf("%4d ", i);
        h0 = (double)rand()/RAND_MAX * 20.0 + 60;
        omk = (double)rand()/RAND_MAX *2.0 -1.0;

        printf("%12.6f %9.6f ", h0, omk);
        chisq = margin_of_MB_TDSL_loglike(omk,h0);
        printf("%20.8e ", chisq);
        if(chisq > 0)
            printf("**************");
        if(chisq < -1e29)
            printf("+++++++++++++");

        printf("\n");
    }
    
    free_TDSL();
    return _SUCCESS_;
}

int test_return_Pantheon_dz(void)
{
    int N = 51;
    char * filename = "/home/ekli/myworks/cosmodata/OHD_51.txt";
    double sigf=170.75260325, lenf=2.61947963;

    initial_gapp(N, filename);
    setup_gapp(sigf, lenf);
    
    printf("Start out put pantheon dz\n");
    return_Pantheon_dz();
    
    free_gapp();
    return _SUCCESS_;
}

int test_return_gp_ohd_loglike(void)
{
    int i;
    double omk, chisq;

    setup_cosmo_gp(0.0, 0);
    //initial_cosmo_gp(0);
    //GP_OHD_Pantheon_output_dc(0);
    
    srand((unsigned)time(NULL));
    printf("\nOutput gp ohd loglike\n");

    for (i=0; i< 50; i++) {
        omk = (double)rand()/RAND_MAX *2-1.0;
        printf("omk = %9.6f", omk);

        chisq = return_gp_ohd_loglike(omk);
        printf("like = %20.8e", chisq);
        if(chisq < -2.9e20)
            printf("+++++++++++++");
        if(chisq > 0.0)
            printf("*************");
        printf("\n");
    }
    return _SUCCESS_;
}

int test_GP_OHD_pantheonout(void)
{
    int num=1;

    setup_cosmo_gp(0.0, num);
    initial_cosmo_gp(num);
    GP_OHD_Pantheon_output_dc(num);
    
    return _SUCCESS_;
}

int test_GP_OHD_stronglensout(void)
{
    int N = 205;
    char * file_SL = "/home/ekli/myworks/cosmodata/strong_lens_obs.txt";
    char * file_dc = "./data/SL_obs.txt";
    
    GP_OHD_StrongLens_output_dc(0, N, file_SL, file_dc);
    return _SUCCESS_;
}

int main(void)
{
    //test_gapp_main();
    //test_gapp_fs8_loglike();
    //test_gapp_loglike();
    //test_gapp_pantheon();
    //test_gapp_pantheon_like();
    //test_gapp_sl_loglike();
    //test_output_mean_cov_SL();
    //test_output_mean_cov_TDSL();
    //test_pantheon_like();
    //test_lcdm_stronglens_like();
    //test_cosmography_stronglens_like();
    //test_return_SLlike();
    //test_return_TDSLlike();
    //test_return_Pantheon_dz();
    //test_margin_of_MB_SLlike();
    //test_margin_of_MB_TDSLlike();
    //initial_SL_usedat(4);
    test_GP_OHD_pantheonout();
    //test_return_gp_ohd_loglike();
    //test_GP_OHD_stronglensout();
    return _SUCCESS_;
}
