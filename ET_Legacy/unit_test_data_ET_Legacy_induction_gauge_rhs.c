#include "unit_tests.h"

int main(int argc, char **argv) {


  const int gridmin     = 0;
  const int gridmax     = 21;
  const int arraylength = gridmax*gridmax*gridmax;

  const int gauss_center = 10;
  const double dX[3] = {0.1, 0.1, 0.1};

  const double Lorenz_damping_factor = 0.1;

  double *gupxx = (double*) malloc(sizeof(double)*arraylength);
  double *gupxy = (double*) malloc(sizeof(double)*arraylength);
  double *gupxz = (double*) malloc(sizeof(double)*arraylength);
  double *gupyy = (double*) malloc(sizeof(double)*arraylength);
  double *gupyz = (double*) malloc(sizeof(double)*arraylength);
  double *gupzz = (double*) malloc(sizeof(double)*arraylength);

  double *psi = (double*) malloc(sizeof(double)*arraylength);
  double *lapse = (double*) malloc(sizeof(double)*arraylength);
  double *betax = (double*) malloc(sizeof(double)*arraylength);
  double *betay = (double*) malloc(sizeof(double)*arraylength);
  double *betaz = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde = (double*) malloc(sizeof(double)*arraylength);
  double *Ax = (double*) malloc(sizeof(double)*arraylength);
  double *Ay = (double*) malloc(sizeof(double)*arraylength);
  double *Az = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shiftx_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shifty_interp = (double*) malloc(sizeof(double)*arraylength);
  double *shiftz_interp = (double*) malloc(sizeof(double)*arraylength);

  double *alpha_Phi_minus_betaj_A_j_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ax_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Ay_interp = (double*) malloc(sizeof(double)*arraylength);
  double *alpha_sqrtg_Az_interp = (double*) malloc(sizeof(double)*arraylength);

  double *phitilde_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ax_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Ay_rhs = (double*) malloc(sizeof(double)*arraylength);
  double *Az_rhs = (double*) malloc(sizeof(double)*arraylength);

  // Data which should be written before it is used is poisoned to validate behavior.
  // RHSs for A are set to 0 because they are assumed to already contain the zero-gauge
  // contribution to the RHS before entering this function.
  const double poison = 1e200;
  // This cannot be in parallel because the randomized quantities need to happen in the
  // same order every time.
  for(int k=gridmin; k<gridmax; k++)
    for(int j=gridmin; j<gridmax; j++)
      for(int i=gridmin; i<gridmax; i++) {
        const int index = indexf(gridmax,i,j,k);

        metric_quantities metric;
        double agxx, agxy, agxz, agyy, agyz, agzz;
        randomize_metric(&lapse[index], &agxx, &agxy, &agxz, &agyy, &agyz,
                         &agzz, &betax[index], &betay[index], &betaz[index]);

        initialize_metric(lapse[index], agxx, agxy, agxz, agyy, agyz, agzz,
                          betax[index], betay[index], betaz[index], &metric);

        psi[index]   = sqrt(metric.psi2);
        gupxx[index] = metric.psi4*metric.adm_gupxx;
        gupxy[index] = metric.psi4*metric.adm_gupxy;
        gupxz[index] = metric.psi4*metric.adm_gupxz;
        gupyy[index] = metric.psi4*metric.adm_gupyy;
        gupyz[index] = metric.psi4*metric.adm_gupyz;
        gupzz[index] = metric.psi4*metric.adm_gupzz;

        const int x = abs(i-gauss_center)*dX[0];
        const int y = abs(j-gauss_center)*dX[1];
        const int z = abs(k-gauss_center)*dX[2];
        const double r2 = x*x + y*y + z*z;

        phitilde[index] = exp(-r2/(2.0*1.0));;
        Ax[index]       = exp(-r2/(2.0*4.0));
        Ay[index]       = exp(-r2/(2.0*5.0));;
        Az[index]       = exp(-r2/(2.0*2.0));;

        alpha_interp[index]  = poison;
        shiftx_interp[index] = poison;
        shifty_interp[index] = poison;
        shiftz_interp[index] = poison;

        alpha_Phi_minus_betaj_A_j_interp[index] = poison;
        alpha_sqrtg_Ax_interp[index]            = poison;
        alpha_sqrtg_Ay_interp[index]            = poison;
        alpha_sqrtg_Az_interp[index]            = poison;

        Ax_rhs[index]       = 0.0;
        Ay_rhs[index]       = 0.0;
        Az_rhs[index]       = 0.0;
        phitilde_rhs[index] = poison;
  }

  FILE* outfile = fopen("ET_Legacy_induction_gauge_rhs_input.bin", "wb");
  check_file_was_successfully_open(outfile, "unit_test_ET_Legacy_induction_gauge_rhs_input.bin");

  fwrite(gupxx, sizeof(double), arraylength, outfile);
  fwrite(gupxy, sizeof(double), arraylength, outfile);
  fwrite(gupxz, sizeof(double), arraylength, outfile);
  fwrite(gupyy, sizeof(double), arraylength, outfile);
  fwrite(gupyz, sizeof(double), arraylength, outfile);
  fwrite(gupzz, sizeof(double), arraylength, outfile);

  fwrite(psi,   sizeof(double), arraylength, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);

  fwrite(phitilde, sizeof(double), arraylength, outfile);
  fwrite(Ax,       sizeof(double), arraylength, outfile);
  fwrite(Ay,       sizeof(double), arraylength, outfile);
  fwrite(Az,       sizeof(double), arraylength, outfile);

  fclose(outfile);

  for(int index=0; index<arraylength; index++) {
    
    gupxx[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupxy[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupxz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupyy[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupyz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    gupzz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    lapse[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    psi[index]      *= (1.0 + randf(-1,1)*1.0e-14);
    betax[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    betay[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    betaz[index]    *= (1.0 + randf(-1,1)*1.0e-14);
    phitilde[index] *= (1.0 + randf(-1,1)*1.0e-14);
    Ax[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Ay[index]       *= (1.0 + randf(-1,1)*1.0e-14);
    Az[index]       *= (1.0 + randf(-1,1)*1.0e-14);
  }

  outfile = fopen("ET_Legacy_induction_gauge_rhs_input_pert.bin", "wb");
  check_file_was_successfully_open(outfile, "unit_test_ET_Legacy_induction_gauge_rhs_input_pert.bin");

  fwrite(gupxx, sizeof(double), arraylength, outfile);
  fwrite(gupxy, sizeof(double), arraylength, outfile);
  fwrite(gupxz, sizeof(double), arraylength, outfile);
  fwrite(gupyy, sizeof(double), arraylength, outfile);
  fwrite(gupyz, sizeof(double), arraylength, outfile);
  fwrite(gupzz, sizeof(double), arraylength, outfile);

  fwrite(psi,   sizeof(double), arraylength, outfile);
  fwrite(lapse, sizeof(double), arraylength, outfile);
  fwrite(betax, sizeof(double), arraylength, outfile);
  fwrite(betay, sizeof(double), arraylength, outfile);
  fwrite(betaz, sizeof(double), arraylength, outfile);

  fwrite(phitilde, sizeof(double), arraylength, outfile);
  fwrite(Ax,       sizeof(double), arraylength, outfile);
  fwrite(Ay,       sizeof(double), arraylength, outfile);
  fwrite(Az,       sizeof(double), arraylength, outfile);

  fclose(outfile);
}
