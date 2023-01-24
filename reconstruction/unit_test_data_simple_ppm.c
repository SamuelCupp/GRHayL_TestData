#include "unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

static double eos_gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in);

int main(int argc, char **argv) {
  const double poison = 1e200;

  const double W_max = 10.0;
  const int neos = 1;
  const double rho_ppoly[1] = {0.0};
  const double Gamma_ppoly[1] = {2.0};
  const double k_ppoly0 = 1.0;
  const double Gamma_th = 2.0;

  eos_parameters eos;
  initialize_hybrid_eos_functions_and_params(W_max,
                                             poison, poison, poison,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &eos);

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

  double *rho = (double*) malloc(sizeof(double)*arraylength);
  double *press = (double*) malloc(sizeof(double)*arraylength);
  double *vx = (double*) malloc(sizeof(double)*arraylength);
  double *vz = (double*) malloc(sizeof(double)*arraylength);

  double *rhor   = (double*) malloc(sizeof(double)*arraylength);
  double *pressr = (double*) malloc(sizeof(double)*arraylength);
  double *vxr    = (double*) malloc(sizeof(double)*arraylength);
  double *vzr    = (double*) malloc(sizeof(double)*arraylength);

  double *rhol   = (double*) malloc(sizeof(double)*arraylength);
  double *pressl = (double*) malloc(sizeof(double)*arraylength);
  double *vxl    = (double*) malloc(sizeof(double)*arraylength);
  double *vzl    = (double*) malloc(sizeof(double)*arraylength);

  // There are two different reconstruction functions. One where
  // rho and P are reconstructed and one where they are not.
  double *vxr2 = (double*) malloc(sizeof(double)*arraylength);
  double *vzr2 = (double*) malloc(sizeof(double)*arraylength);
  double *vxl2 = (double*) malloc(sizeof(double)*arraylength);
  double *vzl2 = (double*) malloc(sizeof(double)*arraylength);

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  // Discontinuities or v>1 might give physically strange
  // output, but as long as cmin+cmax != 0 the function will
  // produce numbers.
  for(int k=1; k<dirlength; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        rho[index] = randf(0.0,100.0);

        double P_cold;
        eos.hybrid_compute_P_cold(&eos, rho[index], &P_cold);
        press[index] = randf(0.5*P_cold, 100*P_cold);

        vx[index] = randf(-1.0,1.0);
        vz[index] = randf(-1.0,1.0);
  }

  FILE* outfile;
  outfile = fopen("simple_ppm_initial_data.bin", "wb");
  check_file_was_successfully_open(outfile, "simple_ppm_initial_data.bin");
  fwrite(rho  , sizeof(double), arraylength, outfile);
  fwrite(press, sizeof(double), arraylength, outfile);
  fwrite(vx   , sizeof(double), arraylength, outfile);
  fwrite(vz   , sizeof(double), arraylength, outfile);
  fclose(outfile);

  const int num_vars = 2;

  for(int perturb=0; perturb<2; perturb++) {
#pragma omp parallel for
    // we are setting the reconstruction direction to x
    for(int k=1; k<dirlength; k++)
      for(int j=1; j<dirlength; j++)
        for(int i=3; i<dirlength-2; i++) {
          const int index = indexf(dirlength,i,j,k);
  
          double rho_stencil[6], press_stencil[6], v_flux_dir[6];
          double var_data[num_vars][6], var_datar[num_vars], var_datal[num_vars];
  
          for(int ind=0; ind<6; ind++) {
            const int stencil  = indexf(dirlength, i+ind-3, j, k); // PPM needs indices from -3 to +2
            v_flux_dir[ind]    = vx[stencil]; // Could be smaller; doesn't use full stencil
            rho_stencil[ind]   = rho[stencil];
            press_stencil[ind] = press[stencil];
            var_data[0][ind]   = vx[stencil];
            var_data[1][ind]   = vz[stencil];
          }
          const double gamma_eff = eos_gamma_eff(&eos, rho[index], press[index]);
  
          simple_ppm(rho_stencil, press_stencil, var_data,
                     num_vars, v_flux_dir, gamma_eff,
                     &rhor[index], &rhol[index],
                     &pressr[index], &pressl[index],
                     var_datar, var_datal);

          vxr[index] = var_datar[0];
          vzr[index] = var_datar[1];

          vxl[index] = var_datal[0];
          vzl[index] = var_datal[1];
    }

#pragma omp parallel for
    // we are setting the reconstruction direction to z
    for(int k=3; k<dirlength-2; k++)
      for(int j=1; j<dirlength; j++)
        for(int i=1; i<dirlength; i++) {
          const int index = indexf(dirlength,i,j,k);
  
          double press_stencil[6], v_flux_dir[6];
          double var_data[num_vars][6], var_datar[num_vars], var_datal[num_vars];
  
          for(int ind=0; ind<6; ind++) {
            const int stencil  = indexf(dirlength, i+ind-3, j, k); // PPM needs indices from -3 to +2
            v_flux_dir[ind]    = vz[stencil]; // Could be smaller; doesn't use full stencil
            press_stencil[ind] = press[stencil];
            var_data[0][ind]   = vx[stencil];
            var_data[1][ind]   = vz[stencil];
          }
          const double gamma_eff = eos_gamma_eff(&eos, rho[index], press[index]);
          simple_ppm_no_rho_P(press_stencil, var_data,
                     num_vars, v_flux_dir, gamma_eff,
                     var_datar, var_datal);

          vxr2[index] = var_datar[0];
          vzr2[index] = var_datar[1];

          vxl2[index] = var_datal[0];
          vzl2[index] = var_datal[1];
    }

    char filename[100];
    char suffix[10] = "";
    if(perturb) sprintf(suffix, "_pert");
    sprintf(filename,"simple_ppm%.5s.bin", suffix);
    outfile = fopen(filename,"wb");
    check_file_was_successfully_open(outfile, filename);
  
    fwrite(rhor  , sizeof(double), arraylength, outfile);
    fwrite(pressr, sizeof(double), arraylength, outfile);
    fwrite(vxr   , sizeof(double), arraylength, outfile);
    fwrite(vzr   , sizeof(double), arraylength, outfile);

    fwrite(rhol  , sizeof(double), arraylength, outfile);
    fwrite(pressl, sizeof(double), arraylength, outfile);
    fwrite(vxl   , sizeof(double), arraylength, outfile);
    fwrite(vzl   , sizeof(double), arraylength, outfile);

    fwrite(vxr2  , sizeof(double), arraylength, outfile);
    fwrite(vzr2  , sizeof(double), arraylength, outfile);
    fwrite(vxl2  , sizeof(double), arraylength, outfile);
    fwrite(vzl2  , sizeof(double), arraylength, outfile);
    fclose(outfile);
  }
}

static double eos_gamma_eff(const eos_parameters *restrict eos, const double rho_in, const double press_in) {
  double K, gamma;
  eos->hybrid_get_K_and_Gamma(eos, rho_in, &K, &gamma);
  const double P_cold = K*pow(rho_in, gamma);
  return eos->Gamma_th + (gamma - eos->Gamma_th)*P_cold/press_in;
}
