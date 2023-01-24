#include "unit_tests.h"
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))

// each point requires phi from -1 to +2
// c{1,2}{min,max}, v{1,2}{r,l}{r,l} B{1,2}{r,l} from 0 to +1
void A_rhs_dir(const int perturb, const int dirlength, const int A_dir,
               const double *restrict phi_bssn,
               const double *restrict cmin_1, const double *restrict cmax_1,
               const double *restrict cmin_2, const double *restrict cmax_2,
               const double *restrict v1rr, const double *restrict v1rl,
               const double *restrict v1lr, const double *restrict v1ll,
               const double *restrict v2rr, const double *restrict v2rl,
               const double *restrict v2lr, const double *restrict v2ll,
               const double *restrict B1r, const double *restrict B1l,
               const double *restrict B2r, const double *restrict B2l,
               double *restrict A_rhs);

int main(int argc, char **argv) {

  const int dirlength = 20;
  const int arraylength = dirlength*dirlength*dirlength;

  double *phi_bssn = (double*) malloc(sizeof(double)*arraylength);

  double *cmin[3];
  double *cmax[3];
  for(int i=0; i<3; i++) {
    cmin[i] = (double*) malloc(sizeof(double)*arraylength);
    cmax[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *vrr[3];
  double *vrl[3];
  double *vlr[3];
  double *vll[3];
  for(int i=0; i<3; i++) {
    vrr[i] = (double*) malloc(sizeof(double)*arraylength);
    vrl[i] = (double*) malloc(sizeof(double)*arraylength);
    vlr[i] = (double*) malloc(sizeof(double)*arraylength);
    vll[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *Br[3];
  double *Bl[3];
  for(int i=0; i<3; i++) {
    Br[i] = (double*) malloc(sizeof(double)*arraylength);
    Bl[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  double *A_rhs[3];
  for(int i=0; i<3; i++) {
    A_rhs[i] = (double*) malloc(sizeof(double)*arraylength);
  }

  // Initialize random data. Note that for this test,
  // we needn't worry too much with physical reasonableness.
  // Discontinuities or v>1 might give physically strange
  // output, but as long as cmin+cmax != 0 the function will
  // produce numbers.
  for(int k=1; k<dirlength; k++)
    for(int j=1; j<dirlength; j++)
      for(int i=1; i<dirlength; i++) {
        const int index = indexf(dirlength,i,j,k);

        phi_bssn[index] = randf(-10.0,10.0);

        for(int coord=0; coord<3; coord++) {
          Br[coord][index]  = randf(-1.0,1.0);
          Bl[coord][index]  = randf(-1.0,1.0);
          vrr[coord][index] = randf(-1.0,1.0);
          vrl[coord][index] = randf(-1.0,1.0);
          vlr[coord][index] = randf(-1.0,1.0);
          vll[coord][index] = randf(-1.0,1.0);

          bool zero = true;
          while(zero) {
            cmin[coord][index] = randf(-10.0,10.0);
            cmax[coord][index] = randf(-10.0,10.0);
            zero = cmin[coord][index]+cmax[coord][index] == 0.0;
          }
          A_rhs[coord][index] = 0.0;
        }
  }

  FILE* outfile;
  outfile = fopen("A_no_gauge_rhs_initial_data.bin", "wb");
  check_file_was_successfully_open(outfile, "A_no_gauge_rhs_initial_data.bin");

  fwrite(phi_bssn, sizeof(double), arraylength, outfile);

  for(int coord=0; coord<3; coord++) {
    fwrite(Br[coord], sizeof(double), arraylength, outfile);
    fwrite(Bl[coord], sizeof(double), arraylength, outfile);

    fwrite(vrr[coord], sizeof(double), arraylength, outfile);
    fwrite(vrl[coord], sizeof(double), arraylength, outfile);
    fwrite(vlr[coord], sizeof(double), arraylength, outfile);
    fwrite(vll[coord], sizeof(double), arraylength, outfile);

    fwrite(cmin[coord], sizeof(double), arraylength, outfile);
    fwrite(cmax[coord], sizeof(double), arraylength, outfile);
  }
  fclose(outfile);

  for(int perturb=0; perturb<2; perturb++) {
    for(int A_dir=1; A_dir<4; A_dir++) {
    int dir1 = A_dir%3, dir2 = (A_dir+1)%3;
    A_rhs_dir(perturb, dirlength, A_dir, phi_bssn,
              cmin[dir1], cmax[dir1], cmin[dir2], cmax[dir2],
              vrr[dir1], vrl[dir1], vlr[dir1], vll[dir1],
              vrr[dir2], vrl[dir2], vlr[dir2], vll[dir2],
              Br[dir1], Bl[dir1], Br[dir2], Bl[dir2],
              A_rhs[A_dir-1]);
    }


    char filename[100];
    char suffix[10] = "";
    if(perturb) sprintf(suffix, "_pert");
    sprintf(filename,"A_no_gauge_rhs%.5s.bin", suffix);
    outfile = fopen(filename,"wb");
    check_file_was_successfully_open(outfile, filename);
  
    for(int coord=0; coord<3; coord++)
      fwrite(A_rhs[coord], sizeof(double), arraylength, outfile);
    fclose(outfile);
  }
}

void A_rhs_dir(const int perturb,
               const int dirlength,
               const int A_dir,
               const double *restrict phi_bssn,
               const double *restrict cmin_1,
               const double *restrict cmax_1,
               const double *restrict cmin_2,
               const double *restrict cmax_2,
               const double *restrict v1rr,
               const double *restrict v1rl,
               const double *restrict v1lr,
               const double *restrict v1ll,
               const double *restrict v2rr,
               const double *restrict v2rl,
               const double *restrict v2lr,
               const double *restrict v2ll,
               const double *restrict B1r,
               const double *restrict B1l,
               const double *restrict B2r,
               const double *restrict B2l,
               double *restrict A_rhs) {

  const int xdir = (A_dir==1);
  const int ydir = (A_dir==2);
  const int zdir = (A_dir==3);

  // This offsets the index by +1 in the perpendicular directions
  const int v_offset[3] = { !xdir, !ydir, !zdir };

  // This offsets the index by +1 in the permuted direction (x<-y<-z)
  const int B1_offset[3] = { ydir, zdir, xdir };

  // This offsets the index by +1 in the permuted direction (x->y->z)
  const int B2_offset[3] = { zdir, xdir, ydir };

  const int imax = dirlength-2*!xdir;
  const int jmax = dirlength-2*!ydir;
  const int kmax = dirlength-2*!zdir;

#pragma omp parallel for
  for(int k=1; k<kmax; k++)
    for(int j=1; j<jmax; j++)
      for(int i=1; i<imax; i++) {
        const int index    = indexf(dirlength,i,j,k);
        const int index_v  = indexf(dirlength,i+v_offset[0], j+v_offset[1], k+v_offset[2]);
        const int index_B1 = indexf(dirlength,i+B1_offset[0],j+B1_offset[1],k+B1_offset[2]);
        const int index_B2 = indexf(dirlength,i+B2_offset[0],j+B2_offset[1],k+B2_offset[2]);

        A_no_gauge_vars vars;

        // This computes psi6 at the point staggered with respect to the two perpendicular
        // directions using the variable phi, which naturally lives at (i, j, k).
        // E.g. A_x needs phi at (i, j+1/2, k+1/2), so it must be interpolated to that point.
        // With the IPH macro, we first interpolate to the points
        // (i, j+1/2, k-1), (i, j+1/2, k), (i, j+1/2, k+1), (i, j+1/2, k+2) and use
        // those to compute phi at (i, j+1/2, k+1/2).
        vars.psi6 =
          (1.0 + perturb*randf(-1,1)*1.0e-14) *
          exp(6.0*IPH(
            IPH(phi_bssn[indexf(dirlength,i-!xdir  , j-xdir  -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i        , j       -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir  , j+xdir  -zdir,   k-!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir-zdir,   k-!zdir)]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir,          k      )],
                phi_bssn[index],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir,          k      )],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir,        k      )]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir  +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i,         j       +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir  +zdir,   k+!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir+zdir,   k+!zdir)]),

            IPH(phi_bssn[indexf(dirlength,i-!xdir,   j-xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i,         j       +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i+!xdir,   j+xdir  +2*zdir, k+2*!zdir)],
                phi_bssn[indexf(dirlength,i+2*!xdir, j+2*xdir+2*zdir, k+2*!zdir)])));

        vars.v1rr = v1rr[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v1rl = v1rl[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v1lr = v1lr[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v1ll = v1ll[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);

        vars.v2rr = v2rr[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v2rl = v2rl[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v2lr = v2lr[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.v2ll = v2ll[index_v] * (1.0 + perturb*randf(-1,1)*1.0e-14);

        vars.B1r = B1r[index_B1] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.B1l = B1l[index_B1] * (1.0 + perturb*randf(-1,1)*1.0e-14);

        vars.B2r = B2r[index_B2] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.B2l = B2l[index_B2] * (1.0 + perturb*randf(-1,1)*1.0e-14);

        vars.c1_min = cmin_1[index_B2] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.c1_max = cmax_1[index_B2] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.c2_min = cmin_2[index_B1] * (1.0 + perturb*randf(-1,1)*1.0e-14);
        vars.c2_max = cmax_2[index_B1] * (1.0 + perturb*randf(-1,1)*1.0e-14);

        A_rhs[index] = A_no_gauge_rhs_stencil(&vars);
  }
}
