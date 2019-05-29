#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// ciao
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include "functions.h"
#include "parameters.h"

int main(){ // ---  MAIN

    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
    // DECLAR VAR
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
    struct posforce PosForce_Cell[N_CELL], PosForce_Barcode[N_BARCODE];

    FILE  *xt, *xt_b, *out ;
    char  file_namext[350], file_namext_b[350], file_name[350];

    int a, i, k, l;
    int counter;
    int active_cells[N_CELL], active_barcodes[N_BARCODE], count_attachments[N_BARCODE];
    int count_cell, count_barcode, count_mix;
    int count_one, count_double, count_triple, count_quadruple;

    double deltax, deltay;
    double sedim_cell, sedim_barcode;



    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // SETTINGS
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

    printf("dt %lf \n", dt);

    printf("STOKES %lf\t OSEEN  %lf\n", STOKES_CELL, OSEEN );

    printf("Vol Fract: %lf\n", M_PI*(N_CELL*(R_CELL*R_CELL)+N_BARCODE*(R_BARCODE*R_BARCODE))/(L_BOX*L_BOX)  );

    /*  random  gen  alloc */
    gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);

    /* sid setting */
    gsl_rng_set(r,SEED);

    for(a=0;a<20; a++){

        printf("run %d \n", a);

    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // DATA FILES
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

        if(PRINT_CONFIG){

            sprintf(file_namext,"../numerical_data/config_cells_%dcells_%dbarcodes_noise-var%0.3lf.dat",N_CELL, N_BARCODE, LAMBDA_CELL);
            xt=fopen(file_namext, "w");

            sprintf(file_namext_b,"../numerical_data/config_barcodes_%dcells_%dbarcodes_noise-var%0.3lf.dat",N_CELL, N_BARCODE, LAMBDA_CELL);
            xt_b=fopen(file_namext_b, "w");
        }

        sprintf(file_name,"../numerical_data/events_%dcells_%dbarcodes_noise-var%0.3lf_run%d_sin.dat",N_CELL, N_BARCODE, LAMBDA_CELL, a);
        out = fopen(file_name, "w");



    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
    // INIZIALIZATION
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

        for (k=0; k< N_CELL; k++) {

                PosForce_Cell[k].x = 0.0;
                PosForce_Cell[k].y = 0.0 ;

            }

        for (k=0; k< N_BARCODE; k++) {

            PosForce_Barcode[k].x = 0.0 ;
            PosForce_Barcode[k].y = 30.0 ;

            }

    randomize_pos( &PosForce_Cell[0], 1+a, N_CELL, L_BOX);
    randomize_pos( &PosForce_Barcode[0], 3+a, N_BARCODE, L_BOX);

    counter = 0;
    count_barcode=0;
    count_cell = 0;
    count_mix =0;

    count_one = 0 ;
    count_double = 0;
    count_triple = 0;
    count_quadruple = 0;

    for(l=0; l< N_CELL; l++) active_cells[l] = 1;
    for(l=0; l< N_BARCODE; l++){

         active_barcodes[l] = 1;
         count_attachments[l] = 0;

        }

    for (i=0;i<=N_STEP; i++){ //---------------------------------------------- being temporal evolution

        sedim_cell = -V_CELL*sin(freq*(double)i*dt);
        sedim_barcode = -V_BARCODE*sin(freq*(double)i*dt);

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CELL DYN
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        for(l=0; l< N_CELL; l++){

            if( active_cells[l] == 1 ){

                deltax = compute_pos( (STOKES_CELL)*gsl_ran_gaussian(r,LAMBDA_CELL), 0.0 , dt)*dt ;
                deltay = compute_pos( (STOKES_CELL)*gsl_ran_gaussian(r,LAMBDA_CELL), sedim_cell, dt  )*dt ;

                PosForce_Cell[l].x += deltax ;
                PosForce_Cell[l].y += deltay ;

            }
        }

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
        // BARCODE DYN
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
        for(l=0; l< N_BARCODE; l++){

            if( active_barcodes[l] == 1 ){

                deltax = compute_pos( (STOKES_BARCODE)*gsl_ran_gaussian(r,LAMBDA_BARCODE), 0.0, dt )*dt ;
                deltay = compute_pos( (STOKES_BARCODE)*gsl_ran_gaussian(r,LAMBDA_BARCODE), sedim_barcode, dt  )*dt ;

                PosForce_Barcode[l].x += deltax ;
                PosForce_Barcode[l].y += deltay ;

            }
        }

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CONTROL CELL-BARCODE DIST
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        count_mix = control_cell_barcode_dist(count_mix, &active_cells[0], &active_barcodes[0], &count_attachments[0], &PosForce_Cell[0], &PosForce_Barcode[0], N_CELL, N_BARCODE, R_CELL, R_BARCODE);

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // CONTROL MULTIPLE ATTACHMENTS
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        count_one = 0 ;
        count_double = 0;
        count_triple = 0;
        count_quadruple = 0 ;

        for(l=0; l< N_BARCODE; l++){

            if(count_attachments[l] ==1 )count_one++;
            if(count_attachments[l] ==2 )count_double++;
            if(count_attachments[l] ==3 )count_triple++;
            if(count_attachments[l] ==4 )count_quadruple++;

        }

        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
        // PRINT ON FILES
        /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

        //if(i>N_STEP -10*40000)
        if(i%strobo==0){

            fprintf(out,"%lf\t %lf\t %lf\t %lf\t %lf\n", (double)i*dt, (100.*count_one)/((double)N_CELL), (100.*count_double)/((double)N_CELL), (100.*count_triple)/((double)N_CELL),  (100.*count_quadruple)/((double)N_CELL));

        if(PRINT_CONFIG){
            // cell config file print
            fprintf(xt, "%lf\t ", 0.001*dt*((double)i) );
            for(l=0;l<N_CELL; l++) fprintf(xt, "%lf\t %lf\t", PosForce_Cell[l].x, PosForce_Cell[l].y );
            fprintf(xt, "\n");

            // barcode config file print
            fprintf(xt_b, "%lf\t ", 0.001*dt*((double)i) );
            for(l=0;l<N_BARCODE; l++) fprintf(xt_b, "%lf\t %lf\t", PosForce_Barcode[l].x, PosForce_Barcode[l].y );
            fprintf(xt_b, "\n");

            }
        }

    } // ---------------------------------------------------------------- end temporal evolution

    fclose(out);

    }// -- end realizations

    gsl_rng_free(r);

return(0);

} // --- END MAIN
