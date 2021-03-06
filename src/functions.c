#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#include <check.h>


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// NEW DATA STRUCT
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

struct posforce{

    double x;
    double y;
    double fx;
    double fy;

};


/* - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - */
//
//                   --- FUNCTIONS --
//
/* - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ - */


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// Randomize positions
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */


void randomize_pos(struct posforce *X, double fraction, int N, double L_BOX){
    int k;

    gsl_rng * r = gsl_rng_alloc (gsl_rng_ranlxs0);
    gsl_rng_set(r,12345*((unsigned long int)fraction));


    for (k=0; k< N; k++) {

            X[k].x =  L_BOX*gsl_rng_uniform(r) ;
            X[k].y = L_BOX*gsl_rng_uniform(r) ;

        }
    }


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// Compute distance
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

double compute_dist(struct posforce X1, struct posforce X2){

    double dist ;

    dist  = (X1.x - X2.x)*(X1.x - X2.x)  ;
    dist += (X1.y - X2.y)*(X1.y - X2.y)  ;
    dist = sqrt(dist);

    return dist;
}



/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
// Integrate pos
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

double compute_pos(double v , double coeff, double dt ){

    double pos;

    pos  = v + coeff;

    if(pos*dt > 0.3) printf("position = %g troppo!\n", pos), exit(0);
    if(isnan(pos) == 1) printf("position cell NANNA!"), exit(0);

    return pos;
}



/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */
// Check separations
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

// -- CELL
int control_cell_barcode_dist(int count, int  *active_cells, int *active_barcodes, int *count_attachments, struct posforce *X, struct posforce *Y, int N_CELL, int N_BARCODE, double R_CELL, double R_BARCODE){

    double dist_mix[N_CELL][N_BARCODE];
    int l, k;

    for(l=0; l< N_CELL; l++){
        for(k=0; k< N_BARCODE; k++){

            if( active_cells[l]==1 && active_barcodes[k]==1 ){

                dist_mix[l][k] = compute_dist(X[l], Y[k]);

                if(dist_mix[l][k] <  0.5*(R_CELL+R_BARCODE) ){

                    active_cells[l] = 0;
                    //active_barcodes[k] = 0;
                    count_attachments[k] += 1;
                    if(count_attachments[k]> 1); // printf("multiple attachment, N= %d, barcode %d, cell %d \n", count_attachments[k], k, l );
                    else count ++;

                }
            }
        }
    }

    return count;
}
