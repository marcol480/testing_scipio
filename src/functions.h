#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// NEW DATA STRUCT
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

struct posforce{
    
    double x;
    double y;
    double fx;
    double fy;
    
};


/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
// DECLAR FUNCTIONS
/* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */

double compute_pos(double, double, double);

double compute_dist(struct posforce , struct posforce );

void randomize_pos(struct posforce *, double , int, double );

int control_cell_barcode_dist(int , int *, int *, int *, struct posforce *, struct posforce *, int, int, double, double);




