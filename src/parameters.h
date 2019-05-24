#define N_STEP  40000*3
#define N_PASSI_CAMPIONAMENTO_PER_CICLO 40000

#define N_CELL 100
#define N_BARCODE N_CELL
#define VEC_COMP 2

#define R_CELL 10.
#define R_BARCODE 20.
#define L_BOX 3000

#define VISCH2O 0.001 //#define VISC (6.23)*VISCH2O
#define VISC VISCH2O

#define TEMP 300
#define kBT 0.0000138*TEMP


#define OSEEN  1/(4*M_PI*VISC)
#define STOKES_CELL  1/(6*M_PI*VISC*R_CELL)
#define STOKES_BARCODE  1/(6*M_PI*VISC*R_BARCODE)

#define LAMBDA_CELL kBT/STOKES_CELL //fluctuation strength
#define LAMBDA_BARCODE kBT/STOKES_BARCODE //fluctuation strength

//0.004*0.19 //fluctuation strength

#define DENSTIY_WATER 1.
#define DENSTIY_CELL 1.05
#define DENSTIY_BARCODE 1.18

#define V_CELL   30.
#define V_BARCODE   3*4*V_CELL

#define dt  20./((double)N_PASSI_CAMPIONAMENTO_PER_CICLO)
#define strobo  20*5*1754
#define SEED 999999999
#define freq 0.2

#define PRINT_CONFIG 0
