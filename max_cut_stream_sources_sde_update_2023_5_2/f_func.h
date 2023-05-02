
#include <ap_fixed.h>
#include <ap_int.h>
#include <hls_math.h>
#include <hls_stream.h>
#include <string.h>

#include "stdio.h"
#include "fstream"
#include "iostream"
//#include "hls_streamofblocks.h"

const int precision_level=15;

typedef ap_fixed<32,17> FIXED_POINT;//balanced range and precision,
typedef ap_fixed<32,19> FIXED_POINT_EXP;//EXTENDED RANGE FOR DIVISION in COrdic

typedef ap_fixed<32,15> FIXED_POINT_INJECTION;//SPECIAL OPTIMIZATION FOR INJECTION CALC,should be 3 bit higher than FIXED_POINT_PRECISION
typedef ap_fixed<32,13> FIXED_POINT_PRECISION;//used for phi_delta, phi_diff, high precision, low abs value, smaller than 7 bit will not work,,
typedef ap_fixed<32,6> FIXED_POINT_U_PRECISION;// very high precision, sin result. etc,no -1 to 1 value
//typedef ap_fixed<60,30> FIXED_POINT_RANGE;//high range low,precision, used for some large mult,rarely used
typedef ap_fixed<32,13> FIXED_POINT_LOW;//USED FOR k,K, and some constants small value and no to low precision requirement
typedef ap_fixed<20,1> FIXED_POINT_RANDOM;

static int K_mode;


const FIXED_POINT_PRECISION pi=3.141592653589793288;
const FIXED_POINT_U_PRECISION pi_u=3.141592653589793288;
const FIXED_POINT_U_PRECISION pi_inverse_u=0.31830988618379;

const FIXED_POINT_PRECISION pi_r=6.28318530717958;
const FIXED_POINT_PRECISION pi_r_inverse=0.15915494309189;
const FIXED_POINT_PRECISION pi_q1=1.57079632797489;
const FIXED_POINT_PRECISION pi_q2=3.141592653589793288;
const FIXED_POINT_PRECISION pi_q3=4.71238898038468;
const FIXED_POINT_PRECISION sigma_test= 0.1;


const int CORDIC_PI_ITERATION=10;
const int CORDIC_SIN_ITERATION=19;


//const FIXED_POINT_U_PRECISION cordic_lut_step=FIXED_POINT_U_PRECISION();


//parmeter tuning N_size ->F->FINAL_SUM->D2F->F2R->PARTIAO_SUM_HEIGHT_PARTITION

const int N_size=10240;


//N10240 J_valid 10000 F 128 ratio 16

const int J_VALID_SIZE=16000;
const int DATAFLOW_PARALLEL_F=128;//

const int BURST_ratio=16;
const int BURST_d_out=DATAFLOW_PARALLEL_F/BURST_ratio;
const int BURST_d_in=BURST_d_out;
const int BURST_j_in=BURST_d_out;
const int BURST_dwt_in=BURST_d_out;
const int ENTRY_PER_BURST=BURST_j_in/2;


const int J_LOOP_IN_INNER= DATAFLOW_PARALLEL_F*2/BURST_j_in;


const int PARTIAL_SUM_SIZE=(N_size*N_size)/DATAFLOW_PARALLEL_F;
const int PARTIAL_SUM_WIDH=N_size;
const int PARTIAL_SUM_HEIGHT=(N_size/DATAFLOW_PARALLEL_F);




const int DATAFLOW_TOP_DELTA_ITERATION=N_size/DATAFLOW_PARALLEL_F;
const int DATAFLOW_TOP_PRE_SUM_ITERATION=N_size*N_size/DATAFLOW_PARALLEL_F;



//performance metrics
//const int LATENCY_D_MAX=N_size/DATAFLOW_PARALLEL_F*3;
//const int LATENCY_D_MIN=N_size/DATAFLOW_PARALLEL_F*1;

//const int LATENCY_P_MAX=PARTIAL_SUM_HEIGHT/DATAFLOW_PARALLEL_P*3;
//const int LATENCY_P_MIN=PARTIAL_SUM_HEIGHT/DATAFLOW_PARALLEL_P*1;



typedef struct phi_f_block
{
    FIXED_POINT_U_PRECISION phi_delta[DATAFLOW_PARALLEL_F]={0};
}phi_f_block;



const FIXED_POINT_PRECISION stv_sq2=0.0001;//2*0.002*0.002

//customized cordic
void cordic_sin(FIXED_POINT_PRECISION theta, FIXED_POINT_U_PRECISION* sin);
void cordic_sin_default(FIXED_POINT_PRECISION theta, FIXED_POINT_U_PRECISION*sin);
void cordic_sin_piecewise(FIXED_POINT_PRECISION theta, FIXED_POINT_U_PRECISION* sin);



//function f 
void f_func_accumulate(FIXED_POINT_U_PRECISION theta,FIXED_POINT_PRECISION* result);
void f_func_eq(FIXED_POINT_LOW lower_k,FIXED_POINT_U_PRECISION theta, FIXED_POINT_PRECISION* result);// output range -10 +10
void f_func_gaussmf(FIXED_POINT_U_PRECISION x, FIXED_POINT_PRECISION mean, FIXED_POINT_PRECISION* result);
void f_func_lut(FIXED_POINT_U_PRECISION theta, FIXED_POINT_U_PRECISION* result,int k_part);




//dataflow stream version, reduce synthesize difficulty on array partitionos, reduce synchronization memory, so we can support large matrices
void dataflow_delta_phase(int total_block_num,int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_row[J_VALID_SIZE],hls::stream <int>&current_row_single,FIXED_POINT_U_PRECISION phi_diff_in_remap[DATAFLOW_PARALLEL_F+1][N_size], 
hls::stream<phi_f_block>&phi_delta_struct_F,hls::stream <FIXED_POINT_U_PRECISION>&phi_diff_single);

void dataflow_pre_sum_phase(hls::stream<int>&current_row_pre_accum_single, hls::stream<int>&current_row_partial_sum_single,hls::stream <phi_f_block>&phi_delta_struct_F,hls::stream <FIXED_POINT_PRECISION>&sin_term_single,int k_part,int total_block_num);
void dataflow_partial_sum_phase(hls::stream <FIXED_POINT_PRECISION>&sin_partial_single, hls::stream<int>&current_row_single, hls::stream <FIXED_POINT_PRECISION> &kuramoto_single,int total_block_num);
void dataflow_final_phase(hls::stream <FIXED_POINT_PRECISION>&kuramoto_sum_single,hls::stream <FIXED_POINT_U_PRECISION>&phi_diff_final,hls::stream <FIXED_POINT_U_PRECISION>&kuramoto_final_single,FIXED_POINT_PRECISION C_sync_fix, FIXED_POINT_LOW C1_fix,int k_part,int size);
void dataflow_sde_accumulate(hls::stream <FIXED_POINT_U_PRECISION>&kuramoto_final_single,FIXED_POINT_U_PRECISION phi_next[N_size],int size);



//PRE-SUM PHASE
//void dataflow_phi_delta_stream( block_phi_delta_F phi_delta_blocks_F,FIXED_POINT_PRECISION phi_delta_dataflow[DATAFLOW_PARALLEL_F]);
void dataflow_f_func_calc(FIXED_POINT_U_PRECISION phi_delta_in[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION f_result_out[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION phi_delta_sync_out[DATAFLOW_PARALLEL_F],int k_part);
void dataflow_kuramor_pre_accum_calc(FIXED_POINT_U_PRECISION phi_delta_in[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION f_result_in[DATAFLOW_PARALLEL_F], FIXED_POINT_U_PRECISION sin_term_out[DATAFLOW_PARALLEL_F]);
void dataflow_sin_term_out_calc(FIXED_POINT_U_PRECISION sin_term_in[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION sin_term_out[DATAFLOW_PARALLEL_F]);
void sin_calc(FIXED_POINT_U_PRECISION f_result,FIXED_POINT_U_PRECISION phi_delta,FIXED_POINT_U_PRECISION* result_out);
void dataflow_sin_accumulate_calc(FIXED_POINT_U_PRECISION sin_term_in[DATAFLOW_PARALLEL_F],FIXED_POINT_PRECISION* sum);


/*
//SUM_PHASE
void dataflow_sum_sinterm_stream(int dimension,FIXED_POINT_U_PRECISION sin_term_in[N_size][N_size],FIXED_POINT_U_PRECISION sin_term_out[N_size]);
void dataflow_sum_indiv_calc(FIXED_POINT_U_PRECISION sin_term_in[DATAFLOW_PARALLEL_ACCUM],FIXED_POINT_PRECISION* sum);
void dataflow_sum_dim_calc(FIXED_POINT_U_PRECISION sin_term_dim[N_size],FIXED_POINT_PRECISION *kuramoto_sum);
void dataflow_sum_phase_chunk(int dimension,FIXED_POINT_U_PRECISION sin_term_in[N_size][N_size],FIXED_POINT_PRECISION* sum);
void dataflow_sum_block_calc(FIXED_POINT_U_PRECISION sin_term_in[DATAFLOW_PARALLEL_F],FIXED_POINT_PRECISION* sum);//process DATAFLOW_PARALLEL_F elements
void dataflow_sum_group(int dimension,FIXED_POINT_U_PRECISION sin_term_in[N_size][N_size], FIXED_POINT_PRECISION* sum );
*/


//FINAL_PHASE
void dataflow_kuramoto_final_sum_calc(FIXED_POINT_PRECISION kuramoto_sum,FIXED_POINT_U_PRECISION phi_diff,FIXED_POINT_U_PRECISION* phi_result,FIXED_POINT_PRECISION C_sync_fix, FIXED_POINT_LOW C1_fix,int k_part);

//DELTA PHASE
void dataflow_delta_calc(FIXED_POINT_PRECISION phi_diff1,FIXED_POINT_PRECISION phi_diff2, FIXED_POINT_PRECISION *phi_delta);



//TOP DATAFLOW

void kuramoto_dataflow(int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_valid_row[J_VALID_SIZE],int total_block_num, FIXED_POINT_U_PRECISION phi_diff_remap[DATAFLOW_PARALLEL_F+1][N_size], int size,
FIXED_POINT_LOW C1_fix, FIXED_POINT_PRECISION C_sync_fix,int k_part,FIXED_POINT_U_PRECISION phi_next[N_size]);

void kuramoto_dataflow_top(double* data_in, int * J_data_in,double C_1, double C_sync, double ts, int iterations, int k_part,double* dwt_double, double* phi_d_out,int size, int length,int reset_J);//Jij*sin(phi_d+f(phi_d)



void top_dwt_in(double* dwt_double,FIXED_POINT_U_PRECISION dwt_temp[N_size*2],FIXED_POINT_U_PRECISION dwt_pre[N_size], int size);
void top_update_phi(FIXED_POINT_U_PRECISION phi_next[N_size], FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F+1][N_size], double* phi_dout,int size);
void top_update_feedb(FIXED_POINT_U_PRECISION feedb_pre[N_size],FIXED_POINT_U_PRECISION feedb_current[N_size],int size);
void top_loop_j_in(int* J_data_in, int total_length,int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_valid_row[J_VALID_SIZE],int* j_f_count,int reset_j);

void top_sde_update(FIXED_POINT_U_PRECISION sde_accum[4][N_size], FIXED_POINT_U_PRECISION phi_next[N_size],double* phi_d_out, FIXED_POINT_U_PRECISION ts,int size);
void top_update_phi_inter_iteration(FIXED_POINT_U_PRECISION phi_next[N_size],FIXED_POINT_U_PRECISION phi_next_inter_iteration[N_size], FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F+1][N_size], 
FIXED_POINT_U_PRECISION tstep, FIXED_POINT_U_PRECISION dwt[N_size*2], FIXED_POINT_U_PRECISION dwt_pre[N_size], FIXED_POINT_U_PRECISION sde_accum[N_size], int iteration, int size);
void top_phi_load(FIXED_POINT_U_PRECISION phi_next[N_size],FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F][N_size],int n_size);
void top_annealing_update(FIXED_POINT_LOW AC1_fix, FIXED_POINT_PRECISION AC_sync_fix,FIXED_POINT_PRECISION Kt_fix,
FIXED_POINT_PRECISION T_inverse_fix,FIXED_POINT_U_PRECISION tstep_fix,FIXED_POINT_PRECISION tstop_inverse_fix, FIXED_POINT_LOW* C1_fix,
FIXED_POINT_PRECISION* C_sync_fix,int iteration);