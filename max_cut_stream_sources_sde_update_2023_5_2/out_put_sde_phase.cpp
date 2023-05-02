#include "f_func.h"


void dataflow_sde_accumulate(hls::stream <FIXED_POINT_U_PRECISION>&kuramoto_final_single,FIXED_POINT_U_PRECISION phi_next[N_size],int size)
{
	


	for (int i = 0; i < size; i++)
	{	
		#pragma HLS pipeline

		FIXED_POINT_U_PRECISION temp_kuramoto_final;

		kuramoto_final_single.read(temp_kuramoto_final);		
		
		phi_next[i]=temp_kuramoto_final;
	}

}


void top_sde_update(FIXED_POINT_U_PRECISION sde_accum[4][N_size], FIXED_POINT_U_PRECISION phi_next[N_size],double* phi_d_out, FIXED_POINT_U_PRECISION ts,int size)
{


LOOP_TOP_SDE:for (int i = 0; i < size; i+=BURST_d_out)
{
#pragma HLS loop_tripcount min=N_size/BURST_d_in max=N_size/BURST_d_in avg=N_size/BURST_d_in
#pragma HLS pipeline
	FIXED_POINT_U_PRECISION phi_out_buffer[BURST_d_out];
	for (int j = 0; j < BURST_d_out; j++)
	{
		FIXED_POINT_U_PRECISION temp_sde_sum;
		temp_sde_sum=sde_accum[0][i+j]+sde_accum[1][i+j]+sde_accum[2][i+j]+sde_accum[3][i+j];
		phi_out_buffer[j]=phi_next[i+j]+(FIXED_POINT_U_PRECISION)0.166666667*ts*temp_sde_sum;

	}
	
	
	for (int j = 0; j < BURST_d_out; j++)
	{
		phi_d_out[i+j]=phi_out_buffer[j];
		phi_next[i+j]=phi_out_buffer[j];
	}
	

}



}



void top_phi_load(FIXED_POINT_U_PRECISION phi_next[N_size],FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F][N_size],int n_size)
{

	LOOP_PHI_LOAD:for (int i = 0; i < n_size; i++)
	{
#pragma HLS LOOP_TRIPCOUNT min=N_size max=N_size avg=N_size
		#pragma HLS pipeline
		#pragma HLS UNROLL factor=2
		for (int j = 0; j < DATAFLOW_PARALLEL_F; j++)
		{
			
			phi_current_remap[j][i]=phi_next[i];
			
		}
		//printf(" %f ", phi_next[i].to_double());
	}
	
//printf("\n");
}
