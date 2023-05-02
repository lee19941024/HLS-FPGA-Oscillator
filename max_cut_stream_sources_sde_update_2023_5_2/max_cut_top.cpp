#include "f_func.h"


void kuramoto_dataflow_top(double* data_in, int * J_data_in,double C_1, double C_sync, double ts, int iterations, int k_part,double* dwt_double, double* phi_d_out,int size, int length,int reset_J)//Jij*sin(phi_d+f(phi_d)
{

	FIXED_POINT_LOW C1_fix=(FIXED_POINT_LOW)C_1;//coupling strength might change the bit precision req
	FIXED_POINT_PRECISION C_sync_fix=(FIXED_POINT_PRECISION)C_sync;	


	FIXED_POINT_U_PRECISION phi_next[N_size];



	#pragma HLS INTERFACE mode=m_axi port = data_in depth=N_size offset=slave bundle=gmem1 max_read_burst_length=BURST_d_in
	#pragma HLS INTERFACE mode=m_axi port = J_data_in depth=N_size*N_size*2 offset=slave bundle=gmem1 max_read_burst_length=BURST_j_in
	#pragma HLS INTERFACE mode=m_axi port = phi_d_out depth=N_size*8000 offset=slave bundle=gmem1 max_write_burst_length=BURST_d_out
	#pragma HLS INTERFACE mode=m_axi port = dwt_double depth=N_size*2*8000 offset=slave bundle=gmem0 max_read_burst_length=BURST_dwt_in

	#pragma HLS INTERFACE mode=s_axilite port=C_1
	#pragma HLS INTERFACE mode=s_axilite port=C_sync
	#pragma HLS INTERFACE mode=s_axilite port=ts
	#pragma HLS INTERFACE mode=s_axilite port=iterations
	#pragma HLS INTERFACE mode=s_axilite port=k_part
	#pragma HLS INTERFACE mode=s_axilite port=return


	#pragma HLS PIPELINE off
		


	
	#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_d_out type=cyclic variable=phi_next




int total_length=length*2*2;//2 per row, then triangle duplicate
static int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F];
static int j_valid_row[J_VALID_SIZE];
static int j_f_count;
	#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=j_valid//DONT RESHAPE!!

	#pragma HLS bind_storage variable=j_valid type=RAM_T2P impl=URAM

top_loop_j_in(J_data_in, total_length,j_valid,j_valid_row,&j_f_count,reset_J);


/*
printf("total block %d \n",j_f_count);
int j_print_row=-1;
for (int i = 0; i < j_f_count; i++)
{

	if(j_valid_row[i]!=j_print_row)
	{
		j_print_row++;
		printf(" \n");
		printf("%d row  ",j_print_row+1);
	}
	for ( int j = 0; j < DATAFLOW_PARALLEL_F; j++)
	{	
		int temp_j=j_valid[i][j];
		
		
			printf("%d ",temp_j+1);

		
	}

}
printf("\n");
*/




/*
for (int i = 0; i < j_f_count; i++)
{
	for ( int j = 0; j < DATAFLOW_PARALLEL_F; j++)
	{
		int temp_j=j_valid[i][j];
		
		
		printf("%d ",temp_j);
	}
	printf(" \n ");	
}

*/





	loop_phi_in:	for(int i_1=0;i_1<size;i_1+=BURST_d_in)
		{
#pragma HLS loop_tripcount min=N_size/BURST_d_in max=N_size/BURST_d_in avg=N_size/BURST_d_in

		FIXED_POINT_U_PRECISION phi_in_temp_array[BURST_d_in];

		#pragma HLS PIPELINE II=BURST_d_in

			for (int j = 0; j < BURST_d_in; j++)
			{
				phi_in_temp_array[j]=(FIXED_POINT_U_PRECISION)data_in[i_1+j];

			}
			
			for (int j = 0; j < BURST_d_in; j++)
			{

				phi_next[i_1+j]=phi_in_temp_array[j];
				
			}
					

		}





	
	FIXED_POINT_U_PRECISION tstep_U=FIXED_POINT_U_PRECISION(ts);
	


	LOOP_ITERATIVE:for (int i = 0; i < iterations; i++)
	{
		FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F+1][N_size];
		FIXED_POINT_U_PRECISION sde_accum[4][N_size];
		FIXED_POINT_U_PRECISION phi_next_sde[N_size];


		FIXED_POINT_U_PRECISION dwt_temp[N_size*2]={0};
		FIXED_POINT_U_PRECISION dwt_pre[N_size]={0};

		#pragma HLS ARRAY_PARTITION dim=1 factor=BURST_dwt_in type=cyclic variable=dwt_temp
		#pragma HLS ARRAY_PARTITION dim=1 factor=BURST_dwt_in type=cyclic variable=dwt_pre


		#pragma HLS pipeline off//no pipeline
	

		#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=phi_current_remap
		#pragma HLS bind_storage variable=phi_current_remap type=RAM_T2P impl=URAM
    	#pragma HLS ARRAY_RESHAPE dim=2 factor=2 type=cyclic variable=phi_current_remap
	
		#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_d_out type=cyclic variable=phi_next_sde
		#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=sde_accum
		#pragma HLS ARRAY_PARTITION dim=2 factor=BURST_d_out type=cyclic variable=sde_accum


		top_phi_load(phi_next,phi_current_remap,size);
		top_dwt_in(&dwt_double[i*N_size*2],dwt_temp,dwt_pre,size);

		LOOP_ITERATIVE_SDE:for (int j = 0; j < 4; j++)
		{
			kuramoto_dataflow(j_valid,j_valid_row,j_f_count,phi_current_remap,size,C1_fix,C_sync_fix,k_part,phi_next_sde);
			top_update_phi_inter_iteration(phi_next,phi_next_sde,phi_current_remap,tstep_U, dwt_temp,dwt_pre,sde_accum[j],j,size);		
		}

		top_sde_update(sde_accum,phi_next,  &phi_d_out[i*N_size],tstep_U, size);

	
	}
	

}

void kuramoto_dataflow(int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_valid_row[J_VALID_SIZE],int total_block_num, FIXED_POINT_U_PRECISION phi_diff_remap[DATAFLOW_PARALLEL_F+1][N_size], int size,
FIXED_POINT_LOW C1_fix, FIXED_POINT_PRECISION C_sync_fix,int k_part,FIXED_POINT_U_PRECISION phi_next[N_size])
{



	hls::stream <phi_f_block,5>phi_delta_struct_F;
	hls::stream <FIXED_POINT_U_PRECISION,1000> phi_diff_single;

	hls::stream <FIXED_POINT_PRECISION,5> sin_partial_single;
	hls::stream <FIXED_POINT_PRECISION,5> kuramoto_single;

	hls::stream <FIXED_POINT_U_PRECISION,5> kuramoto_final_single;


	hls::stream <int,5>current_row_single;
	hls::stream <int,5>current_row_partial_single;
    #pragma HLS DATAFLOW


	dataflow_delta_phase(total_block_num,j_valid,j_valid_row,current_row_single,phi_diff_remap,phi_delta_struct_F,phi_diff_single);
    dataflow_pre_sum_phase(current_row_single,current_row_partial_single,phi_delta_struct_F,sin_partial_single,k_part,total_block_num);//consume DATAFLOW_PARALLEL_F data and output same amount
	dataflow_partial_sum_phase(sin_partial_single,current_row_partial_single, kuramoto_single,total_block_num);
	dataflow_final_phase(kuramoto_single,phi_diff_single,kuramoto_final_single,C_sync_fix,C1_fix,k_part,size);
	dataflow_sde_accumulate(kuramoto_final_single,phi_next,size);
	
}




void top_dwt_in(double* dwt_double,FIXED_POINT_U_PRECISION dwt_temp[N_size*2],FIXED_POINT_U_PRECISION dwt_pre[N_size], int size)
{

	#pragma HLS pipeline off
	#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_dwt_in type=cyclic variable=dwt_temp
	#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_dwt_in type=cyclic variable=dwt_pre


	LOOP_DWT_PRE_UPDATE:for (int i = 0; i < size; i++)
	{
			#pragma HLS loop_tripcount min=N_size max=N_size avg=N_size

		#pragma HLS UNROLL factor=BURST_dwt_in
		#pragma HLS PIPELINE
		dwt_pre[i]=dwt_temp[N_size+i];
	}
	


	LOOP_DWT_NEW_IN1:for (int i = 0; i < size; i+=BURST_dwt_in)
	{
		#pragma HLS loop_tripcount min=N_size/BURST_d_in max=N_size/BURST_d_in avg=N_size/BURST_d_in


		FIXED_POINT_U_PRECISION dwt_buffer[BURST_dwt_in];
		#pragma HLS ARRAY_PARTITION dim=0 type=complete variable=dwt_buffer

		#pragma HLS pipeline
		LOOP_DWT_DOUBLE:for (int j = 0; j < BURST_dwt_in; j++)
		{
			dwt_buffer[j]=(FIXED_POINT_U_PRECISION)dwt_double[i+j];

		}

		LOOP_DWT_DOUBLE_TEMP:for (int j = 0; j < BURST_dwt_in; j++)
		{
			dwt_temp[i+j]=dwt_buffer[j];

		}

	}

	LOOP_DWT_NEW_IN2:for (int i = 0; i < size; i+=BURST_dwt_in)
	{
		#pragma HLS loop_tripcount min=N_size/BURST_d_in max=N_size/BURST_d_in avg=N_size/BURST_d_in


		FIXED_POINT_U_PRECISION dwt_buffer[BURST_dwt_in];
		#pragma HLS ARRAY_PARTITION dim=0 type=complete variable=dwt_buffer

		#pragma HLS pipeline
		LOOP_DWT_DOUBLE_2:for (int j = 0; j < BURST_dwt_in; j++)
		{
			dwt_buffer[j]=(FIXED_POINT_U_PRECISION)dwt_double[N_size+i+j];

		}

		LOOP_DWT_DOUBLE_TEMP_2:for (int j = 0; j < BURST_dwt_in; j++)
		{
			dwt_temp[N_size+i+j]=dwt_buffer[j];

		}

	}

	

}

void top_update_phi(FIXED_POINT_U_PRECISION phi_next[N_size],FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F+1][N_size], double* phi_dout,int size)
{
	#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_d_out type=cyclic variable=phi_next


	
	#pragma HLS ARRAY_RESHAPE dim=2 factor=2 type=cyclic variable=phi_current_remap
	#pragma HLS ARRAY_PARTITION dim=1  type=complete variable=phi_current_remap
	#pragma HLS bind_storage variable=phi_current_remap type=RAM_T2P impl=URAM



#pragma hls pipeline off
		

	for (int i = 0; i < size; i+=BURST_d_out)
	{
#pragma HLS pipeline II=BURST_d_out
#pragma HLS loop_tripcount min=N_size/BURST_d_in max=N_size/BURST_d_in avg=N_size/BURST_d_in

		FIXED_POINT_U_PRECISION phi_buffer[BURST_d_out];
		#pragma HLS ARRAY_PARTITION dim=0 type=complete variable=phi_buffer
		


		for (int j = 0; j < BURST_d_out; j++)
		{
		#pragma HLS UNROLL
			phi_buffer[j]=phi_next[i+j];
		}

		for (int j = 0; j < BURST_d_out; j++)
		{
		#pragma HLS UNROLL

			//phi_current[i+j]=phi_buffer[j];
			phi_dout[i+j]=phi_buffer[j].to_double();

			for (int k = 0; k < DATAFLOW_PARALLEL_F+1; k++)
			{		
				#pragma HLS UNROLL

				phi_current_remap[k][i+j]=phi_buffer[j];
			}
			


		}
		
	}


}


void top_update_phi_inter_iteration(FIXED_POINT_U_PRECISION phi_next[N_size],FIXED_POINT_U_PRECISION phi_next_inter_iteration[N_size], FIXED_POINT_U_PRECISION phi_current_remap[DATAFLOW_PARALLEL_F+1][N_size], 
FIXED_POINT_U_PRECISION tstep, FIXED_POINT_U_PRECISION dwt[N_size*2], FIXED_POINT_U_PRECISION dwt_pre[N_size], FIXED_POINT_U_PRECISION sde_accum[N_size], int iteration, int size)
{
#pragma HLS array_partition variable=dwt type=cyclic factor=BURST_d_out
#pragma HLS array_partition variable=dwt_pre type=cyclic factor=BURST_d_out
#pragma HLS array_partition variable=sde_accum type=cyclic factor=BURST_d_out

#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_d_out type=cyclic variable=phi_next
#pragma HLS ARRAY_PARTITION dim=0 factor=BURST_d_out type=cyclic variable=phi_next_inter_iteration


	
#pragma HLS ARRAY_RESHAPE dim=2 factor=2 type=cyclic variable=phi_current_remap
#pragma HLS ARRAY_PARTITION dim=1  type=complete variable=phi_current_remap
#pragma HLS bind_storage variable=phi_current_remap type=RAM_T2P impl=URAM



LOOP_INTER_SDE:	for (int i = 0; i < size; i++)
	{
		#pragma HLS loop_tripcount min=N_size max=N_size avg=N_size
		#pragma HLS dependence variable=sde_accum type=inter false
		#pragma HLS UNROLL factor=2
		#pragma HLS pipeline
		FIXED_POINT_U_PRECISION temp_dwt;
		FIXED_POINT_U_PRECISION accum_factor;
		FIXED_POINT_U_PRECISION	sde_factor;
		FIXED_POINT_U_PRECISION sde_next;
		switch (iteration)
		{
		case 0:
			temp_dwt=dwt_pre[i];
			accum_factor=1;
			sde_factor=0.5;
			break;
		case 1:
			accum_factor=2;
			sde_factor=0.5;
			temp_dwt=dwt[i];
			break;
		case 2:
			accum_factor=2;
			sde_factor=1;
			temp_dwt=dwt[i];
			break;
		case 3:
			accum_factor=1;
			sde_factor=0;
			temp_dwt=dwt[N_size+i];
			break;
		default:
			accum_factor=0;
			sde_factor=0;
			temp_dwt=0;
			break;
		}

		#pragma HLS PIPELINE

		FIXED_POINT_U_PRECISION k_iteration=phi_next_inter_iteration[i]+temp_dwt;
		FIXED_POINT_U_PRECISION temp_phi_next= phi_next[i]+sde_factor*tstep*k_iteration;

		sde_accum[i]=k_iteration*accum_factor;

		for (int k = 0; k<DATAFLOW_PARALLEL_F+1; k++)
		{
			#pragma HLS UNROLL
			phi_current_remap[k][i]=temp_phi_next;

		}
		
	}
	

}




/*

void top_j_read()
{

	int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F];
	int j_valid_count[N_size];

	#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=j_valid

	loop_J_in:  

	int j_index=0;
	int j_count=0;
	int j_total=0;
	
	
	for (int i_2=0;i_2<size*size;i_2+=BURST_j_in)
		{	
			#pragma HLS dependence variable=j_valid_count type=inter false
			#pragma HLS dependence variable=j_valid_count type=intra false
			#pragma HLS loop_tripcount min=N_size*N_size/BURST_j_in max=N_size*N_size/BURST_j_in avg=N_size*N_size/BURST_j_in

		
			#pragma HLS loop_tripcount min=N_size/BURST_j_in max=N_size/BURST_j_in avg=N_size/BURST_j_in
				#pragma HLS dependence variable=j_valid type=inter false
				#pragma HLS dependence variable=j_valid type=intra false
			#pragma HLS PIPELINE
			
	
			
				FIXED_POINT_J J_buffer[BURST_j_in];
				FIXED_POINT_J J_valid_buffer[BURST_j_in];
				#pragma HLS ARRAY_PARTITION dim=0  type=complete variable=J_buffer

				for (int k = 0; k < BURST_j_in; k++)
				{
					#pragma HLS UNROLL
					J_buffer[k]=(FIXED_POINT_J)J_data_in[i_2+k];
				}
				for (int k = 0; k < BURST_j_in; k++)
				{
					if(J_buffer[k]!=0)
					{
						
					#pragma HLS UNROLL


						J_valid_buffer[k]=i_2%N_size+k;
						j_index++;
						j_count++;
					}
					else
					{
						j_valid_buffer[k]=-1;
					}

			

		}

}

}

*/





void top_loop_j_in(int* J_data_in, int total_length,int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_valid_row[J_VALID_SIZE],int* j_f_count,int reset_j)
{
	int current_row;
	int current_j_valid_row=0;
	int j_count=0;

	current_row=0;
	#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=j_valid

	#pragma HLS bind_storage variable=j_valid type=RAM_T2P impl=URAM

	#pragma HLS pipeline off


if(reset_j==1)
{

	loop_j_valid_reset:for (int i = 0; i < J_VALID_SIZE; i++)
	{
		#pragma HLS pipeline
		for (int j = 0; j < DATAFLOW_PARALLEL_F; j++)
		{
			#pragma HLS UNROLL
			j_valid[i][j]=-1;
		}
	}
	


	loop_j_in: for (int i = 0; i < total_length; i+=DATAFLOW_PARALLEL_F*2)
	{
		#pragma HLS loop_tripcount min=J_VALID_SIZE max=J_VALID_SIZE avg=J_VALID_SIZE
		
			for (int j = 0; j <DATAFLOW_PARALLEL_F; j++)
			{
				#pragma HLS pipeline 		
				#pragma HLS dependence variable=j_valid type=intra false
				#pragma HLS dependence variable=j_valid_row type=intra false


					int temp_row=J_data_in[i+j*2]-1;
					int temp_j_valid=J_data_in[i+j*2+1]-1;

					if(i+j*2<total_length)
					{
					


					if(temp_row!=current_row)
					{	

						j_valid[current_j_valid_row][j]=-1;
						j_valid_row[current_j_valid_row]=current_row;
						current_row++;

						j_valid[current_j_valid_row+1][j]=temp_j_valid;
						current_j_valid_row++;

						j_valid_row[current_j_valid_row]=current_row;
					}
					else
					{
						j_valid[current_j_valid_row][j]=temp_j_valid;
						j_valid_row[current_j_valid_row]=current_row;

					}

					}

				}
		current_j_valid_row++;


			}

	

	*j_f_count=current_j_valid_row;

}

}



void top_annealing_update(FIXED_POINT_LOW AC1_fix, FIXED_POINT_PRECISION AC_sync_fix,FIXED_POINT_PRECISION Kt_fix,
FIXED_POINT_PRECISION T_inverse_fix,FIXED_POINT_U_PRECISION tstep_fix,FIXED_POINT_PRECISION tstop_inverse_fix, FIXED_POINT_LOW* C1_fix,
FIXED_POINT_PRECISION* C_sync_fix,int iteration)
{




*C1_fix=(FIXED_POINT_LOW)1+(AC1_fix-(FIXED_POINT_LOW)1)*tstop_inverse_fix*(iteration-1)*tstep_fix;

FIXED_POINT_LOW C1_temp=(AC1_fix-(FIXED_POINT_LOW)1)*tstop_inverse_fix;
printf(" C1_temp %f \n", C1_temp.to_double());



FIXED_POINT_PRECISION tanh_input;
FIXED_POINT_PRECISION tanh_result;
FIXED_POINT_PRECISION cos_input;
FIXED_POINT_PRECISION cos_result;

cos_input=(FIXED_POINT_PRECISION)2*pi*(FIXED_POINT_PRECISION)(iteration-1)*tstep_fix*T_inverse_fix;
cos_result=hls::cos(cos_input);
tanh_input=Kt_fix*cos_result;
tanh_result=hls::tanh(tanh_input);
*C_sync_fix=(FIXED_POINT_LOW)1+(FIXED_POINT_LOW)2*tanh_result;

}