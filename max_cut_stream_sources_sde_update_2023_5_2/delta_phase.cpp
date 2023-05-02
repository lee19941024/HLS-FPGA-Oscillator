#include "f_func.h"



void dataflow_delta_phase(int total_block_num,int j_valid[J_VALID_SIZE][DATAFLOW_PARALLEL_F],int j_row[J_VALID_SIZE],hls::stream <int>&current_row_single,FIXED_POINT_U_PRECISION phi_diff_in_remap[DATAFLOW_PARALLEL_F+1][N_size], 
hls::stream<phi_f_block>&phi_delta_struct_F,hls::stream <FIXED_POINT_U_PRECISION>&phi_diff_single)
{
	//#pragma HLS ARRAY_PARTITION dim=0 factor=DATAFLOW_PARALLEL_D type=cyclic variable=phi_delta_out
	
	#pragma HLS ARRAY_PARTITION dim=2 type=complete variable=j_valid
	#pragma HLS bind_storage variable=j_valid type=RAM_T2P impl=URAM

	#pragma HLS ARRAY_RESHAPE dim=2 factor=2 type=cyclic variable=phi_diff_in_remap
	#pragma HLS array_partition variable =phi_diff_in_remap type=complete dim=1// factor=DATAFLOW_PARALLEL_F
	#pragma HLS bind_storage variable=phi_diff_in_remap type=RAM_T2P impl=URAM



	#pragma HLS pipeline off
	#pragma HLS inline off
	int j_valid_offset=0;
	int temp_phi_block_num;

	int current_row;
	current_row=-1;//initial for first phi
	FIXED_POINT_U_PRECISION phi_diff1;
				
	LOOP_DELTA_TOP:for (int i = 0; i < total_block_num; i++)
		{
#pragma HLS latency min=80


			#pragma HLS dependence variable=j_valid type=inter false
			#pragma HLS dependence variable=j_valid type=intra false

			#pragma HLS dependence variable=phi_diff_in_remap type=intra false

			#pragma HLS loop_tripcount min=J_VALID_SIZE max=J_VALID_SIZE avg=J_VALID_SIZE

			#pragma HLS pipeline 
			
			if (current_row!=j_row[i])
			{
				current_row++;

				//uint32_t dim1;
				//uint32_t dim2;
				//uint32_t current_row32=current_row;
				//remap_address(current_row32,&dim1,&dim2);
				phi_diff1=phi_diff_in_remap[DATAFLOW_PARALLEL_F][current_row];
				//phi_diff1=phi_diff_in[current_row];
				phi_diff_single.write(phi_diff1);
		
				//printf( "%d delta row ",current_row);
				//printf(" \n");
			}
			
			phi_f_block phi_f;
			


			for (int j = 0; j < DATAFLOW_PARALLEL_F; j++)
				{
					//#pragma HLS dependence variable=phi_diff_in type=intra false
					//#pragma HLS dependence variable=phi_diff_in type=inter false

					#pragma HLS UNROLL
					FIXED_POINT_U_PRECISION phi_temp;
					int phi_index=j_valid[i][j];
					//printf(" %d ", phi_index+1);
					if(phi_index>=0)
					{
						if(phi_index==current_row)
						{
							phi_temp=0;
						}
						else
						{
							FIXED_POINT_U_PRECISION phi_diff_2;
							//uint32_t dim1;
							//uint32_t dim2;
							//uint32_t phi_index_32=phi_index;
							//remap_address(phi_index_32,&dim1,&dim2);

							//phi_diff_2=phi_diff_in[phi_index];
							phi_diff_2=phi_diff_in_remap[j][phi_index];

							//printf("%f ", phi_diff_2.to_double());
							phi_temp=phi_diff1-phi_diff_2;

						}	

					}
					else
					{
						phi_temp=0;
					}
					
					phi_f.phi_delta[j]=phi_temp;
			
				}
				//printf(" \n");
				
				phi_delta_struct_F.write(phi_f);
				current_row_single.write(current_row);											

		}
}
