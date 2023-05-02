#include"f_func.h"

void dataflow_partial_sum_phase(hls::stream <FIXED_POINT_PRECISION>&sin_partial_single, hls::stream<int>&current_row_single, hls::stream <FIXED_POINT_PRECISION> &kuramoto_single,int total_block_num)
{


     int current_row;
     current_row=0;
     FIXED_POINT_PRECISION temp_sum;
     temp_sum=0;

        for (int i = 0; i < total_block_num; i++)
        {
            #pragma HLS loop_tripcount min=N_size*N_size/DATAFLOW_PARALLEL_F max=N_size*N_size/DATAFLOW_PARALLEL_F avg=N_size*N_size/DATAFLOW_PARALLEL_F
            #pragma HLS pipeline

		        int temp_row;
                FIXED_POINT_PRECISION temp_sin_partial;
                current_row_single.read(temp_row);
                sin_partial_single.read(temp_sin_partial);
                

                if(temp_row==current_row)
                {
                    temp_sum+=temp_sin_partial;			
                }
                else
                {
                    //printf("%d temp sum %f \n",current_row+1,temp_sum.to_double());
                    kuramoto_single.write(temp_sum);
                    temp_sum=temp_sin_partial;
                    current_row++;   

                }
        }
        
        //printf("%d temp sum %f \n",current_row+1,temp_sum.to_double());

        kuramoto_single.write(temp_sum);
}






