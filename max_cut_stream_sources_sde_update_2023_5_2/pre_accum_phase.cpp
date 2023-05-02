#include "f_func.h"

void dataflow_pre_sum_phase(hls::stream<int>&current_row_pre_accum_single, hls::stream<int>&current_row_partial_sum_single,hls::stream <phi_f_block>&phi_delta_struct_F,hls::stream <FIXED_POINT_PRECISION>&sin_term_single,int k_part,int total_block_num)
{

	#pragma HLS pipeline off

	loop_pre_sum_phase:for (int i = 0; i < total_block_num; i++)
	{
		#pragma HLS loop_tripcount min=N_size max=N_size avg=N_size
		#pragma HLS pipeline
		
		#pragma HLS loop_tripcount min=N_size*N_size/DATAFLOW_PARALLEL_F max=N_size*N_size/DATAFLOW_PARALLEL_F avg=N_size*N_size/DATAFLOW_PARALLEL_F

		FIXED_POINT_PRECISION sin_term_partial;
		phi_f_block phi_delta_dataflow_struct;



		int current_row;
		current_row_pre_accum_single.read(current_row);
		current_row_partial_sum_single.write(current_row);


    	FIXED_POINT_U_PRECISION f_result_dataflow[DATAFLOW_PARALLEL_F];
    	FIXED_POINT_U_PRECISION phi_delta_dataflow_sync[DATAFLOW_PARALLEL_F];
		FIXED_POINT_U_PRECISION sin_term_dataflow[DATAFLOW_PARALLEL_F];


    	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=f_result_dataflow
    	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=phi_delta_dataflow_sync

	
		phi_delta_struct_F.read(phi_delta_dataflow_struct);
			
    	dataflow_f_func_calc(phi_delta_dataflow_struct.phi_delta ,f_result_dataflow,phi_delta_dataflow_sync,k_part);

    	dataflow_kuramor_pre_accum_calc(phi_delta_dataflow_sync,f_result_dataflow,sin_term_dataflow );


		

		dataflow_sin_accumulate_calc( sin_term_dataflow,&sin_term_partial);//process DATAFLOW_PRALLEL_ACCUM elements	
		


		sin_term_single.write(sin_term_partial);
		//printf("%f ",sin_term_partial.to_double());
		}

	}
	  





void dataflow_sin_accum_phase(hls::stream<FIXED_POINT_U_PRECISION> sin_term_in_block[DATAFLOW_PARALLEL_F],hls::stream<FIXED_POINT_PRECISION>&sin_term_partial_single,int size)
{


		for (int i = 0; i < size; i++)
		{
			#pragma HLS loop_tripcount min=N_size max=N_size avg=N_size
			#pragma HLS pipeline

			for (int j = 0; j<N_size; j+=DATAFLOW_PARALLEL_F)
			{
				#pragma HLS pipeline
				FIXED_POINT_PRECISION temp_sum=0;;
				#pragma HLS loop_tripcount min=N_size/DATAFLOW_PARALLEL_F max=N_size/DATAFLOW_PARALLEL_F avg=N_size/DATAFLOW_PARALLEL_F
				for (int k = 0; k < DATAFLOW_PARALLEL_F; k++)
				{
				#pragma HLS UNROLL

					FIXED_POINT_U_PRECISION temp_sin_term;
					sin_term_in_block[k].read(temp_sin_term);
					temp_sum+=temp_sin_term;
				}
				sin_term_partial_single.write(temp_sum);

			}
			
		}

}



void dataflow_f_func_calc(FIXED_POINT_U_PRECISION phi_delta_in[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION f_result_out[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION phi_delta_sync_out[DATAFLOW_PARALLEL_F],int k_part)
{
	#pragma HLS pipeline
	//printf("f func \n");
	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=phi_delta_in
	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=f_result_out
	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=phi_delta_sync_out



	dataflow_loop_func_fout:for(int i_51=0;i_51<DATAFLOW_PARALLEL_F;i_51++)
	{
		//#pragma HLS dependence variable=phi_delta_in type=inter false
      	//#pragma HLS dependence variable=phi_delta_in type=intra false
		//#pragma HLS dependence variable=f_result_out type=inter false
      	//#pragma HLS dependence variable=f_result_out type=intra false

		//#pragma HLS dependence variable=phi_delta_sync_out type=inter false
      	//#pragma HLS dependence variable=phi_delta_sync_out type=intra false


		#pragma HLS UNROLL factor=DATAFLOW_PARALLEL_F skip_exit_check

		FIXED_POINT_U_PRECISION f_func_input;
		FIXED_POINT_U_PRECISION f_func_output;

		f_func_input= phi_delta_in[i_51];
		phi_delta_sync_out[i_51]=f_func_input;


		f_func_lut(f_func_input,&f_func_output,k_part);
		//f_func_accumulate( f_func_input,&f_func_output);


		f_result_out[i_51]=f_func_output;
		//if(i_51%N_size==0) printf("\n");

		//printf(" delta %f ",f_result_out[i_51].to_double());


	}

}

void dataflow_kuramor_pre_accum_calc(FIXED_POINT_U_PRECISION phi_delta_in[DATAFLOW_PARALLEL_F],FIXED_POINT_U_PRECISION f_result_in[DATAFLOW_PARALLEL_F], FIXED_POINT_U_PRECISION sin_term_out[DATAFLOW_PARALLEL_F])
{
#pragma HLS pipeline


	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=phi_delta_in
	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=f_result_in
	#pragma HLS ARRAY_PARTITION dim=1 factor=DATAFLOW_PARALLEL_F type=cyclic variable=sin_term_out


	dataflow_loop_kuramoto_pre_accum: for(int i_6=0;i_6<DATAFLOW_PARALLEL_F;i_6++)
		{
			//#pragma HLS dependence variable=phi_delta_in type=inter false
      		//#pragma HLS dependence variable=phi_delta_in type=intra false

			//#pragma HLS dependence variable=f_result_in type=inter false
      		//#pragma HLS dependence variable=f_result_in type=intra false

			//#pragma HLS dependence variable=sin_term_out type=inter false
      		//#pragma HLS dependence variable=sin_term_out type=intra false

			//#pragma HLS dependence variable=J_matrix_in type=inter false
      		//#pragma HLS dependence variable=J_matrix_in type=intra false


			#pragma HLS UNROLL factor=DATAFLOW_PARALLEL_F skip_exit_check

				//if(i_6%N_size==0) printf("\n");

				
				sin_calc(f_result_in[i_6],phi_delta_in[i_6],&sin_term_out[i_6]);
				
		}
}


void sin_calc(FIXED_POINT_U_PRECISION f_result,FIXED_POINT_U_PRECISION phi_delta,FIXED_POINT_U_PRECISION* result_out)//need optimization
{
#pragma HLS pipeline
	FIXED_POINT_PRECISION sin_preterm;
	FIXED_POINT_U_PRECISION sin_var;
	FIXED_POINT_U_PRECISION result;

	FIXED_POINT_U_PRECISION expo;
	FIXED_POINT_U_PRECISION sigmoid;

					sin_preterm=pi*phi_delta+f_result;

					//printf("sinterm pi %f ",sin_term[i_6].to_double());

					cordic_sin_piecewise(sin_preterm,&sin_var);//difference from slide ? line 212 matlab
					//sin_var=hls::sin(sin_preterm);

					//printf(" sinterm_sin %f ",sin_term[i_6].to_double());


					sin_var=(FIXED_POINT_U_PRECISION)0.1*(sin_var);//difference from slide ? line 212 matlab strength of square wave k=0.1
					//printf(" sinterm_Ksin %f ",sin_term[i_6].to_double());


					/*
					result=hls::tanh(sin_var);//why tanh always positive, xilinx ??????,better activation functionn


					if(sin_var>0)
					{
						if (result<0) result=-result;

					}
					else
					{
						if (result>0) result=-result;
					}

					//printf(" sinterm tanh %s ",result.to_string(10).c_str());
					*/


					sin_var=-2*sin_var;
					expo=hls::exp(sin_var);

					sigmoid=(FIXED_POINT_U_PRECISION)1/((FIXED_POINT_U_PRECISION)1+expo);

					result=(sigmoid-(FIXED_POINT_U_PRECISION)0.5)*(FIXED_POINT_U_PRECISION)2;

					

					*result_out=(-1)*result;//J is always negative, so..

					//*result_out=result;

}

void dataflow_sin_accumulate_calc(FIXED_POINT_U_PRECISION sin_term_in[DATAFLOW_PARALLEL_F],FIXED_POINT_PRECISION* sum)//process DATAFLOW_PRALLEL_ACCUM elements
{	
	#pragma HLS PIPELINE
	FIXED_POINT_PRECISION tempsum;
	//printf("sin term in \n ");

    tempsum=0;
  loop_dataflow_indivi_calc: for (int i = 0; i<DATAFLOW_PARALLEL_F ; i++)
    {
		#pragma HLS pipeline
        #pragma HLS UNROLL factor=DATAFLOW_PARALLEL_F skip_exit_check
        tempsum+=sin_term_in[i];

		//printf(" %f ", sin_term_in[i].to_double());

    }

  

  *sum=tempsum;
}
