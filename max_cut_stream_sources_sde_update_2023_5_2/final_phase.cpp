#include "f_func.h"

void dataflow_final_phase(hls::stream <FIXED_POINT_PRECISION>&kuramoto_sum_single,hls::stream <FIXED_POINT_U_PRECISION>&phi_diff_final,hls::stream <FIXED_POINT_U_PRECISION>&kuramoto_final_single,FIXED_POINT_PRECISION C_sync_fix, FIXED_POINT_LOW C1_fix,int k_part,int size)
{
//#pragma HLS pipeline
for (int i = 0; i<size; i++)//consume
    {	

		#pragma HLS loop_tripcount min=N_size max=N_size avg=N_size
		#pragma HLS pipeline


		static FIXED_POINT_U_PRECISION k4[4];
		

		FIXED_POINT_PRECISION sum_f;
		FIXED_POINT_U_PRECISION final_f;
		FIXED_POINT_U_PRECISION phi_diff_f;



		kuramoto_sum_single.read(sum_f);
		phi_diff_final.read(phi_diff_f);	
	

		

		//printf(" phi final %f ",phi_diff_f.to_double());
	
        dataflow_kuramoto_final_sum_calc(sum_f,phi_diff_f,&final_f,C_sync_fix,C1_fix,k_part);

		kuramoto_final_single.write(final_f);


    }

}



void dataflow_kuramoto_final_sum_calc(FIXED_POINT_PRECISION kuramoto_sum,FIXED_POINT_U_PRECISION phi_diff,FIXED_POINT_U_PRECISION* phi_result,FIXED_POINT_PRECISION C_sync_fix, FIXED_POINT_LOW C1_fix,int k_part)
{

			FIXED_POINT_U_PRECISION injection;
			FIXED_POINT_INJECTION pre_injection;
			FIXED_POINT_INJECTION after_injection;
			FIXED_POINT_PRECISION injection_input;
			FIXED_POINT_U_PRECISION sin_var;
			FIXED_POINT_PRECISION after_injection_normalization;

			pre_injection=-C1_fix*kuramoto_sum;






		//self biasing not added Ac * h(c) *sin(pi*x(c)),shows zeros on matalab for h(c)



			injection_input=(FIXED_POINT)k_part*pi;//5 bit
			injection_input=phi_diff*injection_input;//5bitx3bit, 8-10 bit

			//sin_var=hls::sin(injection_input);
			cordic_sin_piecewise(injection_input,&sin_var);
			injection= -C_sync_fix*sin_var;

			after_injection=pre_injection+injection;//pi normalization, do we actually need every cycle?

			after_injection_normalization=(after_injection)*pi_inverse_u;

			*phi_result=after_injection_normalization;

			//printf(" %f", temp_phi_diff[i_7].to_double());
		
}
