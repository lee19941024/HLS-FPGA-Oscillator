#include "f_func.h"
/*
void f_func_accumulate(FIXED_POINT_U_PRECISION theta, FIXED_POINT_PRECISION* result)
{
#pragma HLS INLINE OFF
#pragma HLS pipeline



	FIXED_POINT_PRECISION sum=0;
	FIXED_POINT_PRECISION current_sum[upper_K_int-1];
	FIXED_POINT_LOW lower_k;
	//FIXED_POINT upper_K=5;


f_loop:{

//#pragma HLS LATENCY max=200

	for(int i=1;i<upper_K_int;i++)
	{
		//#pragma HLS pipeline rewind
		#pragma HLS inline off

		lower_k=(FIXED_POINT_LOW)i;
		f_func_eq(lower_k,theta,current_sum+i-1);

	}
}
sum_loop:{

	for (int j=0;j<upper_K_int-1;j++)
	{
		//#pragma HLS PIPELINE rewind

		sum=sum+current_sum[j];
	}
}


	//printf("expo f %f \n",sum.to_double());
	//printf("\n");
	*result=sum;

}


void f_func_eq(FIXED_POINT_LOW lower_k, FIXED_POINT_U_PRECISION theta, FIXED_POINT_PRECISION* result)
{

//#pragma HLS DATAFLOW

	FIXED_POINT_PRECISION common_term1;//(2k-1) max 7 
	FIXED_POINT_PRECISION common_term2;//(2k/K) max 2
	FIXED_POINT_PRECISION gaussmf_1,gaussmf_2;
	FIXED_POINT_PRECISION gauss_coeff1,gauss_coeff2;
	FIXED_POINT_PRECISION gaussian_mean1,gaussian_mean2;
	FIXED_POINT_PRECISION eq1,eq2;
	FIXED_POINT_PRECISION result_f;

	FIXED_POINT_LOW constant_two=2;


	common_term1=(constant_two*lower_k-1); //(2k-1)
	common_term2=(constant_two*lower_k)/upper_K;//2k/K


	gaussian_mean1=common_term2;// 7*pi,different from slide, no multiply to pi?
	gaussian_mean2=-common_term2;



	#pragma HLS inline off
	f_func_gaussmf(theta,gaussian_mean1,&gaussmf_1);
	#pragma HLS inline off

	f_func_gaussmf(theta,gaussian_mean2,&gaussmf_2);//return value within 0-1


	gauss_coeff1=common_term1-common_term2;
	gauss_coeff2=common_term2-common_term1;

	eq1=gauss_coeff1*gaussmf_1;
	eq2=gauss_coeff2*gaussmf_2;


	result_f=eq1+eq2;
	result_f=result_f*pi;





	*result= result_f;

}

void f_func_gaussmf(FIXED_POINT_U_PRECISION x, FIXED_POINT_PRECISION mean, FIXED_POINT_PRECISION* result)
{	
	//gaussian member function e ^-((x-mean)^2/(2*stv^2)
	
	FIXED_POINT_PRECISION mean_term;//7/2*pi most likey, so max around 25 
	FIXED_POINT_U_PRECISION mean_term_sq;// max around 25*25
	FIXED_POINT_EXP exp_pwr ;
	//FIXED_POINT_PRECISION stv_sq_2;
	FIXED_POINT_EXP divisor_flip;
	FIXED_POINT_PRECISION constant_two=2;
	FIXED_POINT_PRECISION constant_one=1;
	
	mean_term=x-mean;
	if(hls::abs(mean_term)>0.9)//anything larger than this saturate to 0 basically
	{
		mean_term_sq=0.9;
	}
	else
	{
		mean_term_sq=mean_term*mean_term;//MAX 1

	}
	//stv_sq=stv*stv;
	//stv_sq_2=constant_two*stv_sq;//USING FIXED STV_SQ
	//stv_sq_2=0.02;//

	//divisor_flip=constant_one/stv_sq2;//50
	divisor_flip=1250;
	
	exp_pwr=-mean_term_sq*divisor_flip; //MAX50
	
	
	*result=hls::exp(exp_pwr);

}

*/
void f_func_lut(FIXED_POINT_U_PRECISION theta, FIXED_POINT_U_PRECISION* result,int k_part)

{
#pragma HLS PIPELINE

	//printf ("theta %f \n", theta. to_double());
	const FIXED_POINT_U_PRECISION lut_theta[201]={0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,1,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.09,1.1,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.2,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.3,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.4,1.41,1.42,1.43,1.44,1.45,1.46,1.47,1.48,1.49,1.5,1.51,1.52,1.53,1.54,1.55,1.56,1.57,1.58,1.59,1.6,1.61,1.62,1.63,1.64,1.65,1.66,1.67,1.68,1.69,1.7,1.71,1.72,1.73,1.74,1.75,1.76,1.77,1.78,1.79,1.8,1.81,1.82,1.83,1.84,1.85,1.86,1.87,1.88,1.89,1.9,1.91,1.92,1.93,1.94,1.95,1.96,1.97,1.98,1.99,2};
	const FIXED_POINT_U_PRECISION lut_k4[201]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0005269435,0.003436103,0.01744997,0.06901597,0.2125842,0.5099629,0.9527361,1.386223,1.570796,1.386223,0.9527361,0.5099629,0.2125842,0.06901597,0.01744997,0.003436103,0.0005269435,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0002517377,0.002107774,0.01374441,0.06979989,0.2760639,0.8503367,2.039852,3.810945,5.544891,6.283185,5.544891,3.810945,2.039852,0.8503367,0.2760639,0.06979989,0.01374441,0.002107774,0.0002517377,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000440541,0.003688604,0.02405272,0.1221498,0.4831118,1.488089,3.56974,6.669153,9.70356,10.99557,9.70356,6.669153,3.56974,1.488089,0.4831118,0.1221498,0.02405272,0.003688604,0.000440541,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	const FIXED_POINT_U_PRECISION lut_k3[201]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000674793,0.004048382,0.01891554,0.06883075,0.195062,0.4305158,0.7400004,0.9906064,1.032754,0.8385302,0.5302334,0.261121,0.1001482,0.02991372,0.00695864,0.001260679,0.0001778735,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0008893675,0.006303396,0.0347932,0.1495686,0.5007408,1.305605,2.651167,4.192651,5.163768,4.953032,3.700002,2.152579,0.9753098,0.3441538,0.09457771,0.02024191,0.003373965,0.0004379817,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	const FIXED_POINT_U_PRECISION lut_k2[201]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	ap_fixed<25,13> entrance_lower_f;
	ap_int<11> entrance_lower;
	int entrance_lower_int;
	FIXED_POINT_U_PRECISION lut_out_lower;
	FIXED_POINT_U_PRECISION lut_out_upper;
	FIXED_POINT_U_PRECISION lut_out;
	FIXED_POINT_U_PRECISION theta_in;

	FIXED_POINT_U_PRECISION lut_out_diff;
	FIXED_POINT_U_PRECISION distance_to_lower;
	FIXED_POINT_U_PRECISION theta_lower;
	FIXED_POINT_U_PRECISION linear_approx_amount;

	ap_fixed <11,10>   f_func_step_inverse=100;
	FIXED_POINT_U_PRECISION f_func_step=0.01;

	theta_in=hls::abs(theta);
	entrance_lower_f=theta_in*f_func_step_inverse;
	//entrance_lower=hls::floor(entrance_f);
	//printf("entrance f %f \n",entrance_f.to_double());
	//entrance_lower=hls::floor(entrance_f);
	entrance_lower=entrance_lower_f.to_ap_int();

	entrance_lower_int=entrance_lower.to_int();

	//printf("entrance lower int %d \n",entrance_lower_int);

	if (entrance_lower_int>199)
	{
		lut_out_lower=0;
		lut_out_upper=0;
		lut_out_diff=0;
		distance_to_lower=0;
		//printf(" >200 \n");
	}
	else
	{
		switch (k_part)
		{
		case 2:
			lut_out_lower=lut_k2[entrance_lower_int];
			lut_out_upper=lut_k2[entrance_lower_int+1];
			break;
		case 3:
			lut_out_lower=lut_k3[entrance_lower_int];
			lut_out_upper=lut_k3[entrance_lower_int+1];
			break;
		case 4:
			lut_out_lower=lut_k4[entrance_lower_int];
			lut_out_upper=lut_k4[entrance_lower_int+1];
			break;
		default:
			lut_out_lower=lut_k2[entrance_lower_int];
			lut_out_upper=lut_k2[entrance_lower_int+1];

			break;
		}
	
	//printf("lut out lower %f \n",lut_out_lower.to_double());


	lut_out_diff=lut_out_upper-lut_out_lower;
	theta_lower=lut_theta[entrance_lower_int];
	distance_to_lower=theta_in-theta_lower;

	}

	if (distance_to_lower>(FIXED_POINT_U_PRECISION)0.0001&&distance_to_lower<(FIXED_POINT_U_PRECISION)0.01)//make sure range is correct
	{
		if(lut_out_diff!=0)//make sure range is correct
		{	
			//printf ("linear approx \n");
			ap_fixed<25,13> linear_approx_ratio_expand;
			FIXED_POINT_U_PRECISION linear_approx_ratio;

			linear_approx_ratio_expand=distance_to_lower*f_func_step_inverse;//use result to controll fix bit allignment
			//printf("ratio approx %f \n",linear_approx_ratio.to_double());	
			linear_approx_ratio=linear_approx_ratio_expand;
			linear_approx_amount=lut_out_diff*linear_approx_ratio;
			//printf("amount approx %f \n",linear_approx_amount.to_double());

			
				lut_out=lut_out_lower+linear_approx_amount;
			
			
		}
		else
		{
		lut_out=lut_out_lower;
		}
		
	}
	else
	{

		lut_out=lut_out_lower;

	}

	
	

	if (theta<0)
	{
		lut_out=-lut_out;
	}
	//printf(" lut lowere %lf \n",lut_out_lower.to_double());
	//printf(" lut upper  %lf \n",lut_out_upper.to_double());
	//printf("lut %f \n",lut_out.to_double());

	*result=lut_out;
	


}



















