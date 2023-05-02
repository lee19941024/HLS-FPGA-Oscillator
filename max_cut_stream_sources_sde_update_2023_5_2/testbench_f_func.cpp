#include "f_func.h"



int main()
{
    
	FIXED_POINT_U_PRECISION phi;
	FIXED_POINT_PRECISION f_result;
	FIXED_POINT_PRECISION f_result_lut;
	phi=-2.5;



	for (int i=0;i<1000;i++)
	{	
		
		//f_func_accumulate(phi,&f_result);
		f_func_lut(phi,&f_result_lut);
		f_func_accumulate(phi,&f_result);
		printf("iteration %d %f  %f \n ",i,f_result.to_double(),f_result_lut.to_double());

		phi=phi+(FIXED_POINT_U_PRECISION)0.005;
	}

}