#include "f_func.h"



int main()
{


    FIXED_POINT_U_PRECISION sin1;
	FIXED_POINT_U_PRECISION sin2;
	FIXED_POINT_U_PRECISION sin3;
	FIXED_POINT_U_PRECISION sin4;
	FIXED_POINT_PRECISION sin5;


	FIXED_POINT_PRECISION degree=0.785398;


	for	(int j=0;j<2048;j++)
	{

		cordic_sin(degree,&sin1);
        sin5=hls::sin(degree);

		printf(" reference %f hybrid %f \n", sin5.to_double(),sin1.to_double());
		degree=degree+(FIXED_POINT_PRECISION)0.000767;
	}




	cordic_sin(-degree-6*pi_r,&sin1);
	cordic_sin(-pi_q2+degree-6*pi_r,&sin2);
	cordic_sin(-pi_q2-degree-6*pi_r,&sin3);
	cordic_sin(-pi_r+degree-6*pi_r,&sin4);
	sin5=hls::sin(degree);

	printf(" cordic sin deg %f \n", sin1.to_double());
	printf(" cordic sin 180-deg %f \n",sin2.to_double());
	printf(" cordic sin deg+180 %f \n",sin3.to_double());
	printf(" cordic sin 360-deg %f \n",sin4.to_double());
	printf(" cordic sin reference %f \n",sin5.to_double());

}
