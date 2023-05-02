#include "f_func.h"
void dwt_generator(FIXED_POINT_RANDOM seed[2*DATAFLOW_PARALLEL_FINAL_SUM], FIXED_POINT_U_PRECISION dwt[2*DATAFLOW_PARALLEL_FINAL_SUM], bool init, FIXED_POINT_U_PRECISION An, FIXED_POINT_U_PRECISION tstep)
{
    static FIXED_POINT_RANDOM lfsr[2*DATAFLOW_PARALLEL_FINAL_SUM];
    
    #pragma HLS ARRAY_PARTITION dim=0 factor=DATAFLOW_PARALLEL_FINAL_SUM*2 type=cyclic variable=lfsr
    #pragma HLS ARRAY_PARTITION dim=0 factor=DATAFLOW_PARALLEL_FINAL_SUM*2 type=cyclic variable=seed
    #pragma HLS ARRAY_PARTITION dim=0 factor=DATAFLOW_PARALLEL_FINAL_SUM*2 type=cyclic variable=dwt


    for (int i = 0; i <DATAFLOW_PARALLEL_FINAL_SUM*2 ; i++)
    {
        #pragma HLS UNROLL
        if(init)
        {
            lfsr[i]=seed[i];
        }
        bool b_18=lfsr[i][18];
        bool b_13=lfsr[i][13];
        bool b_2=lfsr[i][2];
        bool b_1=lfsr[i][1];
        bool new_bit=b_18^b_13^b_2^b_1;

        lfsr[i]=lfsr[i]>>1;
        lfsr[i][20]=new_bit;
        dwt[i]=An*sqrt((FIXED_POINT_U_PRECISION)0.5*tstep);
        
    }
    
    


}
