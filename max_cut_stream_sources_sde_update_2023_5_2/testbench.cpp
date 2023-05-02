#include "f_func.h"

#include <cstring>
#include <cstdlib>
#include <random>
int j_sparse_read(char* J_table_file_name,int*J_table,int* j_length_out);

int main()
{


	char J_sparse_file_name1[100];
	char J_sparse_file_name2[100];
	char J_sparse_file_name3[100]="/home/centos/Ascent/n30_4_.6_edges.csv";
	char J_sparse_file_name800[100]="/home/centos/Ascent/G1_edges.csv";
	char J_sparse_file_name2000[100]="/home/centos/Ascent/G22_edges.csv";



	//char line_buffer[N_size*2];
	FIXED_POINT result;

	double theta_double_all[N_size]={0,-0.1821,1.3755,-0.1973,-0.3647,1.4344,1.0725,0.5357,-0.2756,-0.3003, 1.6294,1.8116,0.2540,1.8268,1.2647,0.1951,0.5570,1.0938,1.9150, 1.9298,1.4344,1.0725,0.5357,-0.2756,-0.3003,1.6294,1.8116,0.2540,1.8268,1.5};
	int* J_table_int=(int*) malloc(2*N_size * N_size*sizeof(int));

	int n_size_n=2000;
	int j_length;


	
	double C_1=5;
	double C_sync=0.5;

	
	double J_d_out[100];

    double tstop=0.005;
    double tstep=0.005;

    double norm_rand_mean=0;
    double norm_rand_deviation=1.0;


    int iteration_size=ceil(tstop/tstep);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(norm_rand_mean,norm_rand_deviation);

	switch (n_size_n)
	{
	case 10 :


		if(j_sparse_read(J_sparse_file_name1,J_table_int,&j_length)==0)
		{
			return 0;
		}
		/* code */
		break;
	case 20 :

		if(j_sparse_read(J_sparse_file_name2,J_table_int,&j_length)==0)
		{
			return 0;
		}

		break;

	case 30:
		if(j_sparse_read(J_sparse_file_name3,J_table_int,&j_length)==0)
		{
			return 0;
		}

		break;

	case 800:

		for (int i = 0; i < N_size; i++)
		{
			theta_double_all[i]=distribution(generator);
		}
		
		if(j_sparse_read(J_sparse_file_name800,J_table_int,&j_length)==0)
		{
			return 0;
		}
		break;
	case 2000:

		for (int i = 0; i < N_size; i++)
		{
			theta_double_all[i]=distribution(generator);
		}
		
		if(j_sparse_read(J_sparse_file_name2000,J_table_int,&j_length)==0)
		{
			return 0;
		}
		break;

	default:
		return 0;
		break;
	}

/*


	for (int i = 0; i < j_length*2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			printf("%d ", J_table_int[i*2+j]);
		}

		printf("\n");
		
	}
	
*/



	double* phi_d_out = (double*) malloc(8000 * N_size*sizeof(double));
	//f_func_eq(lower_k,upper_K,theta,sigma,result);

	//result_double=result.to_double();

	//printf(" f kernel result \n");



	//f_func_accumulate(theta,sigma,&result;



	//result_double=result.to_double();


	//printf(" f  accumulate result %lf end  \n ",result_double);
  
  	double* dwt = (double*) malloc(8000 * N_size*2*sizeof(double));

  double An=0.1;
  for(int i=0;i<iteration_size*N_size*2;i++)
    {
    
            dwt[i]=An*sqrt(0.5*tstep)*distribution(generator);
			//printf("%f ",dwt[i]);

    }

int k_part=2;


//kuramoto_dataflow_top(theta_double_all,J_table_int,C_1,C_sync,tstep,iteration_size,k_part,dwt,phi_d_out,n_size_n,j_length,1);

FIXED_POINT_LOW Ac_fix=20;
FIXED_POINT_PRECISION C_sync_fix=5;
FIXED_POINT_PRECISION Kt_fix=10;
FIXED_POINT_PRECISION T_fix=1;
FIXED_POINT_PRECISION t_stop_fix=10;
FIXED_POINT_U_PRECISION t_step_fix=0.005;

FIXED_POINT_LOW C1_result;
FIXED_POINT_PRECISION C_sync_result;

FIXED_POINT_PRECISION t_stop_inverse_fix=(FIXED_POINT_PRECISION)1/t_stop_fix;
FIXED_POINT_PRECISION T_inverse_fix=(FIXED_POINT_PRECISION)1/T_fix;

top_annealing_update(Ac_fix,C_sync_fix,Kt_fix,T_inverse_fix,t_step_fix,t_stop_inverse_fix,&C1_result,&C_sync_result,5);

printf(" C_1 %f C_s %f \n",C1_result.to_double(),C_sync_result.to_double());

		 	
/*
	for (int i=0;i<200;i++)
	{	
		
		f_func_accumulate(phi,&f_result);
		printf("%f,",f_result.to_double());

		if (f_result>(FIXED_POINT_PRECISION)0.0001)
		{
			printf("  start index %d ",i);
		}
		phi=phi+(FIXED_POINT_PRECISION)0.01;
	}
*/



FILE* output_phi;



    char output_file_name[100]="/home/centos/Ascent/HLS/max_cut_stream_sources_sde_update/BLOCK_K2_out.txt";
    output_phi=fopen(output_file_name,"w");


    for (int i=0;i<iteration_size;i++)
    {
        for(int j=0;j<N_size;j++)
        {
            fprintf(output_phi,"%f  ",phi_d_out[i*N_size+j]);
			//printf("%f  ",phi_d_out[i*N_size+j]);

        }
        fprintf(output_phi," \n");
		//printf("\n");
    }

    fclose(output_phi);



	return 0;

}
/*
int j_map_read(char* J_table_file_name,char*J_table)
{
char c='s';
	int n=N_size;

	FILE * J_table_stream;

	
	fopen(J_table_file_name,"r");
	
	
	int m=0;

	if(J_table_stream==NULL)
	{
		printf("error reading file \n");
		return 0;
	}
	else
	{
		while(c!=EOF)

		{
			c=fgetc(J_table_stream);
			if (c>'/'&&c<'2')
			{
				J_table[m]=c-'0';
				m++;


			}
		}
	}

	printf("J readed %d\n",m);



	for (int i=0;i<N_size;i++)
		{
			for(int j=0;j<N_size;j++)
			{
				//printf(" %f ",J_table_double[i*N_size+j]);
			}
			//printf("\n");
		}

}
*/

int j_sparse_read(char* J_table_file_name,int*J_table,int* j_length_out)
{

	FILE* read_file_stream=fopen(J_table_file_name,"r");
	int read_length=100;
	char read[100];
	int n_size;
	int j_length;
	if(read_file_stream==NULL)
{
    printf("error reading fopen file \n");
    return 0;
}
else{
    fgets(read, read_length, read_file_stream);
  //  printf("%s\n",read);
}






char *token=strtok(read,",");

n_size=std::atoi(token);

token = strtok(NULL,",");

j_length=std::stoi(token);

token = strtok(NULL,",");


char buffer[100];
for (int i = 0; i < j_length*2; i++)
{
	fgets(buffer,read_length,read_file_stream);
	
	char* token=strtok(buffer,",");
	J_table[i*2]=std::atoi(token);

	token=strtok(NULL,",");
	J_table[i*2+1]=std::atoi(token);


	token=strtok(NULL,",");


	}


*j_length_out=j_length;

return 1;

}
