#include "mex.h"
#include "matrix.h"
#include <memory.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>



#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define lp 696


#define	N_IN1	prhs[0]
//#define	N_IN2	prhs[1]


#define	InterleavingArray_OUT	plhs[0]


void interleav_matrix(double * addN, double intInterleavingArray[])//double rand_list[]//interleav_matrix//double * addN, double intInterleavingArray[]
{
	
	//int tt=9984;
	//int *addN=&tt ;
	int intInterArrayLen = 0;
	int subscript;
	//int xx = 0;
	//int *intInterleavingArray = (int *)malloc(sizeof(int)*(*addN));//9984
	 int pilot_num_used,subcarriers_in_use,Lp;
	
    
	int intNumberArray[3*(lp+4)] ;//************************************89984//= (int *)malloc(sizeof(int)*(*addN))=2(lp+8)
	int intRefuseArray;		//拒绝数组
	int intRandNum;			//随机数
	int intNumArrayLen;		//数字数组长度
	int Flag;			//标志位
	int I = 0;
	int N ;
	int s = 0;

	
	
	double U;	//true.
	double x;	//true.
	static int i;
	static __int64 J,j;
	static double b ;
	static double z[32] = {-1};
	static double ulp;
	__int64 mask;
	double f;
	int e;	//true.
	int k,d;

	
	//FILE *fp,*fp_2;

	//srand( (unsigned)time( NULL ) );
	//srand(11);
	N = *addN;

	//将intInterleavingArray置空；
	for ( I = 0 ; I < N ; I++ )
	{
		intInterleavingArray[I] = -1;        
	}
	//对int_interleaving_array进行判断，如果长度 <= N,进行int_interleaving_array的设计
	while ( intInterArrayLen != N )
	{
		//对各个数组重新进行初始化,将intInterleavingArray置空;
		for ( I = 0 ; I < N ; I++ )
		{
			intInterleavingArray[I] = -1;
		}
		intInterArrayLen = 0;
		//生成intInterleavingArray，长度为N，存储1--N这N个整数；
		for ( I = 0 ; I < N ; I++ )
		{
			intNumberArray[I] = I+1;
		}
		//生成intRefuseArray，内容为空；
		intRefuseArray = -1;
		intNumArrayLen = N;

		//对数字数组的长度进行判断，不为空时，继续对交织器进行设计；
		while  ( intNumArrayLen > 0 )
		{
			//subscript = (int)( intNumArrayLen * rand_list[rand_index]+1);
			//rand_index++;
			

			//产生0--1的随机数；	
			if ( z[0] == -1 )
			{
				j = (__int64)pow( 2 , 31 );
				for ( k = 0; k < 32; k++)
				{
					x = 0;
					for ( d = 0; d < 53 ; d++)
					{
						j = j ^ (j << 13)  & 4294967295;
						j = j ^ (j >> 17)  & 4294967295;
						j = j ^ (j << 5 )  & 4294967295;//32bit
						x = x * 2 + (double)((j >> 19)  & 1 ); //53bit
					}
					z[k] =x * pow(2,-53);
				}
				j = (__int64)pow( 2 , 31 );
				b = 0;
				i = 0;
				ulp = pow( 2,-53);
			}


			U = 0;
			x = z[(i+20)%32] - z[(i+5)%32] - b;
			if ( x < 0 )
			{
				x = x + 1;
				b = ulp;
			}
			else
			{
				b = 0;
			}
			z[i] = x;
			i = i + 1;
			if (i == 32)
			{
				i = 0;
			}
			J = j;
			j = j ^ (j << 13 )  & 4294967295;
			j = j ^ (j  >> 17 )  & 4294967295;
			j = j ^ (j  << 5  )  & 4294967295;//32bit

			mask = ((j << 32)  & 4503599627370495) ^ J;//52bit
			f = frexp( x, &e );//[f,e] = log2(x);
			x = ((__int64)(f * pow( 2 , 53 )) ^ mask)  * pow ( 2 , e-53);//x = pow2(bitxor(pow2(f,53),mask),e-53);
			U = x;				//U既是随机数，此算法采用和matlab随机数产生相同的算法


			//从1至N这N个数字中随机获取一个下标subscipt；
			subscript = (int)(intNumArrayLen * U + 1);
			/*fp_2 = fopen("RandNum.txt","a");
			if (fp_2 != NULL)
			{
				fprintf(fp_2,"%d %d \n",intNumArrayLen,subscript);
			}
			fclose(fp_2);*/
			//从这个下标所表示的位置提取所选的随机数intRandNum；
			intRandNum = intNumberArray[subscript-1];
			//将intNumberArray中最后一个值存入该数组的subscript位置；
			intNumberArray[subscript-1] = intNumberArray[intNumArrayLen-1];
			//长度减1；
			intNumArrayLen--;
	
			//测试其S-距离特性,设定判断距离S = [（N/2）^0.5]-1
			s = (int)pow( (N/2),0.5 ) - 1;
			//初始化flag的值;
			Flag = 1;
			//如果intInterArrayLen的长度小于s，将s设为intInterArrayLen的长度;
			if  (intInterArrayLen < s)
			{
				s = intInterArrayLen;
			}
			//如果int_number_array的长度不大于2*s，s的值减半;
			if ( (intNumArrayLen+1) <= (2 * s))
			{
				s = (int)( (intNumArrayLen+1)/2 );
			}
			//将intRandNum与intInterleavingArray中的最后S个值进行比较
			for ( I = 1 ; I <= s ;I++ )				
			{
				//对标志位flag进行设定；
				if ( abs( intInterleavingArray [intInterArrayLen - i] - intRandNum ) <= s )	
				{
					// 如果距离小于s,则flag置0，跳出循环;
					Flag = 0;
					break;
				}
			}
			if  (Flag == 0)
			{
				intRefuseArray = intRandNum;
			}
			else
			{
				//将此随机数记入交织数组；
				intInterleavingArray[intInterArrayLen] = intRandNum;	
				//长度加1；
				intInterArrayLen++;
			}
			/*for (i = 0 ; intInterleavingArray[i] != -1; )
			{
				i++;
			}
			intInterArrayLen = i;*/
			//对拒绝数组进行处理
			if  ( intRefuseArray != -1 )
			{
				//将拒绝数组中的数字记入intNumberArray;
				intNumberArray[intNumArrayLen] = intRefuseArray;
				//清空拒绝数组;
				intRefuseArray = -1;
				//数字数组加1；
				intNumArrayLen++;
			}
		
		}
		/*fp = fopen("interleavfirst.txt","w");
		if (fp != NULL)
		{
			for (xx=0;xx<N;xx++)
			{
				fprintf(fp,"%d ",intInterleavingArray[xx]);
			}
		}
		fclose(fp);*/
	}




	/*fp = fopen("interleavfirst.txt","w");
	if (fp != NULL)
	{
		for (x=0;x<(* addN);x++)
		{
			fprintf(fp,"%d ",intInterleavingArray[x]);
		}
	}
	
	fclose(fp);*/
	//free (intInterleavingArray);
	
	//free (intNumberArray);
	
	_CrtDumpMemoryLeaks();
	
	return;
	
}


//转换成mex文件的接口函数
//mexFunction
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{ 
	//double * rand_list;
	double * intInterleavingArray; 
	double * addN; 
    
    /*Check for proper number of arguments */
    
    if (nrhs != 1) 
	{
		mexErrMsgTxt("one input arguments required."); 
    }
	else if (nlhs > 1) 
	{
		mexErrMsgTxt("Too many output arguments."); 
	} 
    
    
	addN = mxGetPr(N_IN1); 
	//rand_list = mxGetPr(N_IN2);
    /* Create a matrix for the return argument */
	InterleavingArray_OUT = mxCreateDoubleMatrix(1, * addN, mxREAL); 
    
    /* Assign pointers to the various parameters */
	intInterleavingArray = mxGetPr(InterleavingArray_OUT);
    
	
    /* Do the actual computations in a subroutine */
	interleav_matrix(addN,intInterleavingArray); //rand_list
    
}







