#include "mex.h"
#include "matrix.h"

#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>

#define	IN_PARAMETER_1	prhs[0]
#define	IN_PARAMETER_2	prhs[1]

#define	OUT_PARAMETER_1	plhs[0]

#define Lp 696
#define DECODE_SOURCE_LENGTH  3*(Lp+4)
#define DECODE_LENGTH (DECODE_SOURCE_LENGTH-12)/3

double max_2(double in1, double in2)
{
	if(in1>in2)
		return in1;
	else
		return in2;
}

double max_8(double *in)
{
	double max_temp;
	int max_i;
	max_temp=in[0];
	for(max_i=1;max_i<8;max_i++)
	{
		if (in[max_i]>max_temp)
			max_temp=in[max_i];	
	}
	return max_temp;
}


void sub_decoder(double *z, double *x, double *y, int dec_length, double *result, double *L)
{
	int i,gama_i,beta_i,alfa_i;	
	
	double branch1[DECODE_LENGTH+3];
	double branch2[DECODE_LENGTH+3];
	
	double gama_00[DECODE_LENGTH+3];
	double gama_04[DECODE_LENGTH+3];
	double gama_14[DECODE_LENGTH+3];
	double gama_10[DECODE_LENGTH+3];
	double gama_25[DECODE_LENGTH+3];
	double gama_21[DECODE_LENGTH+3];
	double gama_31[DECODE_LENGTH+3];
	double gama_35[DECODE_LENGTH+3];
	double gama_42[DECODE_LENGTH+3];
	double gama_46[DECODE_LENGTH+3];
	double gama_56[DECODE_LENGTH+3];
	double gama_52[DECODE_LENGTH+3];
	double gama_67[DECODE_LENGTH+3];
	double gama_63[DECODE_LENGTH+3];
	double gama_73[DECODE_LENGTH+3];
	double gama_77[DECODE_LENGTH+3];
	
	double beta1[DECODE_LENGTH+4];
	double beta2[DECODE_LENGTH+4];
	double beta3[DECODE_LENGTH+4];
	double beta4[DECODE_LENGTH+4];
	double beta5[DECODE_LENGTH+4];
	double beta6[DECODE_LENGTH+4];
	double beta7[DECODE_LENGTH+4];
	double beta8[DECODE_LENGTH+4];
	
	double alfa[8]={0,-100000,-100000,-100000,-100000,-100000,-100000,-100000};
	
	
	double s[8]={0,0,0,0,0,0,0,0};
	double s0,s1;
	
	double out_00,out_01,out_10,out_11;
	
	for(i=0;i<dec_length;i++)
	{
		branch1[i]=(x[i]+z[i])/2;
		branch2[i]=y[i]/2;		
	}
	branch1[dec_length]=x[dec_length]/2;
	branch1[dec_length+1]=x[dec_length+1]/2;
	branch1[dec_length+2]=x[dec_length+2]/2;
	branch2[dec_length]=y[dec_length]/2;
	branch2[dec_length+1]=y[dec_length+1]/2;
	branch2[dec_length+2]=y[dec_length+2]/2;
	
	for(gama_i=0;gama_i<(dec_length+3);gama_i++)
	{
		out_00=-branch1[gama_i]-branch2[gama_i];
		out_01=-branch1[gama_i]+branch2[gama_i];
		out_10=branch1[gama_i]-branch2[gama_i];
   		out_11=branch1[gama_i]+branch2[gama_i];
   		
   		gama_00[gama_i]=out_00; 	
		gama_04[gama_i]=out_11;
		gama_14[gama_i]=out_00;
		gama_10[gama_i]=out_11;
		gama_25[gama_i]=out_01;
		gama_21[gama_i]=out_10;
		gama_31[gama_i]=out_01;
		gama_35[gama_i]=out_10;
		gama_42[gama_i]=out_01;
		gama_46[gama_i]=out_10;
		gama_56[gama_i]=out_01;
		gama_52[gama_i]=out_10;
		gama_67[gama_i]=out_00;
		gama_63[gama_i]=out_11;
		gama_73[gama_i]=out_00;
		gama_77[gama_i]=out_11;		
	}
	
	for(i=0;i<(dec_length+4);i++)
	{
		if(i==(dec_length+3))
		{
			beta1[i]=0;	
	        	beta2[i]=-100000;
                	beta3[i]=-100000;
                	beta4[i]=-100000;
                	beta5[i]=-100000;
                	beta6[i]=-100000;
                	beta7[i]=-100000;
                	beta8[i]=-100000;		
		}
		else
		{
			beta1[i]=0;	
	        	beta2[i]=0;
                	beta3[i]=0;
                	beta4[i]=0;
                	beta5[i]=0;
                	beta6[i]=0;
                	beta7[i]=0;
                	beta8[i]=0;	
		}      
        }

	for(beta_i=(dec_length+2);beta_i>-1;beta_i--)
        {
        	beta1[beta_i]=max_2(gama_00[beta_i]+beta1[1+beta_i],gama_04[beta_i]+beta5[1+beta_i]);
        	beta2[beta_i]=max_2(gama_10[beta_i]+beta1[1+beta_i],gama_14[beta_i]+beta5[1+beta_i]);
      		beta3[beta_i]=max_2(gama_21[beta_i]+beta2[1+beta_i],gama_25[beta_i]+beta6[1+beta_i]);
      		beta4[beta_i]=max_2(gama_31[beta_i]+beta2[1+beta_i],gama_35[beta_i]+beta6[1+beta_i]);
      		beta5[beta_i]=max_2(gama_42[beta_i]+beta3[1+beta_i],gama_46[beta_i]+beta7[1+beta_i]);
      		beta6[beta_i]=max_2(gama_52[beta_i]+beta3[1+beta_i],gama_56[beta_i]+beta7[1+beta_i]);
      		beta7[beta_i]=max_2(gama_63[beta_i]+beta4[1+beta_i],gama_67[beta_i]+beta8[1+beta_i]);
      		beta8[beta_i]=max_2(gama_73[beta_i]+beta4[1+beta_i],gama_77[beta_i]+beta8[1+beta_i]);
       }
        
        for(alfa_i=0;alfa_i<dec_length;alfa_i++)
        {
        	s[0]=alfa[0]+gama_00[alfa_i]+beta1[1+alfa_i];
        	s[1]=alfa[1]+gama_14[alfa_i]+beta5[1+alfa_i];
   		s[2]=alfa[2]+gama_25[alfa_i]+beta6[1+alfa_i];
   		s[3]=alfa[3]+gama_31[alfa_i]+beta2[1+alfa_i];
   		s[4]=alfa[4]+gama_42[alfa_i]+beta3[1+alfa_i];
   		s[5]=alfa[5]+gama_56[alfa_i]+beta7[1+alfa_i];
   		s[6]=alfa[6]+gama_67[alfa_i]+beta8[1+alfa_i];
   		s[7]=alfa[7]+gama_73[alfa_i]+beta4[1+alfa_i];
   		s0=max_8(s);
   		
   		s[0]=alfa[0]+gama_04[alfa_i]+beta5[1+alfa_i];
        	s[1]=alfa[1]+gama_10[alfa_i]+beta1[1+alfa_i];
   		s[2]=alfa[2]+gama_21[alfa_i]+beta2[1+alfa_i];
   		s[3]=alfa[3]+gama_35[alfa_i]+beta6[1+alfa_i];
   		s[4]=alfa[4]+gama_46[alfa_i]+beta7[1+alfa_i];
   		s[5]=alfa[5]+gama_52[alfa_i]+beta3[1+alfa_i];
   		s[6]=alfa[6]+gama_63[alfa_i]+beta4[1+alfa_i];
   		s[7]=alfa[7]+gama_77[alfa_i]+beta8[1+alfa_i];
   		s1=max_8(s);
   		
   		result[alfa_i]=s1-s0;	
   		
   		s[0]=max_2(alfa[0]+gama_00[alfa_i],alfa[1]+gama_10[alfa_i]);
		s[1]=max_2(alfa[3]+gama_31[alfa_i],alfa[2]+gama_21[alfa_i]);
		s[2]=max_2(alfa[4]+gama_42[alfa_i],alfa[5]+gama_52[alfa_i]);
		s[3]=max_2(alfa[7]+gama_73[alfa_i],alfa[6]+gama_63[alfa_i]);
		s[4]=max_2(alfa[1]+gama_14[alfa_i],alfa[0]+gama_04[alfa_i]);
		s[5]=max_2(alfa[2]+gama_25[alfa_i],alfa[3]+gama_35[alfa_i]);
		s[6]=max_2(alfa[5]+gama_56[alfa_i],alfa[4]+gama_46[alfa_i]);
		s[7]=max_2(alfa[6]+gama_67[alfa_i],alfa[7]+gama_77[alfa_i]);
		
		alfa[0]=s[0];
		alfa[1]=s[1];
		alfa[2]=s[2];
		alfa[3]=s[3];
		alfa[4]=s[4];
                alfa[5]=s[5];
                alfa[6]=s[6];
                alfa[7]=s[7];
        }
        for(i=0;i<dec_length;i++)
        {
        	L[i]=result[i]-x[i]-z[i];
        	if(result[i]>0)
        		result[i]=1;
        	else
        		result[i]=0;	
        }
}       
        
void turbo_decode(double decode_source[], double interleave_order[], double out_put[])
{       
	int dec_length=(DECODE_SOURCE_LENGTH-12)/3;
	int branch_i, decode_time,z_i,i,j;
	
	double result[DECODE_LENGTH]; 
	double L[DECODE_LENGTH];
	
	double branch1[DECODE_LENGTH+3];      
	double branch2[DECODE_LENGTH+3];
	double branch3[DECODE_LENGTH+3];
	double branch4[DECODE_LENGTH+3];
	
	double x[DECODE_LENGTH+3];
	double y[DECODE_LENGTH+3];
	double z[DECODE_LENGTH+3];
	
	
	
	for(branch_i=0;branch_i<dec_length;branch_i++)
	{
		branch1[branch_i]=-decode_source[branch_i*3];
	        branch2[branch_i]=-decode_source[branch_i*3+1];
               	branch4[branch_i]=-decode_source[branch_i*3+2]; 
	}
	
	for(branch_i=0;branch_i<dec_length;branch_i++)
	{
		i=interleave_order[branch_i];
		branch3[branch_i]=branch1[i]; 
	}
	
	for(branch_i=0;branch_i<3;branch_i++)
	{
		branch1[dec_length+branch_i]=-decode_source[3*dec_length+2*branch_i];
   		branch2[dec_length+branch_i]=-decode_source[3*dec_length+2*branch_i+1];
   		branch3[dec_length+branch_i]=-decode_source[3*dec_length+2*branch_i+6];
   		branch4[dec_length+branch_i]=-decode_source[3*dec_length+2*branch_i+7];	
	}
	
	for(i=0;i<dec_length;i++)
	{
		result[i]=0;
		L[i]=0;
		z[i]=0;	
	}
	
	for(decode_time=0;decode_time<5;decode_time++)
	{
		for(z_i=0;z_i<dec_length;z_i++)
   		{
   			i=interleave_order[z_i];
   			z[i]=L[z_i];
   		}
   		for(i=0;i<(dec_length+3);i++)
   		{
   			x[i]= branch1[i];
   			y[i]= branch2[i];
   		}
   		
   		sub_decoder(z,x,y,dec_length,result,L);
   		
   		for(z_i=0;z_i<dec_length;z_i++)
   		{	
   			i=interleave_order[z_i];
   			z[z_i]=L[i];
   		}
   		for(i=0;i<(dec_length+3);i++)
   		{
   			x[i]= branch3[i];
   			y[i]= branch4[i];
   		}
   		sub_decoder(z,x,y,dec_length,result,L);
	}
	for(i=0;i<dec_length;i++)
	{
		j=interleave_order[i];
		out_put[j]=result[i];	
	}
}       
        
//mexFunction
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{       
  double *out_put,*out_put2,*out_put3; 
  double *decode_source;
  double *interleave_order;
  int i;
        
    /*Check for proper number of arguments */
    
    if (nrhs != 2) 
	{
		mexErrMsgTxt("two input arguments required."); 
    }
	else if (nlhs > 1) 
	{
		mexErrMsgTxt("Too many output arguments."); 
	} 
    
    
    /* Create a matrix for the return argument */
   OUT_PARAMETER_1 = mxCreateDoubleMatrix(1, ((DECODE_SOURCE_LENGTH-12)/3), mxREAL); 

    
    /* Assign pointers to the various parameters */
   out_put = mxGetPr(OUT_PARAMETER_1);
    
   decode_source = mxGetPr(IN_PARAMETER_1);
   interleave_order = mxGetPr(IN_PARAMETER_2);    
        
    /* Do the actual computations in a subroutine */
   turbo_decode(decode_source,interleave_order,out_put);
    
}

  