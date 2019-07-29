/******************************************/
/***        Author: Pei-Wei Tsai        ***/
/***   E-mail: pwtsai@bit.kuas.edu.tw   ***/
/***         Version: 2.0.0.0           ***/
/***   Released on: November 1, 2009    ***/
/******************************************/

/* The "Benchmark.h" and the "Benchmark.cpp" composes a library contains i_FuncNum test functions. */
#include "Benchmark.h"			// Include its header file.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;

double **p;
int NumDimofP;
int NumDataofP;
double* T;
int numNeuron;
int epoch;

double MF(double point , double center , double spread){
  double b = 1;
  double r3;
  //cout<<"speard:"<<(double)(point - center) / ( spread + 1)<<endl;
  r3 = exp(-1 * ( pow ( (double)(point - center) / ( spread +  0.1) , 2 * b )) );
  return r3;
}


double g1(int x){
  double r5;
  if(x > 0 && x < 250 )
    r5 = sin( ( x * M_PI ) / 25 );
  else if(x >= 250 && x <= 500)
    r5 = 1;
  else if(x >= 500 && x < 750)
    r5 = -1;
  else if(x >= 750 )
    r5 = 0.3 * sin((M_PI * x) / 25) + 0.1 * sin((x*M_PI)/32) + 0.6 * sin((x*M_PI)/10);
  return r5;
}

double g2(int x1, int x2, int x3, int x4, int x5){
  double r6;
  r6 = ( x1 * x2 * x3 * x5 * ( x3 - 1 ) + x4 ) / (1 + pow(x1 , 2) + pow(x2 , 2));
  return r6;
}

void InitializeF(){
  double numData = 100;
  double* y = new double[105];
  double* u = new double[103];
  
  p = new double* [5];
  for(int i = 0 ; i < 5 ; i++)
    p[i] = new double[(int)numData];
  
  unsigned int biggestInteger = -1;
  srand(time(0));
    
  y[0] = ((double)(rand())/biggestInteger);
  y[1] = ((double)(rand())/biggestInteger);
  y[2] = ((double)(rand())/biggestInteger);
  
  for(int i = 3 ; i < numData + 5 ; i++){
    u[i-2] = g1(i+1);
    u[i-3] = g1(i+1);
    y[i] = g2(y[i-1] , y[i-2] , y[i-3] , u[i-2] , u[i-3]);
  }
  for(int i = 0 ; i < numData ; i++ ){
    p[0][i] = y[i+3];
    p[1][i] = y[i+4];
    p[2][i] = y[i+5];
    p[3][i] = g1(i+5);
    p[4][i] = g1(i+6);
  }
  
  NumDimofP = 5;
  NumDataofP = numData;
  
  T = new double[100]; 
  for(int i = 3 ; i < 102 ; i++ )
    T[i-3] = y[i];
  
  numNeuron = 2;
  epoch = 20;
}

double min(double* numbers , int indices){
  double min = numbers[0];
  for(int i = 0 ; i < indices ; i++)
    if(min > numbers[i]){
      min = numbers[i];
    }
  return min;    
}

double max(double* numbers , int indices){
  double max = numbers[0];
  for(int i = 0 ; i < indices ; i++)
    if(max < numbers[i]){
      max = numbers[i];
    }
  return max;    
}

double sum(double* numbers , int indices){
  double result = 0 ;
  //cout<<"numbers[0]"<<numbers[0]<<endl;
  for(int i = 0 ; i < indices ; i++)
    result += numbers[i];
  
  return result;
}

double Benchmark(int i_fn,int i_dim,double d_pos[])
{
	int i, j;
	double d_tmp[4];	// variables for temporary storage.
	double d_result = 0;	// cotains the result to return to the main program.
  
	/* initialize the parameters - Start */
	for(i=0;i<4;i++)
		d_tmp[i]=0.0;
	/* initialize the parameters - End   */
	if(i_fn==1)	// Anfis Function
	{
	  double z = 0.0001;
	  int numDim = NumDimofP;
	  int numData = NumDataofP;
	  //
#include <stdint.h>
	  uint32_t NumRule = pow(numNeuron , numDim);
	  //
	  
	  double **a;
	  a = new double *[NumRule];
	  for(int i = 0 ; i < NumRule ; i++)
	    a[i] = new double[numDim + 1];
	  
	  unsigned int biggestInteger = -1;
	  srand(time(0));
	    
	  // defination of a;
	  for(int i = 0 ; i < NumRule ; i++)
	    for(int j = 0 ; j < numDim + 1 ; j++)
	      a[i][j] = 2 * ((double)(rand())/biggestInteger) - 1;
	  //
	  double *lowc = new double[numDim];
	  double *highc = new double[numDim];
	  double *NC = new double[numDim];
	  
	  for(int i = 0 ; i < numDim ; i++){
	    lowc[i] = min(p[i] , numDim);
	      highc[i] = max(p[i] , numDim);
	      //cout<<"lowc:"<<lowc[i]<<" "<<"highc:"<<highc[i]<<endl;
	      NC[i] = (highc[i] - lowc[i]) / ( numNeuron + 1 );
	  }   
	  
	  // defination of c
	  double **c;
	  c = new double *[numNeuron];
	  for(int i = 0 ; i < numNeuron ; i++)
	    c[i] = new double[numDim];
	  
	  for(int i = 0 ; i < numNeuron ; i++)
	    for(int j = 0 ; j < numDim ; j++)
	      c[i][j] = lowc[j] + NC[j] * ( i + 1 );
	  //
	  // defination of s
	  double *s;
	  s = new double[numDim];
	  for(int i = 0 ; i < numDim ; i++)
	    s[i] = (highc[i] - lowc[i]) / pow(numNeuron , 2);
	  //
	  int max = numNeuron;
	  
	  // defination of Coff
	  double **Coff;
	  int index = pow(max,numDim);
	  Coff = new double *[index];
	  
	  for(int i = 0 ; i < index ; i++)
	    Coff[i] = new double[numDim];
	  for(int i = 0 ; i < numDim ; i++)
	    Coff[0][i] = 1;
	  int indexI = 0 ;
	  int indexK = 1;
	  int indexJ = numDim - 1;

	  // Calculating Coff
	  while(indexJ >= 0){
	    while( Coff[indexI][indexJ] < max ){
	      for(int counter = 0 ; counter < numDim ; counter++){
		Coff[indexK][counter] = Coff[indexI][counter];
		if(counter == indexJ)
		  Coff[indexK][counter] = Coff[indexI][indexJ] + 1;
	      }
	      indexK++;
	      indexI++;
	    }
	    indexI = 0;
	    indexJ--;
	  }
	  // end Calculating Coff
	  double error;
	  double e;
	  
	  // defination f1
	  double **f1;
	  f1 = new double*[numNeuron];
	  for(int i = 0 ; i < numNeuron ; i++)
	    f1[i] = new double[numDim];
	  
	  // defination f2
	  double* f2;
	  f2 = new double[(int)pow(numNeuron,numDim)];
	  
	  // defination f3;
	  double* f3;
	  f3 = new double[NumRule];
	  
	  // defination f4;
	  double* f4;
	  f4 = new double[NumRule];
	  
	  double sum1 = 0 ;
	  
	  // defination f5
	  double f5;
	  
	  
	  for(int ep = 0 ; ep < epoch ; ep++){
	    error = 0;
	    for(int j = 0 ; j < numData ; j++){
	      for(int i  = 0 ; i < numNeuron ; i++){
		for(int k = 0 ; k < numDim ; k++){
		  f1[i][k] = MF(p[k][j] , c[i][k] , s[k]);
		  //cout<<f1[i][k]<<" ";
		}
		//cout<<endl;
	      }
	      
	      for(int i = 0 ; i < pow(numNeuron,numDim) ; i++){
		f2[i] = 1;
		for(int k = 0 ; k < numDim ; k++){
		  f2[i] = f1[(int)Coff[i][k] - 1][k] * f2[i];
		}
	      }
	      
	      for(int i = 0 ; i < NumRule ; i++){
		//cout<<"man injam:" << sum(f2,pow(numNeuron,numDim));
		f3[i] = f2[i] / sum(f2,pow(numNeuron,numDim));
	      }
	      
	      
	      for(int i = 0 ; i < NumRule ; i++){
		sum1 = 0;
		int k;
		for( k = 0 ; k < numDim ; k++){
		  sum1 += a[i][k] * p[k][j];
		}
		k++;
		sum1 += a[i][k];
		f4[i] = f3[i] * sum1;
	      }
	      
	      //cout<<sum(f4,NumRule)<<endl;
	      f5 = sum(f4,NumRule);
	      
	      e = 0.5 * pow( (T[j] - f5 ), 2);
	      error = error + e;
	      for(int i = 0 ; i < NumRule ; i++){
		for(int k = 0 ; k < numDim ; k++)
		  a[i][k] = a[i][k] - z * d_pos[k];
	      }
	    }
	  }
	  d_result = pow ( ( (double)1 / numData ) * error , 0.5);
	}
	
	
	
	
	
	
	
	else if(i_fn==2)	// Sum of different power functions
	{
		for(j=0;j<i_dim;j++)
		{
			d_result += pow(abs(d_pos[j]),j+1);
		}
	}
	else if(i_fn==3)	// Ackley functions
	{
		for(j=0;j<i_dim;j++)
		{
			d_tmp[0] += pow(d_pos[j],2);
			d_tmp[1] += cos(d_PII*d_pos[j]);
		}
		d_tmp[0] = d_tmp[0]/i_dim;
		d_tmp[0] = -1 * 0.2 * sqrt(d_tmp[0]);
		
		d_tmp[1] = d_tmp[1]/i_dim;

	
		d_result= -1 * 20 * exp(d_tmp[0]) - exp(d_tmp[1]) + 20 + exp(1);
	}
/*
	else if(i_fn==4)	// Test Function 4
	{
	}
	else if(i_fn==5)  // Test Function 5
	{
	}
	else // Test Function 6
	{
	}
*/

  return d_result;
}