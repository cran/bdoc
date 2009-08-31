// bdoc.cpp : main project file.

//#include "stdafx.h"

//using namespace System;

//#include "stdafx.h"
//#include <iostream>
#include "R.h"
#include "Rmath.h"
//#include "stdio.h"
//using namespace std;




void bdoc1( int *spec, int *base, int *pos, int *nrowtr, double *epsilon, double *prior, double *adjlike2, int *bar, int *trdata, double *posteriors, int *q, int *bar2)
{
	
	
	int i,j,k,z,sr,count;     
	
	double spec2=*spec;
	double y[*pos][*base][*spec];  //y[pos][bases][spec]
	int y2[*pos][*nrowtr];  //y2[pos][nrowtr]
	double post[*spec];  //post[spec]
	double postdenom;
	double postmat[*spec][*pos];  //postmat[spec][pos]
	

	//this assigns all of the adjusted likelihoods to the 3-d array y
	count=0;
		for(j=0;j<*pos;j++)
		{
			for(k=0;k<*base;k++)
			{
				for(i=0;i<*spec;i++)
				{
					y[j][k][i]=adjlike2[count];
					count++;
				}
			}
		}
	//this assigns all of the training data to a matrix
	count=0;
		for(j=0;j<*pos;j++)
			{
				for(i=0;i<*nrowtr;i++)  
				{
					y2[j][i]=trdata[count];
					count++;
				}
			}

	//for each position, this calculates the posterior probability
	for(j=0;j<*pos;j++)
	{

		//quick check to see if nucleotide in pos j for test barcode appears in pos j for train data;
		z=0;                              
		for(i=0;i<*nrowtr;i++)  
		{
			if(bar[j]==y2[j][i])      //NEED TO GO BACK AND CHANGE TO == WHEN YOU GET THE CONDITION SORTED OUT
			{
				z++;
			}
			if(z>0) break;
		}

		//moving the posteriors toward 1/s if there are no matches
		sr=0; 
		if(z==0)          //NEEDS TO BE CHANGED TO == WHEN YOU GET IT WORKING
		{
			for(i=0;i<*spec;i++)
			{
				post[i]=prior[i]+(1/spec2-prior[i])*(*epsilon);
				if(post[i]==1) sr=1;                             //this will help with triggering the stopping rule
			}
		}

		else
		{								
			//computing posterior at postion j
			postdenom=0;
			if(bar[j]==65)   //A 
			{
				for(i=0;i<*spec;i++)
				{
					postdenom+=prior[i]*y[j][0][i];
				}
				for(i=0;i<*spec;i++)
				{
					post[i]=prior[i]*y[j][0][i]/postdenom;
					if(post[i]==1) sr=1;
				}
			}

			else if(bar[j]==84)  //T 
			{
				for(i=0;i<*spec;i++)
				{
					postdenom+=prior[i]*y[j][1][i];
				}
				for(i=0;i<*spec;i++)
				{
					post[i]=prior[i]*y[j][1][i]/postdenom;
					if(post[i]==1) sr=1;
				}
			}
			
			else if(bar[j]==67) //C
			{
				for(i=0;i<*spec;i++)
				{
					postdenom+=prior[i]*y[j][2][i];
				}
				for(i=0;i<*spec;i++)
				{
					post[i]=prior[i]*y[j][2][i]/postdenom;
					if(post[i]==1) sr=1;
				}
			}

			else if(bar[j]==71)   //G
			{
				for(i=0;i<*spec;i++)
				{
					postdenom+=prior[i]*y[j][3][i];
				}
				for(i=0;i<*spec;i++)
				{
					post[i]=prior[i]*y[j][3][i]/postdenom;
					if(post[i]==1) sr=1;
				}
			}
			else 
			{
				for(i=0;i<*spec;i++)
				{                                        //IT LOOKS LIKE THE FUNCTION IS ONLY DOING THIS CALCULATION.  PERHAPS BAR[J]!= A OR T OR C OR G... NEED TO LOOK INTO
					post[i]=prior[i];
					
				}
			}
		}											
		
		 
		

		
		
		for(i=0;i<*spec;i++)
		{
			prior[i]=post[i];
		}
		

		/*this is a matrix to keep track of the posterior probabilities*/
		for(i=0;i<*spec;i++)
		{
			postmat[i][j]=post[i];
		}
		//if(*max_element(post,post+spec)==1) { break;}
		//sr=0;
		//for(i=0;i<*spec;i++)
		//{
		//	if(post[i]>=.9) sr=1; break;
		//}
		//if(sr==1)  break;   
	}
	*q=j;
  count=0;
  for(k=0;k<j;k++)
  {
	  for(i=0;i<*spec;i++)
	  {
	  	  posteriors[count]=postmat[i][k];
	  	  count++;
	  }
  }
  
  for(i=0;i<*pos;i++)
  {
  	 bar2[i]=bar[i];
  }
}



