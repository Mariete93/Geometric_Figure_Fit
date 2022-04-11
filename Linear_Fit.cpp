#include <cmath>
#include <iostream>
using namespace std;

double retta(double m, double q, double x)
{
	return m*x+q;
}

int main()
{
	cout << "Hello World\n";
	double points[9][2] = {{0,0},{100,101},{200,202},{300,290},{400,402},{500,298},{600,201},{700,99},{800,2}};// Matrix tha contains the coordinates of points
	double sigma = 4;	//Error in the determination of points 
	double m=1,q=0,  q2=0, err_q=0,err_m=0; //Parameters of the straights that identify the square
	double N = 0; //numbers of points per edge
	double ang_coef=0, err_ang_coeff=0;	//the best angular coefficent
	double Chi2=0;	//Cost Function
	int npoints = sizeof(points)/sizeof(points[0][0])/2;	//numbers of points
	
	//Calculation of a and b coefficient
	double Chi2Min = 10000000;
	for (int k=0; k<2; k++)	//loop on the edge of square
	{	int start=0, end=0;	//start and ends of the sum
		double sum2=0;			// sum of the x^2 coordinates of points
		double sumxy=0;			//sum of the product of xand y coordinates
		double sumx=0;			// sum of the x coordinates
		double sumy=0;			// sum of the y coordinates
		if (k==0)
		{
			end=npoints/2+1;	
		}
		else
		{
			start=npoints/2;
			end=npoints;
		}
		for (int i=start; i<end; i++)	//sums
		{
			sum2+=pow(points[i][0],2);
			sumxy+=points[i][0]*points[i][1];
			sumx += points[i][0];
			sumy += points[i][1];
		}
		N=end-start;
    	m = (N*sumxy-sumx*sumy)/(N*sum2-pow(sumx,2));	//angular coefficent
    	q = (sumy-m*sumx)/N; 									//bias first edge
		q2 = points[4][1]+1.0/m*points[4][0]; 						//bias second edge
    	err_m = 2*sigma*N*sumy/(N*sum2-pow(sumx,2));	//error on m
    	err_q= sigma+err_m*sumx/N;							//error on q
		
    	Chi2=0;	
		if (k==0)
		{													//loss function
			for (int i=start; i<end; i++)
			{	
				int j=i+end;
				Chi2+=pow(points[i][1]-retta(m,q,points[i][0]),2);
				Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
			}	
			Chi2=Chi2/pow(sigma,2);
		}
		else
		{
			for (int i=start; i<end; i++)
			{	
				int j=i-start;
				Chi2+=pow(points[i][1]-retta(m,q,points[i][0]),2);
				Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
			}	
			Chi2=Chi2/pow(sigma,2);	
		}
		cout <<"Coefficiente angolare: " << m << " +- " << abs(err_m) << " q = " << q << " q_2 = " << q2 <<  "\nChi2= " << Chi2/(N-2) << "\n";
		if (Chi2<Chi2Min)										//choose the best fit
		{
			Chi2Min = Chi2;
			ang_coef= m;
			err_ang_coeff = err_m;
		}
	}
	cout <<"Coefficiente angolare: " << ang_coef << " +- " << abs(err_ang_coeff) << "\nChi2= " << Chi2/(N-2) << "\n";
}

