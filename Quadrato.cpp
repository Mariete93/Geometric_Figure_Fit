#include <cmath>
#include <iostream>
using namespace std;

double retta(double m, double q, double x)
{
	return m*x+q;
}

/*double evaluete_Chi2(double data[9][2], double error, double m, double q1, double q2)
{
	double Chi2=0;
	int npoints = sizeof(data)/sizeof(data[0][0])/2;
	cout << npoints/2 << "\n";
	for (int i=0; i<npoints/2+1; i++)
	{
		Chi2+=pow(data[i][1]-retta(m,q1,data[i][0]),2);
	}
	for (int j = npoints/2 ; j<npoints; j++)
	{
		Chi2+=pow(data[j][1]-retta(-1/m,q2,data[j][0]),2);
	}	
	Chi2=Chi2/pow(error,2);
	return Chi2/(npoints-3);
}

double evaluete_q_1(double data[9][2], double m)
{
	int npoints = sizeof(data)/sizeof(data[0][0])/2;
	cout << "Numero di punti" << npoints << "\n";
	double q = 0;
	for (int i=0; i<npoints/2+1; i++)
	{
		q+=data[i][0]*(data[i][1]-m*data[i][0]);
	}
	return 2*q/npoints;
}

double evaluete_q_2(double data[9][2], double m)
{
	int npoints = sizeof(data)/sizeof(data[0][0])/2;
	cout << "Numero di punti" << npoints << "\n";
	double q = 0;
	for (int j = npoints/2 ; j<npoints; j++)
	{
		q+=data[j][0]*(data[j][1]-m*data[j][0]);
	}
	return 2*q/npoints;
}*/


int main()
{
	cout << "Hello World\n";
	double points[9][2] = {{0,0},{10,10},{20,20},{30,30},{40,40},{30,50},{20,60},{10,70},{0,80}};// Matrix tha contains the coordinates of points
	double sigma = 2;	//Error in the determination of points 
	double m=1,q1=0,q2=80, err_q1=0, err_q2=0, err_m=0; //Parameters of the straights that identify the square
	double ang_coef=0;	//the best angular coefficent
	double Chi2=0;	//Cost Function
	int npoints = sizeof(points)/sizeof(points[0][0])/2;	//numbers of points
	//cout << npoints/2 << "\n";
	for (int i=0; i<npoints/2+1; i++)
	{
		Chi2+=pow(points[i][1]-retta(m,q1,points[i][0]),2);
	}
	for (int j = npoints/2 ; j<npoints; j++)
	{
		Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
	}	
	Chi2=0.5*Chi2/pow(sigma,2);
	cout << "Chi2 di riferimento " << Chi2/(npoints-2) << "\n";
	cout << "Numero di punti" << npoints << "\n";
	
	double a = 0, e=0,d=0,b=0, err_a = 0, err_b = 0, err_d =0, err_e=0; // coefficient of the 4th degree's equations for the determination of m
	double sum2=0;			// sum of the x^2 coordinates of points
	double sumxy=0;			//sum of the product of xand y coordinates
	double sumx=0;			// sum of the x coordinates
	double sumy=0;			// sum of the y coordinates
	//Calculation of a and b coefficient
	for (int i=0; i<npoints/2+1; i++)
	{
		sum2+=pow(points[i][0],2);
		sumxy+=points[i][0]*points[i][1];
		sumx += points[i][0];
		sumy += points[i][1];
	}
	a=sum2-2.0/npoints*pow(sumx,2);
	err_a = (1+2.0/npoints)*2*sumx*sigma;
	b =-(sumxy-2.0/npoints*sumx*sumy);
	err_b = (1+2.0/npoints)*(sumx+sumy)*sigma;

	//Calculation of d and e coefficient
	sum2 = 0;
	sumxy=0;
	sumx=0;
	sumy=0;
	for (int j = npoints/2 ; j<npoints; j++)
	{
		sum2+=pow(points[j][0],2);
		sumxy+=points[j][0]*points[j][1];
		sumx += points[j][0];
		sumy += points[j][1];
	}
	d =(sumxy-2.0/npoints*sumx*sumy);
	err_d = -(1+2.0/npoints)*2*sumx*sigma;
	e=(sum2-2.0/npoints*pow(sumx,2));
	err_e = (1+2.0/npoints)*(sumx+sumy)*sigma;

	//definition of coefficent for the resolution of 4th degree equation
	double err_s =0, err_p=0, err_r=0, err_Delta=0, err_Q = 0, err_S=0, err_argQ=0;
	double s = 27*(a*d*d+b*b*e);
	err_s=27*(err_a*d*d+2*err_d*a*d+err_e*b*b+2*err_b*e*b);
	double p = -3.0/8.0*b*b/(a*a);
	err_p = -6.0/8.0*(err_b*b+b*b/(a*a*a)*err_a);
	double r= 12*a*e-3*b*d;
	err_r = 12*(err_a*e+a*err_e+3*err_b*d+3*err_d*b);
	double Delta = pow(0.5*(s+sqrt(s*s-4*pow(r,3))),1.0/3.0);
	err_Delta = 1.0/3.0*(err_s+(s*err_s+6*err_r*r*r)/sqrt(s*s-4*pow(r,3)))/pow(0.5*(s+sqrt(s*s-4*pow(r,3))),2.0/3.0);
	double argQ = -2.0/3.0*p+1.0/(3.0*a)*(Delta+r/Delta);
	double Q = 0;
	double S = d/a+pow(b/(2.0*a),3);
	err_S = err_d/a+d/(a*a)*err_a+3*(err_b+err_a*b)*pow(b/(2.0*a),2)/(2.0*a); 


	double x[4] = {0,0,0,0};	//array that cointain the solution of equations
	double err_x[4]={0,0,0,0};	//array that cointain the errors on the solution of equations
	cout << "argQ " << argQ; 
	if (argQ > 0)
	{	
		err_argQ = -2.0/3.0*err_p+1.0/(3.0*a)*(err_a/a*(Delta+r/Delta)+(err_Delta+err_r/Delta+r*err_Delta/Delta*Delta));
		Q = 0.5*sqrt(argQ);
		err_Q=0.25*err_argQ/sqrt(argQ);
		cout << "\nQ=" << Q << "\nDelta_1 = " << -4*Q*Q-2*p+S/Q << "\nDelta_2 " << -4*Q*Q-2*p-S/Q << "\n";  
		if (-4*Q*Q-2*p+S/Q > 0)
		{
			//cout << "Delta_1 positivo \n";
			x[0] = -b/(4*a)-Q+0.5*sqrt(-4*Q*Q-2*p+S/Q);
			err_x[0] = (err_b+b*err_a)/(4.0*a)+err_Q+0.25*(8*Q*err_Q+2*err_p+err_S/Q+S*err_Q/(Q*Q))/sqrt(-4*Q*Q-2*p+S/Q);
			x[1] = -b/(4*a)-Q-0.5*sqrt(-4*Q*Q-2*p+S/Q);
			err_x[1]=err_x[0];
		}
		if ( -4*Q*Q-2*p-S/Q > 0)
		{
			//cout << "Delta_2 positivo \n";
			x[2] = -b/(4*a)+Q+0.5*sqrt(-4*Q*Q-2*p-S/Q);
			err_x[2] = (err_b+b*err_a)/(4.0*a)+err_Q+0.25*(8*Q*err_Q+2*err_p+err_S/Q+S*err_Q/(Q*Q))/sqrt(-4*Q*Q-2*p-S/Q);
			x[3] = -b/(4*a)+Q-0.5*sqrt(-4*Q*Q-2*p-S/Q);
			err_x[3] = err_x[2];
		}
		cout << "x0=" << x[0] << "\n";
		cout << "x1=" << x[1] << "\n";
		cout << "x2=" << x[2] << "\n";
		cout << "x3=" << x[3] << "\n";
		cout << "Risultato\n";
		double Chi2Min=1000000;
		for (int k=0; k<4; k++)
		{
			if (x[k] != 0)
			{ 
				m=x[k];
				sumx=0;
				sumy=0;
				for (int i=0; i<npoints/2+1; i++)
				{	
					sumx+=points[i][0];
					sumy+=points[i][1];
				}
				q1= 2/npoints*(sumy-m*sumx);
				err_q1=2/npoints*(npoints*sigma+err_x[k]*sumx+m*npoints*sigma);
				sumx=0;
				sumy=0;
				for (int j = npoints/2 ; j<npoints; j++)
				{
					sumx+=points[j][0];
					sumy+=points[j][1];
				}
				q1= 2/npoints*(sumy+sumx/m);
				err_q2=2/npoints*(npoints*sigma+err_x[k]*sumx/(m*m)+npoints*sigma/m);
				Chi2 = 0;
				for (int i=0; i<npoints/2+1; i++)
				{
					Chi2+=pow(points[i][1]-retta(m,q1,points[i][0]),2);
				}
				for (int j = npoints/2 ; j<npoints; j++)
				{
					Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
				}	
				Chi2=Chi2/(2*pow(sigma,2));
				Chi2 = Chi2/(npoints-2);
				if (Chi2 < Chi2Min)
				{
					Chi2Min = Chi2;
					ang_coef= m;
					err_m = err_x[k];
				}
			}
		}
		cout <<"Coefficiente angolare: " << ang_coef << " +- " << abs(err_m) << "\nChi2= " << Chi2Min << "\n";
	}
}

