#include <cmath>
#include <iostream>
using namespace std;

double retta(double m, double q, double x)
{
	return m*x+q;
}

double evaluete_Chi2(double data[9][2], double error, double m, double q1, double q2)
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
}


int main()
{
	cout << "Hello World\n";
	double points[9][2] = {{0,0},{11,10},{20,22},{30,31},{39,41},{30,50},{21,58},{11,72},{0,81}};
	double sigma = 2;
	double m=1,q1=0,q2=80;
	double ang_coef=0;
	double Chi2=0;
	int npoints = sizeof(points)/sizeof(points[0][0])/2;
	cout << npoints/2 << "\n";
	for (int i=0; i<npoints/2+1; i++)
	{
		Chi2+=pow(points[i][1]-retta(m,q1,points[i][0]),2);
	}
	for (int j = npoints/2 ; j<npoints; j++)
	{
		Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
	}	
	Chi2=Chi2/pow(sigma,2);
	cout << "Chi2 di riferimento " << Chi2/(npoints-2) << "\n";
	cout << "Numero di punti" << npoints << "\n";
	
	double a = 0, e=0,d=0,b=0;
	double sum =0;
	double sum2=0;
	double sumxy=0;
	double sumx=0;
	double sumy=0;
	for (int i=0; i<npoints/2+1; i++)
	{
		sum2+=pow(points[i][0],2);
		sum += points[i][0];
		sumxy+=points[i][0]*points[i][1];
		sumx += points[i][0];
		sumy += points[i][1];
	}
	a=sum2-2/npoints*pow(sum,2);
	b =-(sum-2/npoints*sumx*sumy);
	sum2 = 0;
	sum = 0;
	sumxy=0;
	sumx=0;
	sumy=0;
	for (int j = npoints/2 ; j<npoints; j++)
	{
		sum2+=pow(points[j][0],2);
		sum += points[j][0];
		sumxy+=points[j][0]*points[j][1];
		sumx += points[j][0];
		sumy += points[j][1];
	}
	d =-(sum-2/npoints*sumx*sumy);
	e=-(sum2-2/npoints*pow(sum,2));

	double s = 27*(a*d*d+b*b*e);
	double p = -3.0/8.0*b*b/(a*a);
	double r= 12*a*e-3*b*d;
	double Delta = pow(0.5*(s+sqrt(s*s-4*pow(r,3))),1.0/3.0);
	double argQ = -2.0/3.0*p+1.0/(3.0*a)*(Delta+r/Delta);
	double Q = 0;
	double S = d/a+pow(b/(2.0*a),3);
	double x[4] = {0,0,0,0};
	cout << "argQ " << argQ; 
	if (argQ > 0)
	{	
		Q = 0.5*sqrt(argQ);
		cout << "\nQ=" << Q << "\nDelta_1 = " << -4*Q*Q-2*p+S/Q << "\nDelta_2 " << -4*Q*Q-2*p-S/Q << "\n";  
		if (-4*Q*Q-2*p+S/Q > 0)
		{
			cout << "Delta_1 positivo \n";
			x[0] = -b/(4*a)-Q+0.5*sqrt(-4*Q*Q-2*p+S/Q);
			x[1] = -b/(4*a)-Q-0.5*sqrt(-4*Q*Q-2*p+S/Q);
		}
		if ( -4*Q*Q-2*p-S/Q > 0)
		{
			cout << "Delta_2 positivo \n";
			x[2] = -b/(4*a)+Q+0.5*sqrt(-4*Q*Q-2*p-S/Q);
			x[3] = -b/(4*a)+Q-0.5*sqrt(-4*Q*Q-2*p-S/Q);
		}
		cout << "x0=" << x[0] << "\n";
		cout << "x1=" << x[1] << "\n";
		cout << "x2=" << x[2] << "\n";
		cout << "x3=" << x[3] << "\n";
		cout << "Risultato\n";
		double Chi2Min=1000000;
		for (int k=0; k<4; k++)
		{
			if (x[k] == 0)
			{ 
				cout << "soluzione immaginaria\n";
			}
			else
			{ 
				m=x[k];
				double q = 0;
				for (int i=0; i<npoints/2+1; i++)
				{
					q+=(points[i][1]-m*points[i][0]);
				}
				q1= 2*q/npoints;
				//q1 = evaluete_q_1(points, m);
				q=0;
				for (int j = npoints/2 ; j<npoints; j++)
				{
					q+=(points[j][1]+points[j][0]/m);
				}
				q2 = 2*q/npoints;
				//q2= evaluete_q_2(points,-1/m);
				Chi2 = 0;
				for (int i=0; i<npoints/2+1; i++)
				{
					Chi2+=pow(points[i][1]-retta(m,q1,points[i][0]),2);
				}
				for (int j = npoints/2 ; j<npoints; j++)
				{
					Chi2+=pow(points[j][1]-retta(-1/m,q2,points[j][0]),2);
				}	
				Chi2=Chi2/pow(sigma,2);
				Chi2 = Chi2/(npoints-2);
				if (Chi2 < Chi2Min)
				{
					Chi2Min = Chi2;
					ang_coef= m;
				}
			}
		}
		cout <<"Coefficiente angolare: " << ang_coef << "\nChi2= " << Chi2Min << "\n";
	}
}

