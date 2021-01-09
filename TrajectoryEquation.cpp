#include<stdio.h>
#include<math.h>
#include<stdlib.h>
//Acceleration due to gravity
#define g 9.81
#define SMALL 1e-7
#define LARGE 3000

#define h 0.001

//Typedef Statement for the array returned by the main Function Block
typedef double (*doublepointer)[LARGE];

//Trajectory Equation
doublepointer traject(double x,double alpha,double v0)
{
	static double arr[LARGE][LARGE];
	double y;
	double j= 2*v0*v0*(pow((cos(alpha)),2));
	int i;
	arr[0][1] = 0;
	arr[0][0] = 0;
	i=0;
	do{
		y = ((x*(tan(alpha)))-(g*(x*x/j)));
		arr[i][1] = y;
		x += h;
		arr[i][0] = x;
		i++;
	}while(x<=0.41);
	return arr;
}

int main(void)
{
	FILE *output;
	output = fopen("dtrajectory45deg.csv", "w");
	//1:30deg
	//2:45deg
	//3:60deg
	double x0, alpha=45, v0=2;
	int t=0;
	double (*m1)[LARGE];
	x0 = 0;
	alpha = alpha*M_PI/180;
	m1 = traject(x0, alpha, v0);
	while(t<=420)
	{
		fprintf(output, "%.5lf,%.20lf\n", m1[t][0], m1[t][1]);	
		t++;
	}
	fclose(output);
	return 0;
}
