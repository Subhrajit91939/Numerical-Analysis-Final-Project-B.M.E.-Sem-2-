#include<stdio.h>
#include<math.h>
#include<stdlib.h>

//Acceleration due to gravity
#define g 9.81
#define SMALL 1e-6
#define LARGE 3000

#define h 0.01

//Trajectory Equation
double Org(double x,double alpha,double v0)
{
	double j= 2*v0*v0*(pow((cos(alpha)),2));
	return ((x*(tan(alpha)))-(g*(x*x/j)));
}

//Trajectory Equation in Differential form
double F(double x, double y,double alpha,double v0)
{
	double m = v0*v0*(pow((cos(alpha)),2));
	return (tan(alpha) - (g*(x/m)));
}

//Time taken for the particle 
double time(double x,double alpha,double v0)
{
	return x/(v0*cos(alpha));
	
}

double height(double x,double alpha,double v0)
{
	return (v0*v0*sin(alpha)*sin(alpha))/(2*g);
	
}
double range(double alpha,double v0)
{
	return ((v0*v0*sin(2*alpha))/g);
}

//Typedef Statement for the array returned by the main Function Block
typedef double (*doublepointer)[LARGE];

//euler method
doublepointer euler(double x0, double y0, double x_f,double alpha,double v0)
{
	double y;
	static double arr1[LARGE][LARGE];
	int i=1;
	arr1[0][0] = x0;
	arr1[0][1] = y0;
	printf("\nEuler Method for Solving ODE:\n\n");
	printf("X\t\ty\t\tY(direct)\n");
	printf("-------\t\t-------\t\t-------\n");
	while(fabs(x_f-x0)>SMALL)
	{
		y = y0 + h*F(x0, y0,alpha,v0);
		y0 = y;
		arr1[i][1] = y0;
		x0 = x0 + h;
		arr1[i][0] = x0;
		double z=Org(x0,alpha,v0);
		printf("%.5lf\t%.20lf\t%.20lf\n", x0, y,z);
		i++;
	}
	//return y;
    return arr1;
    //return 0;
}

//RK2
doublepointer runge_kutta2(double x0, double y0, double x_f,double alpha,double v0)
{
	double y,k1,k2;
	static double arr2[LARGE][LARGE];
	int j=1;
	arr2[0][0] = x0;
	arr2[0][1] = y0;
	printf("\nRunge-Kutta-2 for Solving ODE:\n\n");
	printf("X\t\ty\t\tY(direct)\n");
	printf("-------\t\t-------\t\t-------\n");
    while(fabs(x_f-x0)>SMALL)
	{
        k1=h*F(x0,y0,alpha,v0);
        k2=h*F(x0+h/2.0,y0+k1/2.0,alpha,v0);
        y=y0+k2;
        y0=y;
        arr2[j][1] = y0;
        x0=x0+h;
        arr2[j][0] = x0;
        double z=Org(x0,alpha,v0);
        printf("%.5lf\t%.20lf\t%.20lf\n", x0, y, z);
        j++;
    }
    return arr2;
    //return 0;
}

//RK4
doublepointer runge_kutta4(double x0, double y0, double x_f,double alpha,double v0)
{
	double y,k1,k2,k3,k4;
	static double arr3[LARGE][LARGE];
	int k=1;
	arr3[0][0] = x0;
	arr3[0][1] = y0;
	printf("\nRunge-Kutta-4 for Solving ODE:\n\n");
	printf("X\t\ty\t\tY(direct)\n");
	printf("-------\t\t-------\t\t-------\n");
    while(fabs(x_f-x0)>SMALL)
	{
        k1=h*F(x0,y0,alpha,v0);
        k2=h*F(x0+h/2.0,y0+k1/2.0,alpha,v0);
        k3=h*F(x0+h/2.0,y0+k2/2.0,alpha,v0);
        k4=h*F(x0+h,y0+k3,alpha,v0);
        y=y0+(1/6.0)*(k1+2*k2+2*k3+k4);
        y0=y;
        arr3[k][1] = y0;
        x0=x0+h;
        arr3[k][0] = x0;
        double z=Org(x0,alpha,v0);
        printf("%.5lf\t%.20lf\t%.20lf\n", x0, y, z);
        k++;
    }
    printf("\n");
    return arr3;
    //return 0;
}
//main function
int main(void)
{
	//file pointers for Displaying Graphs using Python matplotlib
	FILE *eulerout, *rk2out, *rk4out;
	eulerout = fopen("euler_output1.csv", "w");
	rk2out = fopen("rk2_output1.csv", "w");
	rk4out = fopen("rk4_output1.csv", "w");
	
	int c, t, Ra=36;
	//We have taken the distance to be 0.3m.
	double x0=0, y0=0, x=0.360,alpha,z,v0;
	double (*m1)[LARGE], (*m2)[LARGE], (*m3)[LARGE];
	//user input for alpha , v0
	printf("Enter the angle of projection in degrees:");
	scanf("%lf", &alpha);
	alpha = alpha*M_PI/180;
	
	printf("Enter the projected velocity V0:");
	scanf("%lf", &v0);
	printf("The X-coordinate taken is 0.3m\n");
	printf("Press 1 For Euler Method:\nPress 2 For Runge Kutta Order 2:\nPress 3 For Runge Kutta Order 4:\nPress 4 For Comparision Among All Three:\n\n ");
	printf("Enter your choice:");
	
	scanf("%d", &c);
	switch(c)
	{
		case 1:
			/*Euler Method*/
			
			m1 = euler(x0, y0, x,alpha,v0);
			
			t=0;
			while(t<=Ra)
			{
				fprintf(eulerout, "%.5lf,%.20lf\n", m1[t][0], m1[t][1]);
				t++;
			}
			
			printf("The time taken for the particle to reach at given x:\n");
			printf("%lf\n",time(x,alpha,v0));
			printf("The maximum height the particle can reach x:\n");
			printf("%lf\n",height(x,alpha,v0));
			printf("the maximum range covered by the particle :\n");
			printf("%lf\n",range(alpha,v0));
			break;
		
		case 2:
			/*Runge-Kutta-2 Method for Solving*/
			
			m2 = runge_kutta2(x0, y0, x,alpha,v0);
			
			t=0;
			while(t<=Ra)
			{
				fprintf(rk2out, "%.5lf,%.20lf\n", m2[t][0], m2[t][1]);
				t++;
			}
			
			printf("The time taken for the particle to reach at given x:\n");
			printf("%lf\n",time(x, alpha,v0));
			printf("The maximum height the particle can reach x:\n");
			printf("%lf\n",height(x,alpha,v0));	
			printf("the maximum range covered by the particle :\n");
			printf("%lf\n",range(alpha,v0));		
			break;
			
		case 3:
			/*Runge-Kutta-4 Method for Solving*/
		
			m3 = runge_kutta4(x0, y0, x, alpha,v0);
			
			t=0;
			while(t<=Ra)
			{
				fprintf(rk4out, "%.5lf,%.20lf\n", m3[t][0], m3[t][1]);
				t++;
			}
			
			printf("The time taken for the particle to reach at given x:\n");
			printf("%lf\n",time(x,alpha,v0));
			printf("The maximum height the particle can reach x:\n");
			printf("%lf\n",height(x,alpha,v0));
			printf("the maximum range covered by the particle :\n");
			printf("%lf\n",range(alpha,v0));						
			break;
			
		case 4:
			/*Comparison between All the Methods*/
			printf("COMPARISON\n");
			m1 = euler(x0, y0, x,alpha,v0);
			m2 = runge_kutta2(x0, y0, x,alpha,v0);
			m3 = runge_kutta4(x0, y0, x,alpha,v0);
			
			
			t=0;
			while(t<=Ra)
			{
				fprintf(eulerout, "%.5lf,%.20lf\n", m1[t][0], m1[t][1]);
				t++;
			}
			
			t=0;
			while(t<=Ra)
			{
				fprintf(rk2out, "%.5lf,%.20lf\n", m2[t][0], m2[t][1]);
				t++;
			}
			
			t=0;
			while(t<=Ra)
			{
				fprintf(rk4out, "%.5lf,%.20lf\n", m3[t][0], m3[t][1]);
				t++;
			}
			
			printf("The time taken for the particle to reach at given x:\n");
			printf("%lf\n",time(x,alpha,v0));
			printf("The maximum height the particle can reach x:\n");
			printf("%lf\n",height(x,alpha,v0));	
			printf("the maximum range covered by the particle :\n");
			printf("%lf\n",range(alpha,v0));					
			//printf("%lf\t%lf\t%lf\n", m1, m2, m3);
			break;
		
		default:
			printf("Invalid Input.\n");
	}
	
	fclose(eulerout);
	fclose(rk2out);
	fclose(rk4out);
	return 0;
}
