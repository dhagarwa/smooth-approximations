#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "symbolicc++.h"
#include <Eigen/Dense>

#define max_M 15
#define min_M 1
#define d 2
#define alpha 0.1
#define ERROR_GRID 100

using namespace Eigen;
using namespace std;
 

double **get_points_random(int N, int n, int M , int row, int column)
{
	double **points = (double **)malloc(N * sizeof(double *));
	int i, j, k;
	double a, b, c;
	float r1, r2;
	for (i = 0; i < N; i++) 
	{
        	points[i] = (double *)malloc(d * sizeof(double));
	}
	a = (column-1.0)/M; b = column/M;
	c = (row -1.0)/M; 
	for(k=1; k <= N; k++)
	{
		r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		points[k-1][0] = a + r1/M; 
		points[k-1][1] = c + r2/M;
	}
	return points;
}

double **get_points_regular(int N, int n, int M, int row, int column)
{

	double **points = (double **)malloc(N * sizeof(double *));
	int i, j, k;
	double diff; 
	for (i = 0; i < N; i++) 
	{
        	points[i] = (double *)malloc(d * sizeof(double));
	}
	double a, b, c;
	a = (column-1.0)/M; b = column/M;
	c = (row -1.0)/M; 
	if(n == 0)
		diff =0.0;
	else
		diff = 1.0/(n*M);
	//cout << "Diff = " << diff << endl;
	k = 1;
	for(i=0; i<=n; i++)
	{
		for(j=0; j <= n - i; j++)
		{
			points[k-1][0] = a + i * diff;
			points[k-1][1] = c + j * diff;
			k++;
		}
	}
	return points;
}

//get_values calculates values of function f to be interpolated at interpolation nodes

double *get_values(int N, double **points)
{
	float delta;
	int i;
	double *values;
	values = (double *)malloc(N * sizeof(double));
	for(i=1; i <= N; i++)
	{
		//values[i-1] = 2*points[i-1][0] + 3*points[i-1][1] +  points[i-1][1]* points[i-1][1] + 1; // f = 2x + 3y + y^2 + 1. Need to find way to represent using SymbolicC++. 
		values[i-1] = sin(points[i-1][0]) + sin(points[i-1][1]) ;  //f = sin(x) + sin(y)
		//values[i-1] = exp(points[i-1][0] + points[i-1][1]); // f= exp(x+y)
		//values[i-1] = (points[i-1][0] + points[i-1][1]) * sin(points[i-1][0]); //f = (x+y)sin(x)
	}
	return values;
}



Symbolic get_mollifier(Symbolic var, double a, double b )
{
	Symbolic mollifier("mollifier");
	//cout << "hello " << var << endl;
	//mollifier = (var^2);
	mollifier = (exp(1/(1 - ((var - a - 1)^2))))/ (exp(1/(1 - ((var - a -1)^2))) + exp(1/(1 - ((b + 1 - var)^2))));
	return mollifier;
}

Symbolic weave_1( Symbolic f, Symbolic g, double a , double b)
{
	Symbolic weaved_func("weaved_func"), mollifier("mollifier"), x("x");
	mollifier = get_mollifier(x, a, b);
	weaved_func = (1 - mollifier)*f + mollifier*g;
	return weaved_func;
}

Symbolic weave_2(Symbolic f, Symbolic g, double a, double b)
{
	Symbolic weaved_func("weaved_func"), mollifier("mollifier"), y("y");
	mollifier = get_mollifier(y, a, b);
	//cout << mollifier << endl;
	weaved_func = (1 - mollifier)*f + mollifier*g;
	return weaved_func;
}


Symbolic vandermonde(int N, int n, int l, double **points)
{
	Symbolic V("V", N, N), x("x"), y("y"), det("det"), z;
	int i, j, power_sum, x_power, y_power;
	//z = (x^3) + (y^3); cout << z << endl;
	for(j = 1; j <= N; j++)
	{
		power_sum = 0;
		i = 1;
		for(power_sum =0; power_sum <=n; power_sum++)
		{
			for(x_power=0; x_power<= power_sum; x_power++)
			{
				y_power = power_sum - x_power;
				if (j == l)
				{
					V(i-1, j-1) = (Symbolic(x)^x_power) * (Symbolic(y)^y_power);
					//cout << V(i-1, j-1);
				}
				else
				{
					V(i-1, j-1) = pow(points[j-1][0],x_power) * pow(points[j-1][1], y_power);
					//cout << V(i-1, j-1) << endl;
				}
				i++;
			}
		}
	}


	//if(l==2) cout << V;
	det = V.determinant();
	//if (l==2) cout << det;
	return det;
}


Symbolic interpolation(int N, int n, double **points, double *values)
{
	Symbolic interp_poly("interp_poly"), den("den"), num("num");
	interp_poly = 0;
	int i = 1;
	den = vandermonde(N, n, 0, points);
	//cout << den;
	for(i = 1; i <= N; i++)
	{
		num = vandermonde(N, n, i, points);
		//cout << num << endl;
		interp_poly += values[i-1] * num/ den;
		//cout << interp_poly << endl;
	}

	//cout << interp_poly;
	return interp_poly;
}

Symbolic interpolation_fast(int N, int n, double **points, double *values)
{
	Symbolic interp_poly("interp_poly"), x("x"), y("y");
	interp_poly = 0;
	MatrixXd V(N,N);
	VectorXd val(N), sol(N);
	int i, j, power_sum, x_power, y_power;
	for(j = 1; j <= N; j++)
	{
		power_sum = 0;
		i = 1;
		for(power_sum =0; power_sum <=n; power_sum++)
		{
			for(x_power=0; x_power<= power_sum; x_power++)
			{
				y_power = power_sum - x_power;
				V(j-1, i-1) = pow(points[j-1][0],x_power) * pow(points[j-1][1], y_power);
				//cout << V(i-1, j-1) << endl;
				i++;
			}
		}
		val(j-1) = values[j-1];
	}  

	sol = V.partialPivLu().solve(val);
	power_sum = 0;
	i = 1;
	for(power_sum =0; power_sum <=n; power_sum++)
	{
		for(x_power=0; x_power<= power_sum; x_power++)
		{
			y_power = power_sum - x_power;
			interp_poly += sol(i-1)* (x^x_power) * (y^y_power);
			i++;
		}
	}
  	//cout << "Left hand matrix: "<< V << endl;
	//cout << "Right hand vector: " << val << endl;
	//cout << "Solution vector first element: " << sol(0) << endl;
	//cout << "Interpolation polynomial: " << interp_poly << endl;
	return interp_poly;
}


void calculate_error(Symbolic interp_func, int M)
{
	Symbolic x("x"), y("y");
	//cout << interp_func << endl;
	double	max_error = 0.0;
	double c1, c2, error, max_c1, max_c2;
	double ***error_grid_points = (double ***)malloc(ERROR_GRID * sizeof(double **));
	int i, j, row, column, third_index, max_row, max_column;
	for (i = 0; i < ERROR_GRID ; i++) 
	{
		error_grid_points[i] = (double **)malloc(ERROR_GRID * sizeof(double *));
		for(j = 0; j < ERROR_GRID; j++)
		{
        		error_grid_points[i][j] = (double *)malloc(d * sizeof(double));
		}
	}	
	for(i =0; i < ERROR_GRID; i++)
	{
		for(j =0; j < ERROR_GRID; j++)
		{
			error_grid_points[i][j][0] = j/ (ERROR_GRID - 1.0);
			error_grid_points[i][j][1] = i/ (ERROR_GRID - 1.0);  
		}
	}
	
	for(i=0; i < ERROR_GRID - 1; i++)
	{
		for(j=0; j < ERROR_GRID -1; j++)
		{
			c1 = error_grid_points[i][j][0]; c2 = error_grid_points[i][j][1];
			column = floor(c1 * M) + 1;
			row = floor(c2 * M) + 1;
			//cout << "x = " << c1 << " y= " << c2; 
			//cout << " Row: " << row << " Column: " << column << endl;
			//cout << "poly: " << interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4) << endl; //<< " value: "<< interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4)[x == c1, y == c2] << endl;
			//error = abs(double(interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4)[x == c1, y == c2]) - (2*c1 + 3*c2 +  c2*c2 + 1));
			error = abs(double(interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4)[x == c1, y == c2]) - (sin(c1) + sin(c2)));
			//error = abs(double(interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4)[x == c1, y == c2]) - (exp(c1+c2)));
			//error = abs(double(interp_func(pow(M, 0)* (row-1) + pow(M, 1)* (column-1), 4)[x == c1, y == c2]) - ((c1+c2)* sin(c1)));
			if(error > max_error)
			{
				max_error = error;
				max_c1 = c1; max_c2 = c2; 
				max_row = row; max_column = column;
			}
		}
	}

	cout << " & " << max_error ;
	//cout << "Max error: " << max_error << endl;
	//cout << "Max error at coordinates: " << max_c1 << "," << max_c2 << endl;
	//cout << "Coordinates fall in subregion" << max_row << "," << max_column << endl;
	//cout << "Interpolation polynomial in this subregion: " << interp_func(pow(M, 0)* (max_row-1) + pow(M, 1)* (max_column-1), 4) << " , Value: " << interp_func(pow(M, 0)* (max_row-1) + pow(M, 1)* (max_column-1), 4)[x == max_c1, y== max_c2] << endl;
	//cout << "Actual value: " << 2*max_c1 + 3*max_c2 +  max_c2*max_c2 + 1 << endl;
}


int main(void)
{
/*	Symbolic a("a"), b("b"), y, z;
	y = (Symbolic(7)^3); cout<< "y =" << y <<  endl;
	y = (a^3); cout << "y = " << y <<  endl;

	z = (a+b)*(a-b); cout << "z = " << z << endl;
	Symbolic A("A", 2,2);
	A(0,0) = 1.0; A(0, 1) = 0.89;
	A(1, 0) = 2.9; A(1, 1) = 9.0;

	cout << A;
	cout << "Determinant of A " << A.determinant() << endl;

	double points[3][2] = {{0, 0},{1, 1}, {1, 3}};

	double values[3] = {0, 3, 5} ;
	interpolation(3, 1, points, values);

	double **regular_points;
	regular_points = get_points_regular(3, 1, 2, 1, 1);
	int i;
	cout << "Points are :"<< endl;
	for(i =1; i <=3; i++)
	{
		cout << "(" << regular_points[i-1][0] << "," << regular_points[i-1][1] << ")" << endl; 
	}

*/

	/*Symbolic A("A", 3,3), x("x"), y("y"), det("det");
	A(0,0) = 1; A(0, 1) = 0.4; A(0, 2) = 0;
	A(1, 0) = 1; A(1, 1) = y; A(1, 2) = x;
	A(2, 0) = 1; A(2, 1) = 0.4; A(2, 2) = 0.2;
	det = A.determinant();
	cout << det << endl;*/

	/*double points[6][2] = {{0, 0},{1, 1}, {1, 3}, {4, 8}, {2, 9}, {3, 7}};

	double values[6] = {1, 7, 11, 45, 29, 33} ;
	interpolation_fast(6, 2, points, values);*/

	//For each M, we need to divide the domain into M^2 congruent subregions, then get points in each subregion
	// and then create  an interpolation
	int N, n;
	N = 1; n=0;
	int M, row, column, i, j, k;

	srand (static_cast <unsigned> (time(0))); // Seed the random number generator
	cout << "Piecewise Polynomial Interpolation Begins:" << endl;
	//Begin iterating for each M
	for(M=min_M; M <= max_M; M++)
	{
		cout<< M ;
		for(n=0, N =1; n < 4; n = n+1, N = N + n + 1)
		{
		Symbolic interp_func("interp_func", M*M, 9);
		//cout << "########M = " << M << "########" << endl;
		//Allocate memory to nodes 4d array which contains interpolation points in all subregions
		double ****nodes = (double ****)malloc(M * sizeof(double ***));
		double ***values = (double ***)malloc(M * sizeof(double **));
		for (i = 0; i < M; i++) 
		{
	        	nodes[i] = (double ***)malloc(M * sizeof(double**));
			values[i] = (double **)malloc(M * sizeof(double*));
			for(j =0 ; j< M; j++)
			{
				nodes[i][j] = (double **)malloc(N * sizeof(double *));
				values[i][j] = (double *)malloc(N * sizeof(double ));
				for(k=0; k< N; k++)
				{
					nodes[i][j][k] = (double *)malloc(2 * sizeof(double));
				}
			}
		}

		for(row = 1; row <= M; row++)
		{
			for(column=1; column <= M; column++)
			{
				//cout << "calculating interpolation for: Row: " << row << " Column: " << column << endl; 
				nodes[row-1][column-1] = get_points_regular(N, n, M, row, column);
				values[row-1][column-1] = get_values(N, nodes[row-1][column-1]);
				//cout << "Row: " << row << " Column: " << column << ": polynomial: "  ;
				interp_func(pow(M, 0)*(row-1) + pow(M, 1)*(column-1), 4) = interpolation_fast(N, n, nodes[row -1][column-1], values[row-1][column-1]);
				//cout << endl;	
			}
		}

		/*for(row = 1; row <= M; row++)
		{
			for(column=1; column <= M; column++)
			{

				//cout << "Row: " << row << " Column: " << column << ": polynomial: "  ;
				if(row != 1)
					interp_func(pow(M, 0)*(row-1) + pow(M, 1)*(column-1), 1) = weave_2(interp_func(pow(M, 0)*(row - 2) + pow(M, 1)*(column-1), 4), interp_func(pow(M, 0)*(row -1)+ pow(M, 1)* (column -1), 4), (row-1 - alpha)/M, (row-1+alpha)/M );
				else
					interp_func(pow(M, 0)*(row-1) + pow(M, 1)*(column-1), 1) = weave_2(interp_func(pow(M, 0)*(row - 1) + pow(M, 1)*(column-1), 4), interp_func(pow(M, 0)*(row -1)+ pow(M, 1)* (column -1), 4), (row-1 - alpha)/M, (row-1+alpha)/M );
				if(column != 1)
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 3) = weave_1(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-2), 4), interp_func(pow(M, 0)*(row -1)+ pow(M, 1)* (column -1), 4), (column-1 - alpha)/M, (column-1+alpha)/M );
				else
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 3) = weave_1(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-1), 4), interp_func(pow(M, 0)*(row -1)+ pow(M, 1)* (column -1), 4), (column-1 - alpha)/M, (column-1+alpha)/M );
				if(column != M)				
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 5) = weave_1(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-1), 4), interp_func(pow(M, 0)* (row -1)+ pow(M, 1)* (column), 4), (column - alpha)/M, (column+alpha)/M );
				else
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 5) = weave_1(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-1), 4), interp_func(pow(M, 0)* (row -1)+ pow(M, 1)* (column-1), 4), (column - alpha)/M, (column+alpha)/M );
				if(row != M)
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 7) = weave_2(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-1), 4), interp_func(pow(M, 0)* row + pow(M, 1) *(column -1), 4), (row - alpha)/M, (row+alpha)/M );
				else
					interp_func(pow(M, 0)*(row-1)+ pow(M, 1)* (column-1), 7) = weave_2(interp_func(pow(M, 0)*(row - 1)+ pow(M, 1)* (column-1), 4), interp_func(pow(M, 0)* (row-1) + pow(M, 1) *(column -1), 4), (row - alpha)/M, (row+alpha)/M );
				//cout << endl;	
			}
		}

		for(row = 1; row <= M; row++)
		{
			for(column=1; column <= M; column++)
			{

				//cout << "Row: " << row << " Column: " << column << ": polynomial: "  ;
				if(column != 1)
					interp_func(pow(M, 0)*(row-1)+ pow(M , 1)* (column-1), 0) = weave_1(interp_func(pow(M ,0)*(row - 1)+ pow(M, 1)* (column-2), 1), interp_func(pow(M ,0)* (row -1)+ pow(M, 1)* (column -1), 1), (column-1 - alpha)/M, (column-1+alpha)/M );
				else
					interp_func(pow(M, 0)*(row-1)+ pow(M , 1)* (column-1), 0) = weave_1(interp_func(pow(M ,0)*(row - 1)+ pow(M, 1)* (column-1), 1), interp_func(pow(M ,0)* (row -1)+ pow(M, 1)* (column -1), 1), (column-1 - alpha)/M, (column-1+alpha)/M );
				if(column != M)
					interp_func(pow(M ,0)* (row-1)+ pow(M, 1) *(column-1), 2) = weave_1(interp_func(pow(M ,0)* (row - 1)+ pow(M ,1)* (column-1), 1), interp_func(pow(M ,0) *(row -1) + pow(M ,1)*( column) , 1), (column - alpha)/M, (column + alpha)/M );
				else
					interp_func(pow(M ,0)* (row-1)+ pow(M, 1) *(column-1), 2) = weave_1(interp_func(pow(M ,0)* (row - 1)+ pow(M ,1)* (column-1), 1), interp_func(pow(M ,0) *(row -1) + pow(M ,1)*( column-1) , 1), (column - alpha)/M, (column + alpha)/M );
				if(column != 1)
					interp_func(pow(M, 0) * (row-1)+ pow(M , 1)* (column-1), 6) = weave_1(interp_func(pow(M , 0)* (row - 1)+ pow(M ,1)* (column-2), 7), interp_func(pow(M, 0)* (row -1)+ pow(M ,1)* (column-1), 7), (column-1 - alpha)/M, (column - 1 + alpha)/M );
				else
					interp_func(pow(M, 0) * (row-1)+ pow(M , 1)* (column-1), 6) = weave_1(interp_func(pow(M , 0)* (row - 1)+ pow(M ,1)* (column-1), 7), interp_func(pow(M, 0)* (row -1)+ pow(M ,1)* (column-1), 7), (column-1 - alpha)/M, (column - 1 + alpha)/M );
				if(column != M)
					interp_func(pow(M ,0)*(row-1)+ pow(M ,1)* (column-1), 8) = weave_2(interp_func(pow(M, 0)* (row - 1)+ pow(M ,1)* (column-1), 7), interp_func(pow(M, 0)* (row -1) + pow(M, 1)* (column) , 7), (column - alpha)/M, (column + alpha)/M );
				else
					interp_func(pow(M ,0)*(row-1)+ pow(M ,1)* (column-1), 8) = weave_2(interp_func(pow(M, 0)* (row - 1)+ pow(M ,1)* (column-1), 7), interp_func(pow(M, 0)* (row -1) + pow(M, 1)* (column-1) , 7), (column - alpha)/M, (column + alpha)/M );				
				//cout << endl;	
			}
		}*/

		//cout << interp_func << endl;
		//cout << "Now calculating error: " << endl;	
		calculate_error(interp_func, M);

		/*Symbolic w("w"), f("f"), g("g"), x("x"), y("y");
		cout << get_mollifier(x, -1, 1)<< endl;
		f = x;
		g = y;
		w = weave_1(f, g, -1, 1);
		cout << w << endl;*/
		}
		cout << endl << "\\\\" << endl;	
	}
		
	
}
