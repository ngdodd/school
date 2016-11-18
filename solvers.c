/*
*				solvers.c
*
*	Be sure to use the -Ofast option to cc:
*		cc -Ofast -o solvers solvers.c -lm 
*
*	This program solves Laplace's equation on the unit square in
*	finite-difference form using Jacobi, Gauss-Seidel, or SOR iteration.
*
*	The BCs are Dirichlet with u(B) = 0 except u(x = 1, y) = 1.
*
*	In u(x_i,y_j) = u_ij = u[j][i], x[i] runs from xmin to xmax & y[j] runs
*	from ymin to ymax:
*		x = xmin + i*dx;
*		y = ymin + j*dy;
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#define MAX_SWEEPS	1000000

typedef void	*POINTER;	// pointer to an unknown data type 
typedef void	(*PFV)();	// pointer to a function returning void 

#define SOR_METHOD		1
#define GAUSS_SEIDEL_METHOD	2
#define JACOBI_METHOD		3

#define PI		3.14159265358979323846
#define EPSILON_M	2.2204460492503131e-16
#define sq(x)		( (x)*(x) )
#define max(a,b)	( ((a) > (b)) ? (a) : (b) )
#define min(a,b)	( ((a) < (b)) ? (a) : (b) )
#define ABS(x)		( ((x) >= 0) ? (x) : -(x) )
#define YES		(1)
#define NO		(0)
#define TRUE		(1)
#define FALSE		(0)
#define ERROR		(-1)
#define INT		((unsigned) sizeof(int))
#define FLOAT		((unsigned) sizeof(float))
#define DOUBLE		((unsigned) sizeof(double))
#define CHAR		((unsigned) sizeof(char))

typedef struct {
	int N;
	double dx, xmin, xmax;
	double dy, ymin, ymax;
	double h;
	double omega, EPSILON;
} GRID;

void implement_BCs(GRID *grid, double *u[]);
void print_solution(GRID *grid, double *u[]);
void SOR(GRID *grid, double *u[]);
void jacobi(GRID *grid, double *u[]);
POINTER Alloc(unsigned N_bytes);
POINTER alloc_vector(int N, unsigned element_size);
POINTER alloc_matrix(int N_rows, int N_columns, unsigned element_size);
void zero_matrix(double *A[], int N_rows, int N_columns);
void cpu_time(void);

int main(void)
{
	double **u;
	GRID grid;
	double xmin = 0., xmax = 1.;
	double ymin = 0., ymax = 1.;
	double h, mu;
	// for convergence of iterative methods, EPSILON = epsilon*h^2 
	double epsilon = 1.e-5;
	int N;
	int method = ERROR;
	PFV elliptic;
	char c[100];

	printf("Laplace equation solution\n");

	fprintf(stderr, "Choose iterative method.  Current choices are\n");
	fprintf(stderr, "\tJacobi (J),\n");
	fprintf(stderr, "\tGauss-Seidel (GS),\n");
	fprintf(stderr, "\tor SOR (SOR).\n");
	fprintf(stderr, "Enter choice: ");
	scanf("%s", c);
	if (strcmp(c, "J") == 0 || strcmp(c, "j") == 0) {
		method = JACOBI_METHOD;
		printf("iterative method = Jacobi\n");
	}
	else if (strcmp(c, "GS") == 0 || strcmp(c, "gs") == 0) {
		method = GAUSS_SEIDEL_METHOD;
		printf("iterative method = Gauss-Seidel\n");
	}
	else if (strcmp(c, "SOR") == 0 || strcmp(c, "sor") == 0) {
		method = SOR_METHOD;
		printf("iterative method = SOR\n");
	}
	else {
		fprintf(stderr, "ERROR: unrecognized method %s\n", c);
		exit(ERROR);
	}

	fprintf(stderr, "Enter the number of dx (= number of dy): ");
	scanf("%d", &N);
	printf("number of dx = number of dy = %d\n", N);
	grid.N = N;
	grid.dx = (xmax-xmin)/N;
	grid.dy = (ymax-ymin)/N;
	h = grid.h = grid.dx;
	grid.xmin = xmin;
	grid.xmax = xmax;
	grid.ymin = ymin;
	grid.ymax = ymax;
	grid.EPSILON = epsilon*sq(h);
	printf("EPSILON = epsilon*h^2 = %g\n", grid.EPSILON);
	printf("\twhere epsilon = %g\n", epsilon);

	switch (method) {
	case SOR_METHOD:
		elliptic = SOR;
		// calculate omega_opt for SOR 
		mu = cos(PI*h); // Jacobi spectral radius 
		grid.omega = 2.*(1.-sqrt(1.-sq(mu)))/sq(mu);
		printf("SOR omega_opt = %g\n", grid.omega);
		break;
	case GAUSS_SEIDEL_METHOD:
		elliptic = SOR;
		grid.omega = 1.;
		break;
	case JACOBI_METHOD:
		elliptic = jacobi;
		break;
	default:
		printf("ERROR: unrecognized method %d\n", method);
		exit(ERROR);
		break;
	}

	// allocate memory for array u 
	u = (double **) alloc_matrix(N+1, N+1, DOUBLE);
	zero_matrix(u, N+1, N+1);

	implement_BCs(&grid, u);
	(*elliptic)(&grid, u);
	print_solution(&grid, u);
	cpu_time();

	return 0;
}

void implement_BCs(GRID *grid, double *u[])
{
	int N = grid->N, i, j;

	for (j = 0; j <= N; j++) {
		u[j][0] = 0.;
		u[j][N] = 1.;
	}
	for (i = 1; i < N; i++) {
		u[0][i] = 0.;
		u[N][i] = 0.;
	}
}

void print_solution(GRID *grid, double *u[])
{
	int N = grid->N, i, j;
	double x, y;
	FILE *file1;

	file1 = fopen("laplace.sol", "w");

	// x, y, u
	for (i = 0; i <= N; i++) {
		for (j = 0; j <= N; j++) {
			x = grid->xmin + i*grid->dx;
			y = grid->ymin + j*grid->dy;
			fprintf(file1, "%-12g %-12g %-12g\n", x, y, u[j][i]);
		}
	}

	fclose(file1);
}

void SOR(GRID *grid, double *u[])
{
	int N = grid->N, i, j, sweep;
	double EPSILON = grid->EPSILON;
	double omega = grid->omega;
	double sum, residual, norm_residual0, norm_residual = 1.;

	sum = 0;
	for (i = 1; i < N; i++) {
		for (j = 1; j < N; j++) {
			residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
				u[j][i+1]+u[j][i-1]);
			sum += fabs(residual);
		}
	}
	norm_residual0 = sum/sq(N-1);

	for (sweep = 0; sweep<MAX_SWEEPS && 
			norm_residual>EPSILON*norm_residual0; sweep++) {
		sum = 0;
		for (i = 1; i < N; i++) {
			for (j = 1; j < N; j++) {
				residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
					u[j][i+1]+u[j][i-1]);
				sum += fabs(residual);
				u[j][i] += omega*residual/4.;
			}
		}
		norm_residual = sum/sq(N-1);
	}
	printf("number of sweeps = %d\n", sweep);
	printf("norm_residual/norm_residual0 = %g\n", 
		norm_residual/norm_residual0);
	if (norm_residual > EPSILON*norm_residual0)
		printf("WARNING: iterative method failed to converge\n");
}

void jacobi(GRID *grid, double *u[])
{
	static double **new_u;
	int N = grid->N, i, j, sweep;
	double EPSILON = grid->EPSILON;
	double sum, residual, norm_residual0, norm_residual = 1.;

	new_u = (double **) alloc_matrix(N+1, N+1, DOUBLE);
	for (i = 0; i <= N; i++)
		for (j = 0; j <= N; j++) new_u[j][i] = u[j][i];
	sum = 0;
	for (i = 1; i < N; i++) {
		for (j = 1; j < N; j++) {
			residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
				u[j][i+1]+u[j][i-1]);
			sum += fabs(residual);
		}
	}
	norm_residual0 = sum/sq(N-1);

	for (sweep = 0; sweep<MAX_SWEEPS && 
			norm_residual>EPSILON*norm_residual0; sweep++) {
		sum = 0;
		for (i = 1; i < N; i++) {
			for (j = 1; j < N; j++) {
				residual = (-4.*u[j][i]+u[j+1][i]+u[j-1][i]+
					u[j][i+1]+u[j][i-1]);
				sum += fabs(residual);
				new_u[j][i] += residual/4.;
			}
		}
		norm_residual = sum/sq(N-1);
		for (i = 1; i < N; i++)
			for (j = 1; j < N; j++) u[j][i] = new_u[j][i];
	}

	printf("number of sweeps = %d\n", sweep);
	printf("norm_residual/norm_residual0 = %g\n", 
		norm_residual/norm_residual0);
	if (norm_residual > EPSILON*norm_residual0)
		printf("WARNING: iterative method failed to converge\n");
}

POINTER alloc_matrix(int N_rows, int N_columns, unsigned element_size)
{
	int i, space, N_pointers;
	POINTER pA, array_origin;

		// due to alignment requirement 
	N_pointers = (N_rows%2) ? N_rows+1 : N_rows; 

	space = N_pointers*sizeof(POINTER) + N_rows*N_columns*element_size;
	pA = Alloc((unsigned) space);

	array_origin = ((char *) (pA)) + N_pointers*sizeof(POINTER);

	for (i = 0; i < N_rows; i++) 
		((POINTER *) pA)[i] = ((char *) array_origin) +
			i*N_columns*element_size;

	return pA;
}

POINTER alloc_vector(int N, unsigned element_size)
{
	return Alloc((unsigned) (N*element_size));
}

POINTER Alloc(unsigned N_bytes)
{
	POINTER p;

	p = malloc(N_bytes);

	if (p == (POINTER) NULL) {
		perror("ERROR in Alloc(): malloc() returned NULL pointer");
		exit(ERROR);
	}
	return p;
}

void zero_matrix(double *A[], int N_rows, int N_columns)
{
	int i, j;

	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++) A[i][j] = 0.;
}

void cpu_time(void)
{
	struct rusage rusage_buffer;
	long user_sec, user_usec, user_min, sys_sec, sys_usec, sys_min;
	long cpu_sec, cpu_usec, cpu_min;
	long MILLION = 1000000;

	getrusage(0, &rusage_buffer);

	user_sec = rusage_buffer.ru_utime.tv_sec;
	user_usec = rusage_buffer.ru_utime.tv_usec;
	sys_sec = rusage_buffer.ru_stime.tv_sec;
	sys_usec = rusage_buffer.ru_stime.tv_usec;

	cpu_usec = (user_usec + sys_usec)%MILLION;
	cpu_sec = user_sec + sys_sec + (user_usec + sys_usec)/MILLION;
	cpu_min = cpu_sec/60;
	cpu_sec = cpu_sec%60;
	user_min = user_sec/60;
	user_sec = user_sec%60;
	sys_min = sys_sec/60;
	sys_sec = sys_sec%60;

	printf("\nuser cpu time: ");
	if (user_min) printf("%ld min ", user_min);
	printf("%g sec\n", user_sec + user_usec/(1000000.));

	printf("system time: ");
	if (sys_min) printf("%ld min ", sys_min);
	printf("%g sec\n", sys_sec + sys_usec/(1000000.));

	printf("total cpu time: ");
	if (cpu_min) printf("%ld min ", cpu_min);
	printf("%g sec\n", cpu_sec + cpu_usec/(1000000.));
}
