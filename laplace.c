/*
*				laplace.c
*
*	This program solves Laplace's equation on the unit square in
*	finite-difference form using Gaussian elimination (solve() from 
*	solve.c).  
*
*	The BCs are Dirichlet with u(B) = 0 except u(x = 1, y) = 1.
*
*	An "exact" Fourier series solution is calculated in fourier_solution()
*	if DEBUG = YES.
*
*	In u(x_i,y_j) = u_ij = u[j][i], x[i] runs from xmin to xmax & y[j] runs
*	from ymin to ymax:
*		x = xmin + i*dx;
*		y = ymin + j*dy;
*
*	For A & b, a single index 
*		n = j*(N-1) + i 
*	is constructed from (i,j), thus reading interior solution values u_ij 
*	from left to right & bottom to top.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

typedef void	*POINTER;	// pointer to an unknown data type 

typedef struct {
	int N;
	double dx, xmin, xmax;
	double dy, ymin, ymax;
} GRID;

#define DEBUG	YES

void solve(int N, double *A[], double b[]);
void implement_BCs(GRID *grid, double *u[]);
void print_solution(GRID *grid, double *u[]);
POINTER Alloc(unsigned N_bytes);
POINTER alloc_vector(int N, unsigned element_size);
POINTER alloc_matrix(int N_rows, int N_columns, unsigned element_size);
void zero_vector(double v[], int N);
void zero_matrix(double *A[], int N_rows, int N_columns);
double fourier_solution(double x, double y);
void print_matrix(char *string, double *A[], int N_rows, int N_columns);
void print_vector(char *string, double v[], int N);

int main(void)
{
	double **A, **u, *b;
	GRID grid;
	double xmin = 0., xmax = 1.;
	double ymin = 0., ymax = 1.;
	int i, j, n, N, modes;

	printf("Direct solver for Laplace Equation\n");

	fprintf(stderr, "Enter the number of dx (= number of dy): ");
	scanf("%d", &N);
	printf("number of dx = number of dy = %d\n", N);
	grid.N = N;
	grid.dx = (xmax-xmin)/N;
	grid.dy = (ymax-ymin)/N;
	grid.xmin = xmin;
	grid.xmax = xmax;
	grid.ymin = ymin;
	grid.ymax = ymax;

	// allocate memory for arrays 
	modes = sq(N-1); // unknowns = interior u_ij values 
	A = (double **) alloc_matrix(modes, modes, DOUBLE);
	b = (double *) alloc_vector(modes, DOUBLE);
	u = (double **) alloc_matrix(N+1, N+1, DOUBLE);

	zero_matrix(A, modes, modes);
	zero_vector(b, modes);
	zero_matrix(u, N+1, N+1);

	implement_BCs(&grid, u);

	//	set up A & b in Laplace's equation
	//		A u = b

	for (i = 0; i <= N-2; i++) {
		for (j = 0; j <= N-2; j++) {
			n = j*(N-1) + i;
			// Note that n = 0, 1, ..., N(N-2) = (N-1)^2 - 1
			A[n][n] = -4.;
			// BEWARE: there are some zeros in sub- & superdiagonals
			if (i > 0) A[n][n-1] = 1.;       // subdiagonal
			if (i < N-2) A[n][n+1] = 1.;     // superdiagonal
			if (j > 0) A[n][n-(N-1)] = 1.;   // sub-fringe
			if (j < N-2) A[n][n+(N-1)] = 1.; // super-fringe
			// or equivalently for the fringes (all ones):
			// if (n >= N-1) A[n][n-(N-1)] = 1.;
			// if (n <= sq(N-1)-N) A[n][n+(N-1)] = 1.;
		}
	}
	for (j = 0; j <= N-2; j++) {
		n = j*(N-1);
		b[n] -= u[j+1][0]; // West BC 
		n = j*(N-1) + (N-2);
		b[n] -= u[j+1][N]; // East BC 
	}
	for (i = 0; i <= N-2; i++) {
		n = i;
		b[n] -= u[0][i+1]; // South BC 
		n = (N-2)*(N-1) + i;
		b[n] -= u[N][i+1]; // North BC 
	}
	if (DEBUG && modes <= 25) {
		print_matrix("matrix A", A, modes, modes);
		print_vector("vector b", b, modes);
	}

	solve(modes, A, b);
	for (i = 0; i <= N-2; i++) {
		for (j = 0; j <= N-2; j++) {
			n = j*(N-1) + i;
			u[j+1][i+1] = b[n];
		}
	}
	print_solution(&grid, u);

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

	if (DEBUG) {
		// x, y, u, fourier u
		for (i = 0; i <= N; i++)
		for (j = 0; j <= N; j++) {
			x = grid->xmin + i*grid->dx;
			y = grid->ymin + j*grid->dy;
			fprintf(file1, "%-12g %-12g %-12g %-12g\n",
				x, y, u[j][i],
				((i == N) ? 1. : fourier_solution(x, y)));
		}
	}
	else {
		// x, y, u
		for (i = 0; i <= N; i++)
		for (j = 0; j <= N; j++) {
			x = grid->xmin + i*grid->dx;
			y = grid->ymin + j*grid->dy;
			fprintf(file1, "%-12g %-12g %-12g\n", x, y, u[j][i]);
		}
	}

	fclose(file1);
}

void solve(int N, double *A[], double b[])
{
	int i, j, k, pivot_row;
	double pivot, multiplier, *rowswap, swap, sum;
	double EPSILON = N*pow(2., -52.); // lower bound for pivots 

		// transform A to upper triangular form 
	// loop over columns of A 
	for (k = 0; k < N-1; k++) {
		pivot = A[k][k];
		pivot_row = k;
		// loop over rows of A to find pivot 
		for (i = k+1; i < N; i++) {
			if (ABS(pivot) < ABS(A[i][k])) {
				pivot = A[i][k];
				pivot_row = i;
			}
		}
		if (ABS(pivot) < EPSILON) {
			printf("ERROR in solve(): pivot too small\n");
			exit(ERROR);
		}
		if (pivot_row != k) {
			// swap rows k & pivot_row of A 
			rowswap = A[k];
			A[k] = A[pivot_row];
			A[pivot_row] = rowswap;
			// swap rows of b 
			swap = b[k];
			b[k] = b[pivot_row];
			b[pivot_row] = swap;
		}
		// loop over rows of A 
		for (i = k+1; i < N; i++) {
			// compute multiplier 
			multiplier = A[i][k] = A[i][k]/A[k][k];
			// compute new row elements of A 
			for (j = k+1; j < N; j++)
				A[i][j] -= multiplier*A[k][j];
			// compute new elements of b 
			b[i] -= multiplier*b[k];
		}
	}

		// back substitution 
	// loop over rows in reverse order 
	for (i = N-1; i >= 0; i--) {
		sum = 0.;
		for (j = i+1; j < N; j++)
			sum += A[i][j]*b[j];
		sum = b[i] - sum;
		if (ABS(A[i][i]) > EPSILON) b[i] = sum/A[i][i];
		else if (ABS(sum) < EPSILON) {
			b[i] = 0.;
			printf("WARNING: sum may be too small in solve()\n");
		}
		else {
			printf("ERROR in solve()\n");
			exit(ERROR);
		}
	}
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

void zero_vector(double v[], int N)
{
	int i;

	for (i = 0; i < N; i++) v[i] = 0.;
}

void zero_matrix(double *A[], int N_rows, int N_columns)
{
	int i, j;

	for (i = 0; i < N_rows; i++)
		for (j = 0; j < N_columns; j++) A[i][j] = 0.;
}

//			fourier_solution()
//
//	Calculates the solution to Laplace's equation based on the Fourier 
//	series expansion.

// number of Fourier modes in fourier_solution() 
#define N_fourier	64

double fourier_solution(double x, double y)
{
	static double f[N_fourier];
	double sum;
	int n;
	static int first = TRUE;

	if (first) {
		first = FALSE;
		for (n = 0; n < N_fourier; n++)
			f[n] = 4./((2.*n+1.)*PI);
	}

	sum = 0.;
	for (n = 0; n < N_fourier; n++)
		sum += f[n]*sin((2.*n+1.)*PI*y)*
			(sinh((2.*n+1.)*PI*x)/sinh((2.*n+1.)*PI));

	return sum;
}

void print_matrix(char *string, double *A[], int N_rows, int N_columns)
{
	int i, j;

	printf("\n%s\n", string);
	for (i = 0; i < N_rows; i++) {
		for (j = 0; j < N_columns; j++) 
			printf("%3d", (int) A[i][j]);
		printf("\n");
	}
}

void print_vector(char *string, double v[], int N)
{
	int i;

	printf("\n%s\n", string);
	for (i = 0; i < N; i++) printf("%3d", (int) v[i]);
	printf("\n");
}
