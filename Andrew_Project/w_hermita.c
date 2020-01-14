#include "makespl.h"
#include "gaus/piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>


#define DX_MULTIPLY 10


/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe: n - liczba funkcji a,b - granice przedzialu aproksymacji i
 * - numer funkcji x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */


/*double
value_spl (spline_t * spl, double x)
{
  int i;
  double dx;

  for (i = spl->n - 1; i > 0; i--)
    if (spl->x[i] < x)
      break;

  dx = x - spl->x[i];

  return spl->f[i]
	+ dx * spl->f1[i]
	+ dx * dx / 2 *  spl->f2[i]
	+ dx * dx * dx / 6 * spl->f3[i];
}*/

double lg_n(int s,double x)
{
	if(s==0)
		return 1;
	else if(s==1)
		return 2*x;
	else if (s==2)
		return 4*x*x - 2;
	else if (s==3)
		return 8*x*x*x - 12*x;
	else if (s==4)
		return 16*x*x*x*x - 48 *x*x +12;
	else if (s==5)
		return 32*x*x*x*x*x - 160*x*x*x + 120*x;
	else if (s==6)
		return 64*x*x*x*x*x*x -480*x*x*x*x + 720*x*x -120;
	else if (s==7)
		return 128*x*x*x*x*x*x*x -1344*x*x*x*x*x + 3360*x*x*x -1680*x;
	else 
		return 1;
}

double fi(int i, double x)
{
	return lg_n(i,x);
}

/* Pierwsza pochodna fi */
double dfi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double dfi_po_dx = fi(i, x + dx) - fi(i, x - dx);
	dfi_po_dx /= 2 * dx;

	return dfi_po_dx;
}

double d2fi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double d2fi_po_dx = dfi(a, b, n, i, x + dx) - dfi(a, b, n, i, x - dx);
	d2fi_po_dx /= 2 * dx;

	return d2fi_po_dx;

}




/* Trzecia pochodna fi */
double d3fi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double d3fi_po_dx = d2fi(a, b, n, i, x + dx) - d2fi(a, b, n, i, x - dx);
	d3fi_po_dx /= 2 * dx;

	return d3fi_po_dx;
}


/*void draw_base(int i)
{
	FILE *out = fopen("debug_base", "w");
	if (out == NULL) {
		fprintf(stderr, "Nie dziala!\n");
		return;
	}
	double x = 0.0;
	double dx = 0.01;

	for (int j = 0; j <= i; j++)
		for (int k = 0; k < 101; k++)
			fprintf(out, "%lf %lf\n", x + k * dx, fi(j, x + k * dx));

	fclose(out);
}*/

void make_spl(points_t * pts, spline_t * spl)
{
	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
	char		*nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );


	eqs = make_matrix(nb, nb + 1);

// Filling matrix

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi(j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(j, x[k]));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	// Creating and filling splines struct
	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0 * DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double		ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (k, xx);
				spl->f1[i] += ck * dfi (a, b, DX_MULTIPLY * pts->n, k, xx);
				spl->f2[i] += ck * d2fi(a, b, DX_MULTIPLY * pts->n, k, xx);
				spl->f3[i] += ck * d3fi(a, b, DX_MULTIPLY * pts->n, k, xx);
			}
		}
	}

	free_matrix(eqs);
}
