#include "makespl.h"
#include "piv_ge_solver.h"

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
double
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
}

double lg_n(int s,double x)
{
	if(s==0)
		return 1;
	if(s==1)
		return 1-x;
	return ((2*(s-1)+1-x)*lg_n(s-1,x)-(s-1)*lg_n(s-2,x))/s;
}

double fi(int i, double x)
{
	return lg_n(i,x);

	/*double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx[5];
	int		j;
	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];
	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 1 / h3 * (x - hx[0]) * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (h3 + 3 * h * h * (x - hx[1]) + 3 * h * (x - hx[1]) * (x - hx[1]) - 3 * (x - hx[1]) * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (h3 + 3 * h * h * (hx[3] - x) + 3 * h * (hx[3] - x) * (hx[3] - x) - 3 * (hx[3] - x) * (hx[3] - x) * (hx[3] - x));
	else
		return 1 / h3 * (hx[4] - x) * (hx[4] - x) * (hx[4] - x); */
}

/* Pierwsza pochodna fi */
double dfi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double dfi_po_dx = fi(i, x + dx) - fi(i, x - dx);
	dfi_po_dx /= 2 * dx;

	return dfi_po_dx;

	/* double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx[5];
	int		j;
	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];
	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 3 / h3 * (x - hx[0]) * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (3 * h * h + 6 * h * (x - hx[1]) - 9 * (x - hx[1]) * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (-3 * h * h - 6 * h * (hx[3] - x) + 9 * (hx[3] - x) * (hx[3] - x));
	else
		return -3 / h3 * (hx[4] - x) * (hx[4] - x);*/
}

/* double dfi(int i, double x)
{
	if (i == 0)
		return 0.0;
	else if (i == 1)
		return 1.0;
	else if (i == 2)
		return 4 * x;
	else {
		double dTn2 = 1;
		double dTn1 = 4 * x;
		double tmp;
		for (int j = 3; j <= i; j++) {
			double jj = j;
			tmp = dTn1;
			dTn1 = 2 * jj / (jj - 1) * x * dTn1 - jj / (jj - 2) * dTn2;
			dTn2 = tmp;
		}
		
		return dTn1;
	}
} */

/* Druga pochodna fi */
double d2fi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double d2fi_po_dx = dfi(a, b, n, i, x + dx) - dfi(a, b, n, i, x - dx);
	d2fi_po_dx /= 2 * dx;

	return d2fi_po_dx;

	/* double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx[5];
	int		j;
	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];
	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3 * (x - hx[0]);
	else if (x > hx[1] && x <= hx[2])
		return 1 / h3 * (6 * h - 18 * (x - hx[1]));
	else if (x > hx[2] && x <= hx[3])
		return 1 / h3 * (6 * h  -18 * (hx[3] - x));
	else
		return 6 / h3 * (hx[4] - x); */
}

/* double d2fi(int i, double x)
{
	if (i == 0 || i == 1)
		return 0.0;
	else if (i == 2)
		return 4.0;
	else if (i == 3)
		return 24.0 * x;
	else {
		double d2Tn2 = 4.0;
		double d2Tn1 = 24.0 * x;
		double tmp;
		for (int j = 4; j <= i; j++) {
			double jj = j;
			tmp = d2Tn1;
			d2Tn1 *= 2 * jj / (jj - 2) * x;
			d2Tn1 -= jj * jj / (jj - 2) / (jj - 2) * d2Tn2;
			d2Tn2 = tmp;
		}
		
		return d2Tn1;
	}
} */



/* Trzecia pochodna fi */
double d3fi(double a, double b, int n, int i, double x)
{
	double dx = (b - a) / n;
	
	double d3fi_po_dx = d2fi(a, b, n, i, x + dx) - d2fi(a, b, n, i, x - dx);
	d3fi_po_dx /= 2 * dx;

	return d3fi_po_dx;

	/* double		h = (b - a) / (n - 1);
	double		h3 = h * h * h;
	int		hi[5] = {i - 2, i - 1, i, i + 1, i + 2};
	double		hx[5];
	int		j;
	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];
	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return 6 / h3;
	else if (x > hx[1] && x <= hx[2])
		return -18 / h3;
	else if (x > hx[2] && x <= hx[3])
		return 18 / h3;
	else
		return -6 / h3; */
}


/* n-ta pochodna fi 
double dnfi(int n, int i, double x)
{
	if (i < n)
		return 0.0;
	else if (i == n)
		for (int j = 0; j < n; j++)
	else if (i == n + 1)
		return 24.0 * x;
	else {
		double d2Tn2 = 4.0;
		double d2Tn1 = 24.0 * x;
		double tmp;
		for (int j = 4; j <= i; j++) {
			double jj = j;
			tmp = d2Tn1;
			d2Tn1 *= 2 * jj / (jj - 2) * x;
			d2Tn1 -= jj * jj / (jj - 2) / (jj - 2) * d2Tn2;
			d2Tn2 = tmp;
		}
		
		return d2Tn1;
	}
} */

void draw_base(int i)
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
}

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
