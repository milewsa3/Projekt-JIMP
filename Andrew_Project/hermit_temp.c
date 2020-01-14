#include "makespl.h"
#include "gaus/piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>


void  make_spl ( points_t *pts, spline_t *spl)
{
	double suma_x;
	double suma_x_kw;
	double suma_y;
	double suma_x_y;
	matrix_t *eqs;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  	char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	
	eqs = make_matrix(2,3);

	for (int i=0;i<spl->n;i++)
	{
		suma_x += spl->x[i];
		suma_y += spl->f[i];
		suma_x_y += spl->x[i] * spl->f[i];
		suma_x_kw += spl->x[i] * spl -> x[i];
	}

	put_entry_matrix(eqs,0,0,spl->n);
	put_entry_matrix(eqs,1,2,suma_x);
	put_entry_matrix(eqs,1,3,suma_x);
	put_entry_matrix(eqs,2,1,suma_x);
	put_entry_matrix(eqs,2,2,suma_x_kw);
	put_entry_matrix(eqs,2,3,suma_x_y);
	
	#ifdef DEBUG
		write_matrix( eqs, stdout );
		printf("Punktow w danych: %d\n",pts->n);
		for (int i=0;i<pts->n;i++)
		{
			printf("X: %lf Y: %lf\n",pts -> x[i],pts -> y[i]);
		}
	#endif
	
		if( piv_ge_solver( eqs ) ) {
			spl->n = 0;
			return;
		}

	#ifdef DEBUG
		write_matrix( eqs, stdout );
	#endif
}


