#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double fun( double x ) {
				double r= ((double)rand() - RAND_MAX/2.0)/RAND_MAX/5; // +-10%
				return (1+r)*(sin(x));
}

int main( int argc, char **argv ) {
				FILE *out= argc > 1 ? fopen( argv[1], "w" ) : stdout;
				int n= argc > 2 ? atoi( argv[2] ) : 100;
				double a= argc > 3 ? atof( argv[3] ) : 0;
				double b= argc > 4 ? atof( argv[4] ) : 10;

				srand( argc > 5 ? atoi(argv[5]) : time(NULL) );

				int i;
				double dx = (b-a)/(n-1);

				for( i =0; i < n; i++ ) {
								fprintf( out, "%g %g\n", a+i*dx, fun(a+i*dx) );
				}
				fclose( out );

				return 0;
}
