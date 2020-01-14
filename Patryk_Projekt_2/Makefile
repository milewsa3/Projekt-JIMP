poly: main.o splines.o points.o polimator.o gaus/libge.a
	$(CC) -o poly main.o splines.o points.o polimator.o gaus/libge.a -L gaus -l ge

aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o	

polimator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c polimator.c
aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c


rysuj:
	./poly -s spl/$(FUN).spl -p test/$(FUN).gen -g wykresy/$(FUN).wyk -f 0 -t 10 -n 300

baza:
	export APPROX_BASE_SIZE=$(BAZA)

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta poly
