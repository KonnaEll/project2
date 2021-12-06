search: search.o lsh_funcs.o lsh_for_vectors.o lsh_for_frechet.o cube_funcs.o cube_for_vectors.o
	gcc -Wall -Werror -std=c9x search.o lsh_funcs.o lsh_for_vectors.o lsh_for_frechet.o cube_funcs.o cube_for_vectors.o -o search -lm

search.o: search.c
	gcc -Wall -Werror -std=c9x -c search.c

lsh_funcs.o: lsh_funcs.c
	gcc -Wall -Werror -std=c9x -c lsh_funcs.c

lsh_for_vectors.o: lsh_for_vectors.c
	gcc -Wall -Werror -std=c9x -c lsh_for_vectors.c

lsh_for_frechet.o: lsh_for_frechet.c
	gcc -Wall -Werror -std=c9x -c lsh_for_frechet.c

cube_funcs.o: cube_funcs.c
	gcc -Wall -Werror -std=c9x -c cube_funcs.c

cube_for_vectors.o: cube_for_vectors.c
	gcc -Wall -Werror -std=c9x -c cube_for_vectors.c


clean:
	rm -f search *.o lsh_for_vectors *.o search
