cluster: cluster.o cluster_funcs.o lsh_funcs.o cube_funcs.o tree.o lsh_for_frechet.o
	gcc -Wall -Werror -std=c9x cluster.o cluster_funcs.o lsh_funcs.o cube_funcs.o tree.o lsh_for_frechet.o -o cluster -lm

cluster.o: cluster.c
	gcc -Wall -Werror -std=c9x -c cluster.c

cluster_funcs.o: cluster_funcs.c
	gcc -Wall -Werror -std=c9x -c cluster_funcs.c

lsh_funcs.o: lsh_funcs.c
	gcc -Wall -Werror -std=c9x -c lsh_funcs.c

cube_funcs.o: cube_funcs.c
	gcc -Wall -Werror -std=c9x -c cube_funcs.c

tree.o: tree.c
	gcc -Wall -Werror -std=c9x -c tree.c

lsh_for_frechet.o: lsh_for_frechet.c
	gcc -Wall -Werror -std=c9x -c lsh_for_frechet.c

clean:
	rm -f cube*.o lsh*.o cluster*.o tree*.o cluster
