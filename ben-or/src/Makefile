all:
	gcc -Wall -g -pg -O3  -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp ben-or.c

amd:
	gcc -Wall -g -pg -O3 -mtune=native -march=znver2 -ffast-math -funroll-loops  -fopenmp ben-or.c

clang:
	clang -Wall -g -pg -Ofast -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp ben-or.c

clean:
	rm -f a.out

