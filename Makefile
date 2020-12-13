all: snf saturation_vectors

snf: ./src/gap_snf.c
	gcc -fopenmp ./src/gap_snf.c ./src/snf.c -o ./src/gap_snf
saturation_vectors: ./src/gap_saturation_vectors.c
	gcc -fopenmp ./src/gap_saturation_vectors.c ./src/snf.c -o ./src/gap_saturation_vectors
