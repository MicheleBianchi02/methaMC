all: spectrum libspectrum.so

# Compile the standalone CLI executable
spectrum: spectrum.c pcg_basic.c
	$(CC) spectrum.c pcg_basic.c -Wall -W -O3 -march=native -ffast-math -flto -g -lm -o spectrum

# Compile the shared library for Python
libspectrum.so: spectrum.c pcg_basic.c
	$(CC) spectrum.c pcg_basic.c -Wall -W -O3 -march=native -ffast-math -flto -g -shared -fPIC -DBUILD_AS_LIB -lm -o libspectrum.so

clean:
	rm -f spectrum libspectrum.so
