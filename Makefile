# Audio-Resampler Makefile

CC := gcc

utils := art artest

all: $(utils)

art: art.c stretch.c stretch.h resampler.c resampler.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS art.c stretch.c resampler.c decimator.c workers.c biquad.c -lm -pthread -o art

artest: artest.c resampler.c resampler.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS artest.c resampler.c decimator.c workers.c biquad.c -lm -pthread -o artest

clean:
	rm -f art artest
