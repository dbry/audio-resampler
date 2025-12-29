# Audio-Resampler Makefile

CC := gcc

utils := art art64 artest artest64

all: $(utils)

art: art.c stretch.c stretch.h resampler.c resampler.h extrapolator.c extrapolator.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -Wall -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS -DENABLE_EXTRAPOLATION art.c stretch.c resampler.c extrapolator.c decimator.c workers.c biquad.c -lm -pthread -o art

art64: art.c stretch.c stretch.h resampler.c resampler.h extrapolator.c extrapolator.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -DPATH_WIDTH=64 -Wall -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS -DENABLE_EXTRAPOLATION art.c stretch.c resampler.c extrapolator.c decimator.c workers.c biquad.c -lm -pthread -o art64

artest: artest.c resampler.c resampler.h extrapolator.c extrapolator.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -Wall -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS -DENABLE_EXTRAPOLATION artest.c resampler.c extrapolator.c decimator.c workers.c biquad.c -lm -pthread -o artest

artest64: artest.c resampler.c resampler.h extrapolator.c extrapolator.h decimator.c decimator.h workers.c workers.h biquad.c biquad.h
	$(CC) -DPATH_WIDTH=64 -Wall -O3 -mavx2 -fno-signed-zeros -fno-trapping-math -fassociative-math -DENABLE_THREADS -DENABLE_EXTRAPOLATION artest.c resampler.c extrapolator.c decimator.c workers.c biquad.c -lm -pthread -o artest64

clean:
	rm -f art art64 artest artest64
