INCLUDEPATH=~/software/ATLAS/include
LIBPATH=~/software/ATLAS/lib/AC
LIBS=  -lblas -lcblas -latlas -llapack -lcblas  -latlas -lm
FLAGS= -g 
CC=gcc

all: extract ar-gmm_convert ar-gmm


ar-gmm_convert:ar-gmm_convert.c
	${CC} ${FLAGS} -I ${INCLUDEPATH} -L ${LIBPATH} ${LIBS} ar-gmm_convert.c -o ar-gmm_convert
ar-gmm_train:ar-gmm.c
	${CC} ${FLAGS} ar-gmm.c ${LIBS}  -o ar-gmm  -I ${INCLUDEPATH} -L ${LIBPATH}

extract:extract.c
	gcc ${FLAGS}  extract.c -o extract

