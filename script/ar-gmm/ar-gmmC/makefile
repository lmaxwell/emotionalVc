INCLUDEPATH=~/software/ATLAS/include
LIBPATH=~/software/ATLAS/lib/AC
LIBS=  -lblas -lcblas -latlas -llapack -lcblas  -latlas -lm
FLAGS= -g -O2
CC=gcc

all: extract ar-gmm_convert ar-gmm


ar-gmm_convert:ar-gmm_convert.c
	${CC} ${FLAGS} -I ${INCLUDEPATH} -L ${LIBPATH} ${LIBS} ar-gmm_convert.c -o ar-gmm_convert
ar-gmm_train:ar-gmm.c
	${CC} ${FLAGS} ar-gmm.c ${LIBS}  -o ar-gmm  -I ${INCLUDEPATH} -L ${LIBPATH}
gmm_init:gmm_init.c
	${CC} ${FLAGS} -I . -I /home/lixian/speech/SPTK-3.6/include -L /home/lixian/speech/SPTK-3.6/lib -lm gmm_init.c -lSPTK  -o gmm_init
extract:extract.c
	gcc ${FLAGS}  extract.c -o extract

