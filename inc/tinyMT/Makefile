OBJECTS = tinymt32.o 
#tinymt64.o

CFLAGS = -O3 --std=c11

libTinyMT.a : ${OBJECTS}
	ar cr $@ $<

clean :
	rm ${OBJECTS}
