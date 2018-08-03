#Adaptation of Triangle makefile to compile only the triangle.c file, 
#into a standalone library that can then be used to compile mex modules

include ./configure.make

AR = ar

TRILIBDEFS = -DTRILIBRARY

all: triangle.$(LIB_EXT)

OBJECTS=triangle.$(OBJ_EXT)

triangle.$(LIB_EXT): $(OBJECTS)
	lib -out:triangle.$(LIB_EXT) $(OBJECTS)

triangle.$(OBJ_EXT): triangle.c triangle.h
	$(CC) $(CSWITCHES) $(TRILIBDEFS) -c triangle.c

clean: 
	rm -rf *.$(LIB_EXT) *.$(OBJ_EXT) *.LIB *.LST *.$(OBJ_EXT)bj *.BAK
