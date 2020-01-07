
SOURCES= dyn_splay.c dyn_tree.c

OBJS= $(CFILES1:%.c=%.o)

CC= gcc

CFLAGS=-O -c -D"DOUBLE int"  

dyn_tree: ${OBJS} 
	${CC} ${CFLAGS} ${SOURCES}


dyn_splay.o: _dyn_tree.h dyn_tree.h
dyn_tree.o:  _dyn_tree.h dyn_tree.h

