CC = gcc
CFLAG = -g
#CFLAG = -O2

all: DBUX #default setting when use solely "make"

DBUX: main.o functions.o define.h #setting when use "make DBUX"
	$(CC) -o DBUX.exe main.o functions.o

clean: #setting when use "make clean"
	rm -f *.o
