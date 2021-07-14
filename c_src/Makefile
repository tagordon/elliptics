all: ellip.so ellip.so

ellip.so: ellip.o 
	gcc -shared -o ellip.so ellip.o

ellip.o: ellip.c
	gcc -c -Wall -Werror -fpic ellip.c
