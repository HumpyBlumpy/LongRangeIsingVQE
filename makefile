default : 2local-qaoa.c 2local-qaoa.h
	gcc -g -Wall -shared -o lib2local-qaoa.so -fPIC 2local-qaoa.c
