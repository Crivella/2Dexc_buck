#Env variables

PREFIX=$(shell pwd)/..

CURDIR=$(shell pwd)
LDIR=$(PREFIX)/lib
IDIR=$(PREFIX)/include
BDIR=$(PREFIX)/bin

ODIR=./obj

LD_FLAGS=

CC=gcc
MPICC=mpiicc
CFLAGS=-g -O2 -Wall -rdynamic -fpic -std=c11 -I$(IDIR) $(LD_FLAGS)
LD=-Wl,-rpath=$(LDIR) -L$(LDIR) -L/user/local/lib -lqepp -lm
