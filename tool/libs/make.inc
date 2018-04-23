#Library env variables

PREFIX=$(shell pwd)/../..

INC =./inc
SDIR=./src
ODIR=./obj
CURDIR=$(shell pwd)
LDIR=$(PREFIX)/lib
IDIR=$(PREFIX)/include

LD_FLAGS=

CC=gcc
MPICC=mpicc
CFLAGS=-g -O2 -Wall -fpic -std=c11 -I$(IDIR) -I$(INC) $(LD_FLAGS)
