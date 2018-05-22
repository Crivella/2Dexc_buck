#Library env variables

PREFIX=$(shell pwd)/../..

INC =./inc
SDIR=./src
ODIR=./obj
CURDIR=$(shell pwd)
LDIR=$(PREFIX)/lib
IDIR=$(PREFIX)/include

LD_FLAGS=

ifndef CC
  CC=gcc
endif
ifndef MPICC
  MPICC=mpicc
endif
CFLAGS=-g -O2 -Wall -fpic -std=c11 -I$(IDIR) -I$(INC) $(LD_FLAGS)
