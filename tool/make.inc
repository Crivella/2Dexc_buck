#Env variables

PREFIX=$(HOME)/.local

CURDIR=$(shell pwd)
LDIR=$(PREFIX)/lib
IDIR=$(PREFIX)/include
BDIR=$(PREFIX)/bin

ODIR=./obj

LD_FLAGS=

CC=gcc
MPICC=mpicc
CFLAGS=-g -O2 -Wall -rdynamic -fpic -std=c11 -I$(IDIR) $(LD_FLAGS)
LD=-L$(LDIR) -L/user/local/lib -lqepp -lcterm -lxml2 -lgsl -lgslcblas -lm -ldl

export LD_LIBRARY_PATH=$(LDIR)
