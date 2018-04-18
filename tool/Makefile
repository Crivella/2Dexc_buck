#Top makefile

include make.inc

DIRS=libs src

.PHONY: libs src clean veryclean install uninstall

default : all install

libs		: 
	@if test -d libs; then \
	(cd libs ; $(MAKE) TLDEPS= all || exit 1) ; fi
	@echo "--------------------------------------------------------------------------------\n\n\n\n"

qepp		:
	@if test -d libs; then \
	(cd libs ; $(MAKE) TLDEPS= $@ || exit 1) ; fi
	@echo "--------------------------------------------------------------------------------"

src		: libs install
	@if test -d src; then \
	(cd src ; cp ../make.inc . ; $(MAKE) TLDEPS= all || exit 1) ; fi
	@echo "--------------------------------------------------------------------------------"

%.x :
	@if test -d src; then \
	(cd src ; $(MAKE) TLDEPS= $*.abc || exit 1) ; fi
	@echo "--------------------------------------------------------------------------------"


all: libs src

clean : 
	rm -f *~ core;
	@for dir in \
		$(DIRS) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done

veryclean : clean uninstall
	rm -f *~ core;
	@for dir in \
		$(DIRS) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= veryclean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done
	

install : 
	@for dir in \
		$(DIRS) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= install ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done

uninstall :
	@for dir in \
		$(DIRS) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= uninstall ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done








