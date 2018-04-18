#Library makefile

include make.inc

LIBDIR=$(PREFIX)/lib
DIRS=cjson cmpi cterm filemanager qepp

default		: all

.PHONY		: qepp clean veryclean install uninstall

all		: qepp

qepp		: 
	@if test -d qepp; then \
	(cd qepp ; $(MAKE) TLDEPS= all || exit 1) ; fi
	@echo "--------------------------------------------------------------------------------"

clean : 
	@for dir in \
		$(DIRS) \
	; do \
		if test -d $$dir ; then \
			( cd $$dir ; \
			$(MAKE) TLDEPS= clean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
		 fi \
	done

veryclean : 
	@for dir in \
		$(DIRS) \
	; do \
		if test -d $$dir ; then \
			( cd $$dir ; \
			$(MAKE) TLDEPS= veryclean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
		 fi \
	done

install :
	@for dir in \
		$(DIRS) \
	; do \
		if test -d $$dir ; then \
			( cd $$dir ; \
			$(MAKE) TLDEPS= install ) ; \
		echo "--------------------------------------------------------------------------------" ; \
		 fi \
	done

uninstall :
	@for dir in \
		$(DIRS) \
	; do \
		if test -d $$dir ; then \
			( cd $$dir ; \
			$(MAKE) TLDEPS= uninstall ) ; \
		echo "--------------------------------------------------------------------------------" ; \
		fi \
	done
	#find $(LIBDIR) -type l -delete





