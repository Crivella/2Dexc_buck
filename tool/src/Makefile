#exec makefile

include ../make.inc

DIR = $(shell find . -maxdepth 1 -type d | grep -v obj | grep "/")
C_FILES = $(shell find . -maxdepth 1 -name "*.c")
OBJ_FILES := $(addprefix $(ODIR)/,$(notdir $(C_FILES:.c=.o)))
X_FILES := $(notdir $(C_FILES:.c=.x))

DIRECTORY = $(shell echo "$(DIR)" | grep -v obj)

.PHONY : fold all_subdir all_exe all clean veryclean install uninstall

fold : 
	mkdir -p $(ODIR);

all_subdir:
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		cp ../Makefile .; \
		cp ../make.inc .; \
		$(MAKE) TLDEPS= all ) \
	done

all_exe : fold all_subdir $(X_FILES)

%.x : $(ODIR)/%.o
	@if test -f $*.c; then \
		echo "--------------------------------------------------------------------------------"; \
		echo "$(MPICC) -o $@ $< $(CFLAGS) $(LD)"; \
		$(MPICC) -o $@ $< $(CFLAGS) $(LD); \
		echo "--------------------------------------------------------------------------------\n\n"; \
	fi

$(ODIR)/%.o : %.c
	@if test -f $^; then \
		echo "$(MPICC) -c -o $@ $< $(CFLAGS)"; \
		$(MPICC) -c -o $@ $< $(CFLAGS); \
	fi

%.abc :
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= $@ ) \
	done
	@if test -f $*.c; then \
		echo "$(MPICC) -o $*.x $*.c $(CFLAGS) $(LD)"; \
		$(MPICC) -o $*.x $*.c $(CFLAGS) $(LD); \
		echo "\n\n"; \
	fi

all : all_exe

clean : 
	rm -f $(ODIR)/*.o *~ core; \
	rm -fr obj;
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done

veryclean: clean uninstall
	rm -f *.x
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= veryclean ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done

install : all
	@mkdir -p $(BDIR)
	ln -sf $(CURDIR)/*.x $(BDIR)/
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= install ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done

uninstall :
	@for file in \
		$(X_FILES) \
	; do \
		if test -h $(BDIR)/$$file; then \
			rm $(BDIR)/$$file; \
		fi ; \
		echo "--------------------------------------------------------------------------------" ; \
	done
	@for dir in \
		$(DIRECTORY) \
	; do \
		( cd $$dir ; \
		$(MAKE) TLDEPS= uninstall ) ; \
		echo "--------------------------------------------------------------------------------" ; \
	done



