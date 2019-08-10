#
#  This file is part of MUMPS 5.0.1, released
#  on Thu Jul 23 17:08:29 UTC 2015
#
topdir = .
libdir = $(topdir)/lib

default:	dexamples

.PHONY: default alllib all c z s d \
	sexamples dexamples cexamples zexamples \
	mumps_lib requiredobj libseqneeded clean

alllib:		c z s d
all:		cexamples zexamples sexamples dexamples

c:
	$(MAKE) ARITH=c mumps_lib
z:
	$(MAKE) ARITH=z mumps_lib
s:
	$(MAKE) ARITH=s mumps_lib
d:
	$(MAKE) ARITH=d mumps_lib


# Is Makefile.inc available ?
Makefile.inc:
	@echo "######################################################################"
	@echo "# BEFORE COMPILING MUMPS, YOU SHOULD HAVE AN APPROPRIATE FILE"
	@echo "# Makefile.inc AVALAIBLE. PLEASE LOOK IN THE DIRECTORY ./Make.inc FOR"
	@echo "# EXAMPLES OF Makefile.inc FILES, AT Make.inc/Makefile.inc.generic"
	@echo "# IN CASE YOU NEED TO BUILD A NEW ONE AND READ THE MAIN README FILE"
	@echo "######################################################################"
	@exit 1

include Makefile.inc

mumps_lib: requiredobj
	(cd src ; $(MAKE) $(ARITH))

cexamples:	c
	(cd examples ; $(MAKE) c)

zexamples:	z
	(cd examples ; $(MAKE) z)

sexamples:	s
	(cd examples ; $(MAKE) s)

dexamples:	d
	(cd examples ; $(MAKE) d)

requiredobj: Makefile.inc $(LIBSEQNEEDED) $(libdir)/libpord$(PLAT)$(LIBEXT)

# dummy MPI library (sequential version)

libseqneeded:
	(cd libseq; $(MAKE))

# Build the libpord.a library and copy it into $(topdir)/lib
$(libdir)/libpord$(PLAT)$(LIBEXT):
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); \
	  $(MAKE) CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC="$(OUTC)" LIBEXT=$(LIBEXT); \
	fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cp $(LPORDDIR)/libpord$(LIBEXT) $@; \
	fi;

clean:
	(cd src; $(MAKE) clean)
	(cd examples; $(MAKE) clean)
	(cd $(libdir); $(RM) *$(PLAT)$(LIBEXT))
	(cd libseq; $(MAKE) clean)
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); $(MAKE) realclean; \
        fi;

