#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#/*                                                                           */
#/*    This file is part of the program BACS                                  */
#/*                                                                           */
#/*    an implementation of a Branch-And-Cut algorithm to solve the           */
#/*    stable set problem.                                                    */
#/*                                                                           */
#/*    Copyright (C) 2024-  Discrete Optimization Group, TU Darmstadt         */
#/*                                                                           */
#/*                                                                           */
#/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
#/*    Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                    */
#/*                                                                           */
#/*    Both are licensed under the Apache License, Version 2.0.               */
#/*                                                                           */
#/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#----------------------------------------------------------------------------
# paths
#----------------------------------------------------------------------------

# SCIP_PATH is an enivironment variable that should be set outside
# this makefile and point to the main directory of SCIP.

SCIPDIR         =       $(SCIP_PATH)
SCIPREALPATH	=	$(realpath $(SCIPDIR))


#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------

include $(SCIPDIR)/make/make.project
-include $(SCIPDIR)/make/local/make.$(HOSTNAME)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP)
-include $(SCIPDIR)/make/local/make.$(HOSTNAME).$(COMP).$(OPT)

#-----------------------------------------------------------------------------
# basic definitions
#-----------------------------------------------------------------------------

VERSION		=	1.0
TEST		=	short

# overwrite flags for dependencies
DFLAGS          =       -MMD

USEBESTSOL	=	false

GXXWARN		+= 	-Wzero-as-null-pointer-constant

# detect with symmetry code has been used by SCIP
SCIPVERSION	:=	$(shell $(SCIPDIR)/bin/scip.$(BASE).$(LPS).$(TPI)$(EXEEXTENSION) -v | sed -e 's/$$/@/')
override SYM:=	$(shell echo "$(SCIPVERSION)" | sed -e 's/.* SYM=\([^@]*\).*/\1/')

# detect symmetry handling code
ifeq ($(SYM),nauty)
USRFLAGS	+=	-DBACS_WITH_NAUTY
ifeq ($(NAUTYEXTERNAL),true)
USRLDFLAGS	+=	$(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)nauty.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif
endif
ifeq ($(SYM),snauty)
USRFLAGS	+=	-DBACS_WITH_NAUTY
ifeq ($(NAUTYEXTERNAL),true)
USRLDFLAGS	+=	$(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)nauty.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif
endif
ifeq ($(SYM),bliss)
USRFLAGS	+=	-DBACS_WITH_BLISS -I$(SCIPDIR)/$(LIBDIR)/include/
USRLDFLAGS	+=	$(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)bliss.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif
ifeq ($(SYM),sbliss)
USRFLAGS	+=	-DBACS_WITH_SBLISS -I$(SCIPDIR)/$(LIBDIR)/include/
USRLDFLAGS	+=	$(LINKCC_L)$(SCIPREALPATH)/$(LIBDIR)/$(LIBEXTTYPE) $(LINKCC_l)bliss.$(OSTYPE).$(ARCH).$(COMP)$(LINKLIBSUFFIX)
endif


#-----------------------------------------------------------------------------
# main program
#-----------------------------------------------------------------------------

MAINNAME	=	bacs

MAINCXXOBJ	=	probdata_bacs.o cliquepartition.o graphpresolve.o symmetry.o  \
			branch_degree.o branch_maxlpneigh.o branch_cliquepartition.o cons_clique.o \
			heur_dynamicdeg.o heur_greedydeg.o heur_greedylp.o heur_greedyrounding.o heur_tabu.o heur_tabuseq.o \
			presol_persistence.o presol_stableprobing.o \
			prop_cliquefixing.o prop_comp.o prop_dominance.o prop_lowdegreenodes.o prop_neighborhoods.o prop_simplicial.o \
			sepa_neigh.o main.o
MAINCOBJ	=	transfer_stats.o changesoltime.o

MAINOBJ		=	$(MAINCXXOBJ) $(MAINCOBJ)

MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINCXXOBJ:.o=.cpp))
MAINSRC		+=	$(addprefix $(SRCDIR)/,$(MAINCOBJ:.o=.c))

MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

# other targets
COMPLEMENTGRAPHOBJ   =       complementgraph.o
COMPLEMENTGRAPHSRC   =       $(addprefix $(SRCDIR)/,$(COMPLEMENTGRAPHOBJ:.o=.cpp))
COMPLEMENTGRAPHOBJFULES =    $(addprefix $(OBJDIR)/,$(COMPLEMENTGRAPHOBJ))
COMPLEMENTGRAPHBIN   =       $(BINDIR)/complementgraph.$(BASE)

#-----------------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------------

# running the makefile without arguments will build everything
.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

complementgraph: $(COMPLEMENTGRAPHBIN)

# create directory for objective files
$(OBJDIR):
		@-mkdir -p $(OBJDIR)

# create directory for binaries
$(BINDIR):
		@-mkdir -p $(BINDIR)

# possibly compile SCIP
.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

# create documention
.PHONY: doc
doc:
		cd doc; $(DOXY) $(MAINNAME).dxy

# call lint code checker on source files
.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/lint main-gcc.lnt +os\(lint.out\) -u -zero \
			$(USRFLAGS) $(FLAGS) -I/usr/include -UNDEBUG -USCIP_WITH_READLINE -USCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'

# call pclint code checker on source files
.PHONY: pclint
pclint:		$(MAINSRC)
		-rm -f pclint.out
ifeq ($(FILES),)
		@$(SHELL) -ec 'echo "-> running pclint ..."; \
			for i in $^; \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
else
		@$(SHELL) -ec  'echo "-> running pclint on specified files ..."; \
			for i in $(FILES); \
			do \
				echo $$i; \
				$(PCLINT) -I$(SCIPREALPATH) -I$(SCIPREALPATH)/pclint main-gcc.lnt +os\(pclint.out\) -b -u -zero \
				$(USRFLAGS) $(FLAGS) $(SDPIINC) -uNDEBUG -uSCIP_WITH_READLINE -uSCIP_ROUNDING_FE -D_BSD_SOURCE $$i; \
			done'
endif

# clean all objective files
.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		@-(rm -f $(OBJDIR)/*.o && rm -f $(OBJDIR)/*.d && rmdir $(OBJDIR));
		@echo "-> remove objective files"
endif

# create tags
.PHONY: tags
tags:
		rm -f TAGS; etags src/*.cpp src/*.h $(SCIPDIR)/src/*/*.c $(SCIPDIR)/src/*/*.h --output=TAGS; sed 's!\#undef .*!!g' TAGS > tags; mv tags TAGS

# run tests
.PHONY: test
test:		$(MAINFILE)
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) $(DISPFREQ) $(CONTINUE) $(LOCK) $(VERSION) $(LPS) $(USEBESTSOL) $(SETCUTOFF);

# run tests on cluster
.PHONY: testcluster
testcluster:	$(MAINFILE)
		cd check; \
		$(SHELL) ./check_cluster.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)).$(HOSTNAME) $(TIME) $(NODES) $(MEM) \
		$(DISPFREQ) $(VERSION) $(LPS) $(SETCUTOFF) $(QUEUE) $(PPN) $(CLIENTTMPDIR) $(EXCLUSIVE);

#-----------------------------------------------------------------------------
# project rules
#-----------------------------------------------------------------------------

# do not attempt to include .d files if there will definitely be any (empty DFLAGS), because it slows down the build on Windows considerably
ifneq ($(DFLAGS),)
-include $(MAINOBJFILES:.o=.d)
else
ifeq ($(VERBOSE),true)
$(info No compilation dependencies. If changing header files, do a make clean before building.)
endif
endif

# turn off output unless VERBOSE is true
ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK) $(COMPLEMENTGRAPHBIN) $(COMPLEMENTGRAPHOBJFULES)
endif

# main target
$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> using SYM=$(SYM)"
		@echo "-> linking $@"
		$(LINKCXX) $(MAINOBJFILES) $(LINKCXXSCIPALL) $(LINKCXX_o)$@

# create short link
$(MAINSHORTLINK): $(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

# target for running tclique
$(COMPLEMENTGRAPHBIN): $(OBJDIR) $(OBJDIR)/complementgraph.o
		@echo "-> linking $@"
		$(LINKCXX) -o $(COMPLEMENTGRAPHBIN) $(OBJDIR)/complementgraph.o $(LINKCXXSCIPALL) $(LDFLAGS)


# rules to build objective files from source code
$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) $(DFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) $(DFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
