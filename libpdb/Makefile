#
#	For machines that don't have lint, remove $(LINTLIB) from the
#	"all" target.
#
#	On BSD machines, LINTOPT should be '-C', RANLIB should be 'ranlib'
#
#	On System V machines, LINTOPT should be '-o $(NOTHING)',
#	RANLIB should be ':'
#
SHELL		= /bin/sh
RANLIB		= :
NOTHING		=
LINTOPT		= -o $(NOTHING)

OPT		= -O -g
DEBUG		=
CFLAGS		= $(OPT) $(DEBUG)
LIBRARY		= pdb

LIBARCH		= lib$(LIBRARY).a
LINTLIB		= llib-l$(LIBRARY).ln
OBJS		= pdb_read.o pdb_sprntf.o pdb_sscanf.o pdb_write.o pdbrun.o ms.o
SRCS		= pdb_read.c pdb_sprntf.c pdb_sscanf.c pdb_write.c pdbrun.c ms.c

all:		$(LIBARCH) $(LINTLIB)

install:	$(LIBARCH)
		install -F /usr/local/lib $(LIBARCH)

$(LIBARCH):     $(OBJS)
		@echo "Loading $(LIBARCH) ... "
		@ar cru $(LIBARCH) $(OBJS)
		@$(RANLIB) $(LIBARCH)
		@echo "done"

$(LINTLIB):	$(SRCS) $(HDRS)
		@echo "Creating $(LINTLIB) ... "
		@lint $(LINTOPT)$(LIBRARY) $(SRCS)
		@echo "done"

clean:;		@rm -f $(OBJS)

spotless:;	@rm -f $(OBJS) $(LIBARCH) $(LINTLIB)

lint:		$(SRCS)
		lint -u $(SRCS)

tags:           $(HDRS) $(SRCS); @ctags -w $(HDRS) $(SRCS)
