
# note the TARFILE_NAME embeds the release version number
TARFILE_NAME	= nifti2clib-0.0.1

USEZLIB         = -DHAVE_ZLIB

## Compiler  defines
CC		= gcc
#IFLAGS          = -I/home/sapjw12/dev/include/gsl
IFLAGS          = -I/home/sapjw12/dev/include
LFLAGS		= -L/home/sapjw12/dev/lib
#CFLAGS          = -g -Q -v -Wall -std=gnu99 -pedantic $(USEZLIB) $(IFLAGS)
CFLAGS		= -g -Wall -std=gnu99 -pedantic $(USEZLIB) $(IFLAGS)

LLIBS 		= -lz -lm -lgsl -lgslcblas

MISC_OBJS	= nifticdf.o znzlib.o
OBJS	   	= nifti2_io.o $(MISC_OBJS)

# List my own programs
#MYXS		= MEfit MEfitNL MEfitLogLin MEfitLM
MYXS		= MEfitLM

# --------------------------------------------------
# default compile for C files
#%.o : %.c %.h
#	$(CC) -c $(CFLAGS) $< -o $@

# --------------------------------------------------

all: $(MYXS)

nifti2objs: $(OBJS)

#MEfit: MEfit.o nifti2objs
#	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $< $(OBJS) $(LLIBS)

#MEfitNL: MEfitNL.o nifti2objs
#	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $< $(OBJS) $(LLIBS)

#MEfitLogLin: MEfitLogLin.o nifti2objs
#	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $< $(OBJS) $(LLIBS)

MEfitLM: MEfitLM.o nifti2objs
	$(CC) -o $@ $(CFLAGS) $(LFLAGS) $< $(OBJS) $(LLIBS)

clean:
	$(RM) *.o $(MYXS) 

