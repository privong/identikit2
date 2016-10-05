# ~/coll/idkit2/export/Makefile: build Identikit2 program.

PRECISION = MIXEDPREC

VAR =

# NOTE: On MACOSX the -lXmu and -lX11 flags are not necessary if the
# native version of libglut is used.

idkit: idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) \
	  -DVARIANT=\"\" -L/usr/X11R6/lib -o idkit \
	  idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 -lgsl -lm -lcblas

idkit_bulgeless: idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) \
      -DBULGELESS \
	  -DVARIANT=\"\" -L/usr/X11R6/lib -o idkit_bulgeless \
	  idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 -lgsl -lm -lcblas

idkit_$(VAR): idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) $(OPTIONS) \
	  -DVARIANT='" $(VAR)"' -L/usr/X11R6/lib -o idkit_$(VAR) \
	  idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 -lgsl -lm

idkit_nw: idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) $(OPTIONS) \
	  -DNO_WEIGHTING -DVARIANT='" nw"' -L/usr/X11R6/lib -o idkit_nw \
	  idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 -lgsl -lm

idkit_nn: idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) $(OPTIONS) \
	  -DNO_NORMALIZE -DVARIANT='" nn"' -L/usr/X11R6/lib -o idkit_nn \
	  idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 -lgsl -lm

idkit_nwn: idkit.c
	$(ZCC) $(ZCCFLAGS) $(ZLDFLAGS) -D$(PRECISION) $(OPTIONS) \
	  -DNO_WEIGHTING -DNO_NORMALIZE -DVARIANT='" nwn"' -L/usr/X11R6/lib \
	  -o idkit_nwn idkit.c -lNBody -lClib -lglut -lGLU -lGL -lXmu -lX11 \
	  -lgsl -lm
