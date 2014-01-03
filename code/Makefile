
CC=g++
CFLAGS=-I. -D__EXPORT="" -Dnullptr="0"

ODIR=obj
LDIR =../lib

LIBS=-lm

#_DEPS = test.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main_closed_loop_float.o estimator.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#$(DEPS)
$(ODIR)/%.o: %.cpp
	mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

#
estimator_closed_loop_test: $(OBJ)
	g++ -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~