
CC=g++
CFLAGS=-I. -D__EXPORT="" -Dnullptr="0"

ODIR=obj
LDIR =../lib

LIBS=-lm

#_DEPS = test.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main_closed_loop_float.o estimator_23states.o estimator_utilities.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

#$(DEPS)
$(ODIR)/%.o: %.cpp
	mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS)

#
estimator_closed_loop_test: $(OBJ)
	g++ -O3 -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean plots

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

plots:
	python plot_attitude.py &
	python plot_position.py &
	python plot_states.py &
