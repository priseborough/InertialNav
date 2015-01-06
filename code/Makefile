
CXX?=g++
EXTRAARGS?=
CFLAGS=-I. -D__EXPORT="" -Dnullptr="0" -fno-strict-aliasing \
			   -fno-strength-reduce \
			   -fomit-frame-pointer \
   			   -funsafe-math-optimizations \
   			   -Wall \
			   -Wextra \
			   -Wshadow \
			   -Wfloat-equal \
			   -Wframe-larger-than=1024 \
			   -Wpointer-arith \
			   -Wmissing-declarations \
			   -Wpacked \
			   -Wno-unused-parameter \
			   -Werror=format-security \
			   -Werror=array-bounds \
			   -Wfatal-errors \
			   -Wno-missing-field-initializers \
			   -std=gnu++0x

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
	$(CXX) -c -o $@ $< $(CFLAGS)

#
estimator_closed_loop_test: $(OBJ)
	$(CXX) -O3 -o $@ $^ $(CFLAGS) $(LIBS) $(EXTRAARGS)

.PHONY: clean plots

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~

plots:
	python plot_attitude.py &
	python plot_position.py &
	python plot_states.py &
	python plot_velocity.py &

tests:
	make clean
	export CXX="g++"
	export EXTRAARGS=""
	make
	unzip -o InputFilesPX4.zip
	./estimator_closed_loop_test
#make clean
#export CXX=arm-none-eabi-gcc
#export EXTRAARGS=--specs=nosys.specs
#make
#export CXX="g++"
#export EXTRAARGS=""
	make clean
	./validate_output.py
