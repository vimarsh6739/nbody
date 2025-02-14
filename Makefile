CILK_C= $(realpath opencilk-build/bin/clang)
CILK_CXX= $(realpath opencilk-build/bin/clang++)
CILK_FLAGS= -fopencilk
LDFLAGS= 
CLFAGS= -O2 -Wall

.PHONY: all
all: fib qsort

qsort: qsort.o
	$(CILK_C) $(CILK_FLAGS) $^ -o $@ $(LDFLAGS)

qsort.o: qsort.c
	$(CILK_C) $(CFLAGS) $(CILK_FLAGS) -c $< -o $@

fib: fib.o
	$(CILK_C) $(CILK_FLAGS) $^ -o $@ $(LDFLAGS)

fib.o: fib.c
	$(CILK_C) $(CFLAGS) $(CILK_FLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o
	rm -f fib qsort
