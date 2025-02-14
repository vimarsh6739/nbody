CILK_C= opencilk-build/bin/clang
CILK_CXX= opencilk-build/bin/clang++
CILK_FLAGS= -fopencilk
LDFLAGS= 
CLFAGS= -O2 -Wall

.PHONY: all
all: fib

fib: fib.o
	$(CILK_C) $(CILK_FLAGS) $^ -o $@ $(LDFLAGS)

fib.o: fib.c
	$(CILK_C) $(CFLAGS) $(CILK_FLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o
	rm -f fib
