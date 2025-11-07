# Compiler setup
CC      ?= cc
CFLAGS  ?= -O2 -Wall -std=c11
LDFLAGS ?=
LDLIBS  ?= -lm

# If you installed OpenBLAS (recommended on macOS M1):
# brew install openblas
BREW := $(shell command -v brew 2>/dev/null)
ifneq ($(BREW),)
  OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null)
  ifneq ($(OPENBLAS_PREFIX),)
    CFLAGS  += -I$(OPENBLAS_PREFIX)/include
    LDFLAGS += -L$(OPENBLAS_PREFIX)/lib -Wl,-rpath,$(OPENBLAS_PREFIX)/lib
    LDLIBS  += -lopenblas
  endif
endif

# Targets
all: Lattice_Diffusion

Lattice_Diffusion: Lattice_Diffusion.o find_neighbor.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o Lattice_Diffusion