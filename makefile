# Compiler and Flags
CXX      := g++
CXXFLAGS := $(shell pkg-config --cflags ibex) -frounding-math -ffloat-store
LIBS     := $(shell pkg-config --libs ibex) -libex -lprim -lsoplex -lz -lClp -lCoinUtils -lm

# Debug Mode
ifeq ($(DEBUG), yes)
    CXXFLAGS += -O0 -g -pg -Wall
else
    CXXFLAGS += -O3 -DNDEBUG -Wno-deprecated -Wno-logical-op-parentheses -g
endif

# Source and Header Files
SRCS     := simulation.cpp C_STL_adaptative.cpp
BINS     := simulation.out

# Default Target
all: $(BINS)

# Rule for building the executable
$(BINS): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS) $(LIBS)

# Clean Up
clean:
	rm -f $(BINS)

