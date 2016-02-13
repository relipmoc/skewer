CXX=g++
DEBUG?=0
CXXFLAGS?=-c
LDFLAGS=-pthread
LDLIBS=-lrt
ifneq ($(OS), Windows_NT)
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S), Darwin)
		LDFLAGS=-lpthread
		LDLIBS=
	endif
endif
SRC=src
SOURCES=$(SRC)/main.cpp $(SRC)/parameter.cpp $(SRC)/matrix.cpp $(SRC)/fastq.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=skewer

ifeq ($(DEBUG), 1)
	CXXFLAGS += -O0 -g -Wall -D_DEBUG
else
	CXXFLAGS += -O2
endif

.PHONY: all debug clean

all:$(EXECUTABLE)

debug:
	$(MAKE) $(MAKEFILE) DEBUG=1

$(EXECUTABLE):$(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@ $(LDLIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

$(SRC)/main.o: $(SRC)/parameter.h $(SRC)/matrix.h $(SRC)/fastq.h $(SRC)/common.h
$(SRC)/parameter.o: $(SRC)/parameter.h $(SRC)/fastq.h $(SRC)/common.h
$(SRC)/matrix.o: $(SRC)/matrix.h $(SRC)/common.h
$(SRC)/fastq.o: $(SRC)/fastq.h $(SRC)/common.h

# Clean
clean:
	rm -rf $(OBJECTS) $(EXECUTABLE)

# Install
install:
	mv -f $(EXECUTABLE) /usr/local/bin
