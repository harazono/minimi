DEPENDS := .depend.mk
CXXFLAGS := -g -Wall -O3 -DNDEBUG -std=c++11 -fopenmp -lgomp
# CXXFLAGS := -g -Wall -O0 -std=c++11 -fopenmp -lgomp
GTESTFLAGS := -lgtest -lpthread
SOURCES := count_kmer.cpp
TEST_SOURCES := kmer_library_test.cpp
CXX := $(HOME)/local/bin/g++

.PHONY: depend clean test build_test

all: count_kmer minimi


count_kmer: count_kmer.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

test: build_test
	./kmer_library_test

build_test: kmer_library_test

kmer_library_test: kmer_library_test.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ $(GTESTFLAGS)

depend:
	$(CXX) -MM $(SOURCES) $(TEST_SOURCES) > $(DEPENDS)

clean:
	-rm *.o

include $(DEPENDS)

minimi: minimi.o
	$(CXX) align_NW.cpp -o minimi -DMAIN

