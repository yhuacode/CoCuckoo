all: ycsb example

ycsb: ../CoCuckoo.h ycsb.cc ../MurmurHash3.cpp
	g++ -std=c++11 -Wunused ycsb.cc ../MurmurHash3.cpp -O3 -o ycsb -lpthread

example: ../CoCuckoo.h example.cc ../MurmurHash3.cpp
	g++ -std=c++11 -Wunused example.cc ../MurmurHash3.cpp -O3 -o example -lpthread

clean:
	rm -f ycsb example thread-*.log

.PHONY: all clean