LDLIBS=-lm
CFLAGS=-g

data_res1.dat: freq_tracer
	./freq_tracer < 20130517_47074.osc > data_res1.dat



#all: data_res.dat
#	./data.gp

data_res.dat: freq_tracer data.dat
	./freq_tracer < data.dat > data_res.dat

data.dat data_test.dat: freq_mktest
	./freq_mktest

freq_tracer:freq_tracer.c
freq_mktest:freq_mktest.c