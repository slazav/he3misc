all: 	host   node

HOBJS =	host.o

NOBJS =	node.o btnlib.o chklib.o d1mach.o

.SUFFIXES: .f .o

host.o:	host.f
	f77 -c -vx host.f

node.o: node.f
	f77 -c -vx node.f

.f.o:
	f77 -c -vx $<

host:	$(HOBJS)
	f77 $(HOBJS) -o host -host -vec

node:	$(NOBJS)
	f77 $(NOBJS) -o node -node -vec
