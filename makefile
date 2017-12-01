Iflags = -I/opt/local/include/
Lflags = -L/opt/local/lib/
cflags = -Wno-deprecated -O3
lflags = -lgsl -lm 

hysim: hysim.o
	g++ $(cflags) $(Iflags) $(Lflags) $(lflags) hysim.o -o hysim

hysim.o: hysim.c
	g++ $(cflags) $(Iflags) -c hysim.c
