CC = g++ -std=c++11 -Wall
CFLAGS =  

ODIR=obj

TARGET = taxdistinct
LIBS = -lz 

_DEPS =  utilities.h options.h taxdistinct.h metataxa_types.h

OBJ =  utilities.o options.o main.o taxdistinct.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.c++ $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)  $(_DEPS)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
