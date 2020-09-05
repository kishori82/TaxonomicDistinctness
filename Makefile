IDIR =../libStatGen/include

CC = g++ -std=c++11
CFLAGS =  

ODIR=obj

TARGET = taxdistinct
LIBS = -lz 

_DEPS =  utilities.h taxdistinct.h tax_distinct_types.h
#DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ =  utilities.o main.o taxdistinct.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)  $(_DEPS)
	$(CC) -o $@ $^ $(CFLAGS)  $(LIBS) 

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
