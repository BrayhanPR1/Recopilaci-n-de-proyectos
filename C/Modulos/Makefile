CC = gcc
GFLAGS = -Wall
LFLAGS = -lm
OBJS = main.o allvars.o funcion1.o
EXEC = simulador_n_cuerpos.out

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LFLAGS) -o $(EXEC)
	rm $(OBJS)

clean:
	rm *~ *.out
