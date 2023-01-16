CC = g++
CFLAGS = -Wall 

main: main.cpp
	$(CC) $(CFLAGS) -o main main.cpp 

run:
	./main

clean:
	$(RM) main



