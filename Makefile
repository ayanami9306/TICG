CC = g++
TARGET = main
CXXFLAGS = -Wall -O2 -std=c++17
OBJS = main.o TICG_adjust.o TICG_engine.o TICG_hopping.o TICG_initialize.o TICG_io.o TICG_objects.o TICG_op.o

%.o: %.cpp
	$(CC) $(CXXFLAGS) -c $<

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $(OBJS) -o $(TARGET) 

clean:
	rm -f $(OBJS)
