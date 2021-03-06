TARGET = rmc
 
SRCS = $(TARGET).cpp
OBJS = $(TARGET).o
 
 
CXXFLAGS   = `root-config --cflags` -g  -std=c++11  -Wl,-rpath,${ROOTSYS}/lib -O1 
CXXLIBS    = `root-config --libs` -g
CC = g++ 
 
$(TARGET): $(OBJS)
	$(CC) $(CXXLIBS) $(OBJS) -o $@
 
# suffix rule
.cc.o:
	$(CC) $(CXXFLAGS) -c $<
 
# clean
clean:
	rm -f $(TARGET) $(OBJS)

