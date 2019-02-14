TARGET = rmc
 
SRCS = $(TARGET).cpp
OBJS = $(TARGET).o
 
 
CXXFLAGS   = `root-config --cflags` -Wl,-rpath,${ROOTSYS}/lib -O1
CXXLIBS    = `root-config --libs`
CC = g++ 
 
$(TARGET): $(OBJS)
	$(CC) $(CXXLIBS) $(OBJS) -o $@
 
# suffix rule
.cc.o:
	$(CC) $(CXXFLAGS) -c $<
 
# clean
clean:
	rm -f $(TARGET) $(OBJS)
