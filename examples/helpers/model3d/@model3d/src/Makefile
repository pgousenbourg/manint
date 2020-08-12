CXX = g++4
CXXFLAGS = -O2 -Wall -DMATLAB_MEX_FILE -pipe

ARCH = $(shell uname -p)
ifeq ($(ARCH),x86_64)
		CXXFLAGS += -fPIC
		SUFFIX = mexa64
else
	SUFFIX = mexglx
endif


TARGET = ../model3d.$(SUFFIX)
TARGET2 = ../serialize.$(SUFFIX)
SOURCES = dxf.cpp mexfunction.cpp model3d.cpp c3ds.cpp
MATDIR = /usr/local/matlab
INCDIR = -I. -I$(MATDIR)/extern/include

all : $(TARGET) $(TARGET2)

$(TARGET2) : serialize.o model3d.o 
	$(CXX) $(CXXFLAGS)  -shared $^ -o $@

$(TARGET) : $(patsubst %.cpp,%.o,$(SOURCES))
	$(CXX) $(CXXFLAGS)  -shared $^ -o $@

%.o : %.cpp %.h
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $<

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCDIR) -c $<

clean : 
	rm -rf $(TARGET) $(patsubst %.cpp,%.o,$(SOURCES)) serialize.o
