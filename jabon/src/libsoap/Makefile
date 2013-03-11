#
#
#
SHELL = /bin/bash
GCC = g++
SRC_DIR = src
BUILD_DIR = ../../build
TARGET = libsoap.a

SOAP_FLAGS= -static -DMAXGAP=3 -DMAXHITS=10000 -DTHREAD -O3
OLIGOLEN= -DREAD_36
# REF_CHR=  -DDB_CHR
REF_CHR=  -DDB_SHORT
CCOPT = -g

LIBS = -lpthread

SOURCE = align dbseq pairs param reads utilities
OBJS = $(patsubst %,%.chr.o,$(SOURCE))

all: init $(BUILD_DIR)/$(TARGET)

init:
	@( [[ ! -e $(BUILD_DIR) ]] && mkdir $(BUILD_DIR) ) || true

$(BUILD_DIR)/$(TARGET): $(OBJS)
	ar cru $(BUILD_DIR)/$(TARGET) $(OBJS) 

%.chr.o: %.cpp
	$(GCC) $(CCOPT) $(SOAP_FLAGS) $(REF_CHR) $(OLIGOLEN) -c $< -o $@
	
clean:
	rm -f $(BUILD_DIR)/$(TARGET)
	rm -f $(OBJS)
