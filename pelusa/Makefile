#
#
#
SHELL = /bin/bash
GCC = g++
SRC_DIR = src
BUILD_DIR = build
TARGET = pelusa

CCOPT = -g
LIBS = -lboost_program_options -lpthread -L$(BUILD_DIR)
# DEFINES = -DDB_SHORT -DREAD_36 -DMAXGAP=3 -DMAXHITS=10000
# INCLUDES = -Isrc/libsoap
ifneq ($(shell uname -s),Darwin)
    STATIC   = -static
endif

OBJS = $(patsubst %.cpp,%.o, $(wildcard $(SRC_DIR)/*.cpp))

all: init $(TARGET)

init:
	@( [[ ! -e $(BUILD_DIR) ]] && mkdir $(BUILD_DIR) ) || true

$(TARGET): $(OBJS)
	$(GCC) $(DEFINES) $(INCLUDES) $(CCOPT) $(OBJS) $(STATIC) -o $(BUILD_DIR)/$@ $(LIBS)

%.o: %.cpp
	$(GCC) $(DEFINES) $(INCLUDES) $(CCOPT) -c $< -o $@

INSTALL_DIR ?= .
install:
	/usr/bin/install -m 555 $(BUILD_DIR)/$(TARGET) $(INSTALL_DIR)/

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(OBJS)
