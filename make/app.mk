CC=g++

SRC_DIR := $(TOPDIR)/src/solver
OBJ_DIR := $(TOPDIR)/obj
BIN_DIR := $(TOPDIR)/app

EXE := $(BIN_DIR)/euler
SRC := $(wildcard $(SRC_DIR)/*.cpp $(SRC_DIR)/*/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o ) 


CXXFLAGS := -std=c++17 -O3  -I$(SRC_DIR) -MMD -MP 
CPFLAGS   := -Wall
LDFLAGS  := -Llib 
LDLIBS   := -lm

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ) | $(BIN_DIR)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) $(CPFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(EXE) $(OBJ_DIR)

-include $(OBJ:.o=.d)