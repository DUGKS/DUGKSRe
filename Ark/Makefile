vpath %.h Header:Lib:Constant
vpath %.cpp Cpp:Cpp/Mesh:Lib

Target = ACDUGKS
QModel = D2Q9AC

CC = g++ 

CFLAGS = -std=c++11 -Wall -Wextra -I./Header -I./Lib -I./Constant -O3

MPFLAG = -fopenmp

ROOT_SRC = ./Cpp

DIR_SRC = $(shell find ${ROOT_SRC} -maxdepth 1 -type d)

SRC = $(foreach dir,$(DIR_SRC),$(wildcard $(dir)/*.cpp))

DIR_OBJ = ./Obj

DIR_OUT = ../FlowField

#../FlowField/global ../FlowField/Convergence

subObjs = $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))

objects = $(subObjs) ${DIR_OBJ}/${QModel}.o

# ALL:

# 	@echo "SRC:" $(notdir ${SRC})

$(Target):$(objects)

	$(CC) $(objects) -o $@ ${MPFLAG}

${DIR_OBJ}/%.o:%.cpp

	$(CC) $(CFLAGS) -c $< -o $@ ${MPFLAG}

${DIR_OBJ}/${QModel}.o:./Lib/${QModel}.cpp

	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	find -name "${Target}" -exec mv {} ~/.trash/ \;
	find ${DIR_OBJ} -name "*.o" -type f -exec mv {} ~/.trash/ \;
	find ${DIR_OUT} -name "*.dat" -type f -exec mv {} ~/.trash/ \; 
	find ${DIR_OUT} -name "*.ark" -type f -exec mv {} ~/.trash/ \;
	find ${DIR_OUT} -name "tecplot-derived" -type d -exec mv -f {} ~/.trash/ \; 

# headers = ZeroDeclaration.h Mesh_2D.h ZeroConstant.h $(QModel).h ZeroReference.h ZeroFlip.h

# main : $(objects)
# 	$(CC) -o main $(objects) -fopenmp

# main.o : main.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 main.cpp

# MeshConstruct.o : MeshConstruct.cpp MeshConstructFunc.h $(headers)
# 	$(CC) -c -O3 -std=c++11 MeshConstruct.cpp

# Preprocess.o : Preprocess.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Preprocess.cpp

# DUGKS2DSolver.o : DUGKS2DSolver.cpp NewtonCotes.h $(headers)
# 	$(CC) -c -O3 -std=c++11 DUGKS2DSolver.cpp -fopenmp

# $(QModel).o : $(QModel).cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 $(QModel).cpp

# Mesh_2D.o : Mesh_2D.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Mesh_2D.cpp

# Output.o : Output.cpp $(headers)
# 	$(CC) -c -O3 -std=c++11 Output.cpp

# .PHONY : clean
# clean :
# 	rm pexe $(objects)
