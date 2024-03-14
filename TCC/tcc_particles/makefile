OUTPUT=dowloading_wall_restitution_OMP

CC=gcc
CC_OPT=-std=c11

CC_PTH=-lgsl -lgslcblas -lm

.PHONY: all
all: $(OUTPUT)

$(OUTPUT): $(OUTPUT).c
	$(CC) $(OUTPUT).c -o $(OUTPUT) $(CC_OPT) $(CC_PTH) 

.PHONY: clean
clean:
	rm $(OUTPUT)

