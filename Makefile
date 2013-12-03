# Makefile

CFLAGS=-c -Wall
C_PATH=c
PYTHON_PATH=python
DATA_PATH=data
TARGET_PATH=target
VITERBI=viterbi4
DATA_FILE=data.tsv

all: clean compile run

compile:
	mkdir $(TARGET_PATH)
	gcc -o $(TARGET_PATH)/$(VITERBI) $(C_PATH)/$(VITERBI).c
	
run:
	python $(PYTHON_PATH)/simulation.py $(TARGET_PATH)/$(VITERBI) $(DATA_PATH)/$(DATA_FILE)

data:
	python $(PYTHON_PATH)/generate_data.py $(DATA_PATH)/$(DATA_FILE)
	
clean:
	rm -rf $(TARGET_PATH)