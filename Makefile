TARGETS = genmat gencrs gencrs_new gendia rk rk_new rk_new_last geninit vis
OBJECTS = $(TARGETS:=.o) graphlib.o

#OMP_NUM_THREADS ?= 16

#export OMP_NUM_THREADS = 1
#export OMP_PLACES="{0,1,2,3,4,5,6,11,16,17,18,19,20,21,22,23}"
#export OMP_PLACES="{7,8,9,10,12,13,14,15}"

#export OMP_PROC_BIND=spread

CC = g++
CCFLAGS = -g -O2 -lm -ldl -Wall -Wextra  -fopenmp -lpthread -fsanitize=address,undefined

.PHONY:	all
all:	$(TARGETS)

%:	%.c
	$(CC) $(CCFLAGS) $^ -o $@

%.o:	%.c
	$(CC) $(CCFLAGS) $^ -c -o $@

.PHONY:	clean
clean:
	-$(RM) $(TARGETS)
	-$(RM) $(OBJECTS)

.PHONY:	realclean
realclean:	clean
	-$(RM) *~
