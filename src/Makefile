QFY = quantify 
PUR = purify 
CC0 = g++
CC1 = $(QFY) $(CC0)
CC2 = $(PUR) $(CC0)
CC = $(CC0)
LD = $(CC0)
FLAGS = -g -fsanitize=address # -O3
C_FLAGS = $(FLAGS) 
LD_FLAGS = $(FLAGS) 
FLAGS1 = #-DREP_COUNT -DREP_COUNT_PRINT # -DPROGRESS
FLAGS2 = -DIMPROVE_UPPER_BOUND -DIMPROVE_LAMBDA_BOUNDS
FLAGS3 = -DDEBUG # -DDFS_RECUR #-DPRINT_SCC
MYD_FLAGS = $(FLAGS1) $(FLAGS2) $(FLAGS3)
D_FLAGS = $(MYD_FLAGS) -Dlint -D__lint -Wall -Winline -Wno-deprecated -Wno-strict-overflow
D_FLAGS2 = $(D_FLAGS) -DADD_SOURCE_NODE

HEADERS = ad_globals.h ad_graph.h ad_queue.h ad_cqueue.h ad_pq.h ad_stack.h ad_util.h
OBJS0 = ad_main.o ad_util.o
OBJS1 = $(OBJS0) ad_graph.o
OBJS2 = $(OBJS1) ad_alg_util.o
OBJS3 = $(OBJS0) ad_graph2.o ad_pq.o

BURNS_BASE = burns
BURNS_PREFIXES = burn bur bu b
HOWARD_BASE = howard 
HOWARD_PREFIXES = howar howa how ho h
KARP_BASE = karp
KARP_PREFIXES = kar ka
KO_BASE = ko 
KO_PREFIXES = k
LAWLER_BASE = lawler
LAWLER_PREFIXES = lawle lawl law la l
SZY_BASE = szymanski 
SZY_PREFIXES = szymansk szymans szyman szyma szym szy sz s
TARJAN_BASE = tarjan
TARJAN_PREFIXES = tarja tarj tar ta t
VALITER_BASE = valiter
VALITER_PREFIXES = valite valit vali val va v
YTO_BASE = yto
YTO_PREFIXES = yt y

BURNS = ad_alg_burns
HOWARD = ad_alg_howard
KARP = ad_alg_karp
KO = ad_alg_ko
LAWLER = ad_alg_lawler
SZY = ad_alg_szymanski
TARJAN = ad_alg_tarjan
VALITER = ad_alg_valiter
YTO = ad_alg_yto

# Excluded Burns' algorithm & Karp's algorithm due to its slowness
EXES = $(HOWARD_BASE) \
  $(KO_BASE) \
  $(LAWLER_BASE) \
  $(SZY_BASE) \
  $(TARJAN_BASE) \
  $(VALITER_BASE) \
  $(YTO_BASE)

all: $(EXES)

# Burns' algorithm:
$(BURNS_PREFIXES): $(BURNS_BASE)

$(BURNS_BASE): $(OBJS1) $(BURNS).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS1) $(BURNS).o

$(BURNS).o: ad_queue.h ad_graph.h $(BURNS).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(BURNS).cc

# Howard algorithm:
$(HOWARD_PREFIXES) : $(HOW_BASE)

$(HOWARD_BASE): $(OBJS1) $(HOWARD).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS1) $(HOWARD).o

$(HOWARD).o: ad_queue.h ad_graph.h $(HOWARD).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(HOWARD).cc

# Karp's algorithm:
$(KARP_PREFIXES) : $(KARP_BASE)

$(KARP_BASE): $(OBJS1) $(KARP).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS1) $(KARP).o

$(KARP).o: ad_pq.h ad_graph.h $(KARP).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(KARP).cc

# Karp-Orlin algorithm:
$(KO_PREFIXES) : $(KO_BASE)

$(KO_BASE): $(OBJS3) $(KO).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS3) $(KO).o

$(KO).o: ad_pq.h ad_graph.h $(KO).cc
	$(CC) $(C_FLAGS) $(D_FLAGS2) -c $(KO).cc

# Lawler's algorithm:
$(LAWLER_PREFIXES) : $(LAW_BASE)

$(LAWLER_BASE): $(OBJS2) $(LAWLER).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS2) $(LAWLER).o

$(LAWLER).o: ad_cqueue.h ad_graph.h $(LAWLER).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(LAWLER).cc

# Szymanski's algorithm:
$(SZY_PREFIXES) : $(SZY_BASE)

$(SZY_BASE): $(OBJS2) $(SZY).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS2) $(SZY).o

$(SZY).o: ad_graph.h $(SZY).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(SZY).cc

# Tarjan's algorithm:
$(TARJAN_PREFIXES) : $(TARJAN_BASE)

$(TARJAN_BASE): $(OBJS2) $(TARJAN).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS2) $(TARJAN).o

$(TARJAN).o: ad_cqueue.h ad_graph.h $(TARJAN).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(TARJAN).cc

# Value iteration (valiter) algorithm
$(VALITER_PREFIXES) : $(VALITER_BASE)

$(VALITER_BASE): $(OBJS1) $(VALITER).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS1) $(VALITER).o

$(VALITER).o: ad_queue.h ad_graph.h $(VALITER).cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c $(VALITER).cc

# Young-Tarjan-Orlin algorithm:
$(YTO_PREFIXES) : $(YTO_BASE)

$(YTO_BASE): $(OBJS3) $(YTO).o
	$(LD) -o $@.x $(LD_FLAGS) $(OBJS3) $(YTO).o

$(YTO).o: ad_pq.h ad_graph.h $(YTO).cc
	$(CC) $(C_FLAGS) $(D_FLAGS2) -c $(YTO).cc

# Files required by all algorithms:
ad_main.o: ad_globals.h ad_util.h ad_graph.h ad_main.cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c ad_main.cc

ad_graph.o: ad_globals.h ad_util.h ad_graph.h ad_graph.cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c ad_graph.cc

ad_graph2.o: ad_globals.h ad_graph.h ad_graph.cc
	$(CC) $(C_FLAGS) $(D_FLAGS2) -c ad_graph.cc -o ad_graph2.o

ad_pq.o: ad_globals.h ad_pq.h ad_pq.cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c ad_pq.cc

ad_util.o: ad_globals.h ad_util.h ad_util.cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c ad_util.cc

ad_alg_util.o: ad_globals.h ad_graph.h ad_alg_util.cc
	$(CC) $(C_FLAGS) $(D_FLAGS) -c ad_alg_util.cc

# Testing:
test:
	./utest.sh

# Cleaning:
clean c cl cle clea: 
	rm -f *.o *~ core *.x
cleano: 
	rm -f *.o *~ core

# End of file
