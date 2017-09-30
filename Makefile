CC=gcc
CFLAGS=-std=gnu11 -Wall -Wextra -g -O3 -Isrc -DNDEBUG -DGSL_RANGE_CHECK_OFF -DHAVE_INLINE $(OPTFLAGS)
#CFLAGS=-std=gnu11 -Wall -Wextra -g -O3 -ffast-math -mtune=native -Isrc -DNDEBUG $(OPTFLAGS)
#CFLAGS= -std=gnu99 -Wall -Wextra -pedantic -g -O2 -Isrc -DNDEBUG -fomit-frame-pointer -mtune=native
LIBS= -lm -lgsl -lgslcblas -lgmp $(OPTLIBS)
PREFIX?=/usr/local
SRCDIR=src
OBJDIR=obj

MAIN_SOUCES:=$(wildcard $(SRCDIR)/*_main.c)
MAINS:=$(subst _,-,$(patsubst %_main.c,%,$(notdir $(MAIN_SOUCES))))

SOURCES:=$(filter-out $(MAIN_SOUCES), $(wildcard $(SRCDIR)/*.c))
OBJECTS:=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(patsubst %.c,%.o,$(SOURCES)))


TEST_SRC:=$(wildcard tests/*_tests.c)
TESTS:=$(patsubst %.c,%,$(TEST_SRC))

TARGET=build/liblattice.a
#SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

# The Target Build
all: $(TARGET) $(MAINS)

dev: CFLAGS=-std=gnu11 -Wall -Wextra -g -Isrc $(OPTFLAGS)
dev: all

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(OBJECTS)
	ar rcs $@ $(OBJECTS)
	ranlib $@

#$(SO_TARGET): $(TARGET) $(OBJECTS)
#	$(CC) -shared -o $@ $(OBJECTS)

an-solve: $(SRCDIR)/an_solve_main.c $(TARGET)
	$(CC) $(CFLAGS) $^ -o an-solve $(LIBS)
	
lat-gen: $(SRCDIR)/lat_gen_main.c $(TARGET)
	$(CC) $(CFLAGS) $^ -o lat-gen $(LIBS)

lat-solve: $(SRCDIR)/lat_solve_main.c $(TARGET)
	$(CC) $(CFLAGS) $^ -o lat-solve $(LIBS)

lat-sim: $(SRCDIR)/lat_sim_main.c $(TARGET)
	$(CC) $(CFLAGS) $^ -o lat-sim $(LIBS)

rnd-point: $(SRCDIR)/rnd_point_main.c $(TARGET)
	$(CC) $(CFLAGS) $^ -o rnd-point $(LIBS)


$(OBJDIR)/%.o: $(SRCDIR)/%.c | build
	$(CC) $(CFLAGS) -c $< -o $@

#@mkdir -p bin
build:
	@mkdir -p build
	@mkdir -p obj

# The Unit Tests
.PHONY: tests
tests: 
tests: CFLAGS += -Wno-unused-parameter
tests: $(TESTS) 
	sh ./tests/runtests.sh

%_tests: %_tests.c $(TARGET)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

valgrind:
	VALGRIND="valgrind --log-file=tests/valgrind.log" $(MAKE) tests

# The Cleaner
.PHONY: clean
clean:
	rm -rf build $(OBJDIR) $(TESTS)
	rm -f tests/tests.log
	find . -name "*.gc*" -exec rm {} \;
	rm -rf `find . -name "*.dSYM" -print`

# The Install
install: all
	install -d $(DESTDIR)/$(PREFIX)/lib/
	install $(TARGET) $(DESTDIR)/$(PREFIX)/lib/

# The Checker
BADFUNCS='[^_.>a-zA-Z0-9](str(n?cpy|n?cat|xfrm|n?dup|str|pbrk|tok|_)|stpn?cpy|a?sn?printf|byte_)'
check:
	@echo Files with potentially dangerous functions.
	@egrep $(BADFUNCS) $(SOURCES) || true

objects:
	@echo $(OBJECTS)

sources:
	@echo $(SOURCES)

main-sources:
	@echo $(MAIN_SOUCES)

mains:
	@echo $(MAINS)

