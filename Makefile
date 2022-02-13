#
# TODO: Move `libmongoclient.a` to /usr/local/lib so this can work on production servers
#
 
CC := cc # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
FF := gfortran
SRCDIR := src
BUILDDIR := build
TARGET := bin/metbp.linux
LINKCC := cc
 
SRCEXT := c
FOREXT := f
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) ) 
SOURCESF := $(shell find $(SRCDIR) -type f -name *.$(FOREXT) )
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))   $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCESF:.$(FOREXT)=.o))
CFLAGS := -g -Wall -static
LIB := -lm -lgfortran -lquadmath #-pthread -lmongoclient -L lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt
# INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<
	
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(FOREXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(FF)  -c -o $@ $<"; $(FF)  -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR)/*.o $(TARGET)"; $(RM) -r $(BUILDDIR)/*.o $(TARGET)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.c $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean



