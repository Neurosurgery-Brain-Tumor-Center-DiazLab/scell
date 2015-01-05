# Aaron Clauset
# created Sept 2001

# Makefile  02.22.2003

#### local macros
# remove without fussing about it
RM = /bin/rm -f

# compiler name and flags
CCC = g++
CCFLAGS = -O3 -fomit-frame-pointer -funroll-loops -fforce-mem -fforce-addr -fexpensive-optimizations

# loader flags
LDFLAGS = 

### local program information
EXEC=fitHRG
SOURCES=fitHRG.cc


### intermediate objects
OBJECTS = $(SOURCES: .cc=.o)

### includes
INCLUDES = 

### headers
HEADERS = 

### targets, dependencies and actions
$(EXEC): $(OBJECTS) Makefile
	$(LINK.cc) -o $(EXEC) $(OBJECTS)

### sort out dependencies
depend:
	makedepend $(INCLUDES) $(HEADERS) $(SOURCES)

### housekeeping

clean:
	$(RM) *.o *~

cleanall: clean
	$(RM) $(EXEC)
