CXX = g++
CXXFLAGS = -std=c++11 -O3
ifeq ($(CXX), icpc)
  CXXFLAGS += -qopenmp
else
  CXXFLAGS += -fopenmp
endif
INCLUDES = -I$(PWD)/include

LDFLAGS = -lm -ldl



SRCDIR := src
SRCS = $(sort $(wildcard $(SRCDIR)/*.cpp))
OBJSDIR := .objs
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJSDIR)/%.o, $(SRCS))
TARGET = main
DEPSDIR := .deps
MAKES := Makefile

DEPGEN = -MT $@ -MMD -MP -MF $(DEPSDIR)/$*.d

.PHONY: all clean distclean
all: $(TARGET)

$(OBJS): $(OBJSDIR)/%.o: $(SRCDIR)/%.cpp $(DEPSDIR)/%.d $(MAKES) | $(DEPSDIR) $(OBJSDIR)
	$(CXX) $(DEPGEN) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

$(DEPSDIR):
	@mkdir -p $@
$(OBJSDIR):
	@mkdir -p $@

DEPSINCS = $(patsubst $(SRCDIR)/%.cpp, $(DEPSDIR)/%.d, $(SRCS))
$(DEPSINCS):

include $(sort $(wildcard $(DEPSINCS)))

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS) -fopenmp

clean:
	-rm -rvf $(OBJSDIR)
	-rm -rvf $(DEPSDIR)
	-rm -f $(TARGET)

distclean: clean
