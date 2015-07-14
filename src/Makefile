libdir=$(prefix)/lib

CPPFLAGS=-Wall -g -fPIC
CPP=g++

RPIHUANG2004_SRC=rpihuang2004.C
RPIHUANG2004_OBJ=$(RPIHUANG2004_SRC:.C=.o)

build: libRpiHuang2004.a libRpiHuang2004.so

libRpiHuang2004.a: $(RPIHUANG2004_OBJ)
	ar rc $@ $^

libRpiHuang2004.so: $(RPIHUANG2004_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: $(libdir)/libRpiHuang2004.a $(libdir)/libRpiHuang2004.so

$(libdir)/libRpiHuang2004.a $(libdir)/libRpiHuang2004.so: libRpiHuang2004.a libRpiHuang2004.so
	cp libRpiHuang2004.a libRpiHuang2004.so $(libdir)

uninstall:
	- rm -f $(libdir)/libRpiHuang2004.a $(libdir)/libRpiHuang2004.so

clean:
	- rm -f libRpiHuang2004.a libRpiHuang2004.so $(RPIHUANG2004_OBJ)
