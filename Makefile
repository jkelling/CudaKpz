CPATH += inc/
export KEPLER = on

SUBDIRS = inc tools kpz

all: $(SUBDIRS)

inc:
	$(MAKE) -C $@ allall

tools: inc
	$(MAKE) -C $@

kpz: inc
	$(MAKE) -C $@ kpzScaling

clean:
	$(MAKE) -C inc clean
	$(MAKE) -C kpz clean

.PHONY: all $(SUBDIRS)
