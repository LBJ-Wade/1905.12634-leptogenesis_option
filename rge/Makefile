VENDOR_DIR=$(shell git rev-parse --show-toplevel)/vendor
MR=$(VENDOR_DIR)/mr/build

rgecheck: rgecheck.cpp
	$(CXX) -o$@ -I$(MR)/include/mr -L$(MR)/lib -lmr -ltsil $<

rgecheck_for_paper: rgecheck_for_paper.cpp
	$(CXX) -o$@ -I$(MR)/include/mr -L$(MR)/lib -lmr -ltsil $<
