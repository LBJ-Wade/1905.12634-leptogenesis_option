VENDOR_DIR=$(shell git rev-parse --show-toplevel)/vendor
MR=$(VENDOR_DIR)/mr/build

rgecheck: rgecheck.cpp
	$(CXX) -o$@ -I$(MR)/include/mr -L$(MR)/lib -lmr -ltsil $<
