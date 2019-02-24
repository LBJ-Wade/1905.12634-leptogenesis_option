.PHONY: relax init freeze

#PYTHON2=`pyenv versions | grep 2\.\[67\]\. | sed 's/.*\(2\.[67]\.[0-9]*\).*/\1/'`
#PYTHIA8=/usr/local/pythia8
#PYTHIA8_VERSION=8235
#ROOTLIB=${ROOTSYS}/lib

relax:
	@echo usage: make init

init: requirements.txt ../.gitmodules
	#pip install -r requirements.txt
	git submodule init
	git submodule update
	# build
	./mr_build.sh

freeze:
	mv -f requirements.txt requirements.bak
	pip freeze -r requirements.bak > requirements.txt
