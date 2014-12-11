# Builds the distributable version

build:
	python setup.py build

build-package: 
	python setup.py sdist -k

clean: 
	rm -rf PathoVar-*
	rm -rf PathoVar.egg*
	rm -rf build/
	rm -rf dist/
