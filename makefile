# Builds the distributable version

build-package: 
	python setup.py sdist -k

clean: 
	rm -r PathoVar-0.4.1
	rm -r build/

