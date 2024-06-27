python:
	fpm install --profile release --prefix build/python
	f2py -m yaeos \
		-L$(PWD)/../build/python/lib/ \
		-I$(PWD)/../build/python/include \
		-c yaeos_c.f90 -lyaeos
	mv yaeos.* yaeos/
