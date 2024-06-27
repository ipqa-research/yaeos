python:
	fpm install --flag "-g" --profile debug --prefix build/python
	f2py -m yaeos \
		-L$(PWD)/../build/python/lib/ \
		-I$(PWD)/../build/python/include \
		-c yaeos_c.f90 -lyaeos # --debug-capi
	mv yaeos.* yaeos/
