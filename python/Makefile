python:
	fpm install --profile release --prefix build/python
	f2py -m yaeos \
		-L ../build/python/lib/libyaeos.a \
		-I ../build/python/include \
		-c yaeos_c.f90
	mv yaeos.* yaeos/
