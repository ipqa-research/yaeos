
WD=/home/ruther/docs/research/thermo/algos/yaeos
all:
	fpm install --prefix "capi" --flag "-fPIC -ftrampoline-impl=heap"
	gfortran \
		-shared -fPIC -ftrampoline-impl=heap -g -Wall \
		-I. -I$(WD)/capi/include -L$(WD)/capi/lib \
		-o libycffi.so yaeos_call_c.f90 -lyaeos
#	fpm install --prefix "capi"
#	gfortran \
#		-fPIC -ftrampoline-impl=heap -g -Wall \
#		-I. -I$(WD)/capi/include -L$(WD)/capi/lib \
#		-o test_capi yaeos_call_c.f90 ar_c.c -lyaeos
#	./test_capi
