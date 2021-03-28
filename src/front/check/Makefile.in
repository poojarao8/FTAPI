.PHONY: clean

INCS=-I../../../include -g

test_fapi: test_fapi.c ../../../include/FTAPI.h
	gcc $(INCS) test_fapi.c  -lcheck -o test_fapi -L../../../lib/x86_64 -lFronTier -lm -lmpi
clean:
	rm test_fapi
