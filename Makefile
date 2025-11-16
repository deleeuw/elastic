

ESRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSElasticEngine.c smacofSSElasticMajorize.c smacofSSElasticMonotone.c

%.o: %.c smacofSSElastic.h
	clang -c $@
	
eshlib: smacofSSElastic.h $(ESRC)
	R CMD SHLIB -o smacofSSElastic.so $(ESRC)

clean:
	rm -f *.o *.so

