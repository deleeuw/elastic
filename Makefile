
SSRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSEngine.c smacofSSMajorize.c smacofSSMonotone.c

ESRC = smacofIsotone.c smacofMPInverseV.c smacofSort.c smacofSSElasticEngine.c smacofSSMajorize.c smacofSSMonotone.c

%.o: %.c smacofSS.h
	clang -c $@

sshlib: smacofSS.h $(SSRC)
	R CMD SHLIB -o smacofSS.so $(SSRC)

eshlib: smacofSSElastic.h $(ESRC)
	R CMD SHLIB -o smacofSSElastic.so $(ESRC)

clean:
	rm -f *.o *.so

