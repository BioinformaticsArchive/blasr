SHELL = /bin/bash
JABON = ../build/jabon 

vanilla:
	@$(JABON) 						\
		--query query.fasta 		\
		--target target.fasta 		\
		> $@.out 2> $@.err
	
CHAIN_PARMS = --windowSize 20 --windowStep 1 --seed 8 --minChainLength 2

highdiag:
	@$(JABON)						\
		--diagonalCutoff 10	 		\
		--joinCutoff     30			\
		--query query_diag.fasta 	\
		--target target_diag.fasta 	\
		$(CHAIN_PARMS)	 			\
		> $@.out 2> $@.err

lowdiag:
	@$(JABON) 						\
		--diagonalCutoff 0	 		\
		--joinCutoff     30			\
		--query query_diag.fasta 	\
		--target target_diag.fasta 	\
		$(CHAIN_PARMS)	 			\
		> $@.out 2> $@.err

highjoin:
	@$(JABON) 						\
		--diagonalCutoff 10 		\
		--joinCutoff 30 			\
		--query query_diag.fasta 		\
		--target target_diag.fasta 		\
		$(CHAIN_PARMS)	 			\
		> $@.out 2> $@.err

lowjoin:
	@$(JABON) 						\
		--diagonalCutoff 10 		\
		--joinCutoff 0 				\
		--query query.fasta 		\
		--target target.fasta 		\
		$(CHAIN_PARMS)	 			\
		> $@.out 2> $@.err

TESTS = vanilla 	\
		highdiag 	\
		lowdiag 	\
		highjoin 	\
		lowjoin

test: $(TESTS)
	@for OUTFILE in `echo *.out` ; { diff $$OUTFILE $$OUTFILE.std ; }	
	@echo All tests ok

standards: $(TESTS)
	@for OUTFILE in `echo *.out` ; { \cp -f $$OUTFILE $$OUTFILE.std ; }	
	@echo Making standards complete

clean:
	rm -f *.out
	rm -f *.err

