Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/secondarytest/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: blasr.
  $ BIN=$TESTDIR/../alignment/bin
  $ EXEC=$BIN/blasr

Test blasr on ecoli.
Test blasr with -sam
#See $STDOUT/ecoli.sam for 1.4 output.
  $ rm -rf $OUTDIR/ecoli.sam
  $ $EXEC $DATDIR/ecoli.fasta $DATDIR/ecoli_reference.fasta -sam -out $OUTDIR/ecoli.sam -nproc 15
  $ tail -n+5 $OUTDIR/ecoli.sam | sort | cut -f 1-11| md5sum
  d2e2b6cfe710b7b6a065e73c3244b0b9  -

#  $ tail -n+5 $STDDIR/ecoli.sam | sort | cut -f 1-11| md5sum
#  d2e2b6cfe710b7b6a065e73c3244b0b9  -

Test blasr with -m 0 ~ 5 
  $ rm -rf $OUTDIR/read.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 0 -out $OUTDIR/read.m0
  $ diff $OUTDIR/read.m0 $STDDIR/read.m0

  $ rm -rf $OUTDIR/read.m1
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 1 -out $OUTDIR/read.m1
  $ diff $OUTDIR/read.m1 $STDDIR/read.m1

  $ rm -rf $OUTDIR/read.m2
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 2 -out $OUTDIR/read.m2
  $ diff $OUTDIR/read.m2 $STDDIR/read.m2

  $ rm -rf $OUTDIR/read.m3
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 3 -out $OUTDIR/read.m3
  $ diff $OUTDIR/read.m3 $STDDIR/read.m3

  $ rm -rf $OUTDIR/read.m4
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -m 4 -out $OUTDIR/read.m4
  $ diff $OUTDIR/read.m4 $STDDIR/read.m4

Test alignment score
  $ rm -rf $OUTDIR/testscore.m0
  $ $EXEC $DATDIR/read.fasta  $DATDIR/ref.fasta -minReadLength 1 -m 0 -out $OUTDIR/testscore.m0
  $ diff $OUTDIR/testscore.m0 $STDDIR/testscore.m0

Test affineAlign
  $ rm -rf $OUTDIR/affineAlign.m0
  $ $EXEC $DATDIR/affineAlign.fofn $DATDIR/substr_with_ins.fasta -m 0 -out $OUTDIR/affineAlign.m0  -affineAlign  -readIndex 493 -insertion 100 -deletion 100
  $ diff $OUTDIR/affineAlign.m0 $STDDIR/affineAlign.m0

