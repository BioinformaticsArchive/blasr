Set up directories
  $ CURDIR=$TESTDIR
  $ REMOTEDIR=/mnt/secondary-siv/secondarytest/testdata/BlasrTestData/ctest
  $ DATDIR=$REMOTEDIR/data
  $ OUTDIR=$CURDIR/out
  $ STDDIR=$REMOTEDIR/stdout

Set up the executable: loadPulses.
  $ BIN=$TESTDIR/../pbihdfutils/bin
  $ EXEC=$BIN/loadPulses

#Test loadPulses: input is a bas.h5 file
#Test -byread and -bymetric
  $ BAS_IN_1=$DATDIR/lambda_lp.fofn
  $ CMP_IN_1=$DATDIR/lambda_lp.cmp.h5
  $ CMP_STDOUT_1=$STDDIR/lambda_lp.cmp.h5
  $ CMP_OUT_byread_1=$OUTDIR/lambda.byread.cmp.h5
  $ CMP_OUT_bymetric_1=$OUTDIR/lambda.bymetric.cmp.h5
  $ METRICS=QualityValue,MergeQV,InsertionQV,DeletionQV,DeletionTag,PulseWidth,SubstitutionQV,SubstitutionTag

  $ rm -rf $CMP_OUT_byread_1
  $ cp $CMP_IN_1 $CMP_OUT_byread_1
  $ $EXEC $BAS_IN_1 $CMP_OUT_byread_1 -metrics $METRICS -byread > $OUTDIR/tmp.log
  $ h5diff $CMP_OUT_byread_1 $CMP_STDOUT_1
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

  $ rm -rf $CMP_OUT_bymetric_1
  $ cp $CMP_IN_1 $CMP_OUT_bymetric_1
  $ $EXEC $BAS_IN_1 $CMP_OUT_bymetric_1 -metrics $METRICS -bymetric > $OUTDIR/tmp.log 
  $ h5diff $CMP_OUT_bymetric_1 $CMP_STDOUT_1 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

#Test loadPulses: input is a pls.h5 file
#Test -byread and -bymetric
  $ PLS_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_2=$DATDIR/ecoli_lp_tiny.cmp.h5
  $ CMP_STDOUT_2=$STDDIR/ecoli_lp_tiny.cmp.h5
  $ CMP_OUT_byread_2=$OUTDIR/ecoli_lp_tiny.byread.cmp.h5
  $ CMP_OUT_bymetric_2=$OUTDIR/ecoli_lp_tiny.bymetric.cmp.h5
  $ METRICS=StartFrame,PulseWidth,WidthInFrames,pkmid,IPD,Light

  $ rm -rf  $CMP_OUT_byread_2
  $ cp $CMP_IN_2 $CMP_OUT_byread_2
  $ $EXEC $PLS_IN $CMP_OUT_byread_2 -metrics $METRICS -byread
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2

  $ h5diff $CMP_OUT_byread_2 $CMP_STDOUT_2
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

  $ rm -rf $CMP_OUT_bymetric_2
  $ cp $CMP_IN_2 $CMP_OUT_bymetric_2
  $ $EXEC $PLS_IN $CMP_OUT_bymetric_2 -metrics $METRICS -bymetric 
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2

  $ h5diff $CMP_OUT_bymetric_2 $CMP_STDOUT_2 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

#Test loadPulses for deep sorted cmp.h5
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_SORTED=$DATDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_STDOUT_SORTED=$STDDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_OUT_SORTED_bymetric=$OUTDIR/ecoli_lp_tiny_sorted_bymetric.cmp.h5
  $ CMP_OUT_SORTED_byread=$OUTDIR/ecoli_lp_tiny_sorted_byread.cmp.h5
  $ METRICS=StartFrame,PulseWidth,WidthInFrames,pkmid,IPD,Light,DeletionQV,InsertionQV,SubstitutionQV,MergeQV,QualityValue,DeletionTag,SubstitutionTag,ClassifierQV,PreBaseFrames,PulseIndex

  $ rm -rf $CMP_OUT_SORTED_bymetric
  $ cp $CMP_IN_SORTED $CMP_OUT_SORTED_bymetric
  $ $EXEC $FOFN_IN $CMP_OUT_SORTED_bymetric -bymetric -metrics $METRICS > $OUTDIR/tmp.log
  $ h5diff $CMP_OUT_SORTED_bymetric $CMP_STDOUT_SORTED 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]


  $ rm -rf  $CMP_OUT_SORTED_byread
  $ cp $CMP_IN_SORTED $CMP_OUT_SORTED_byread
  $ $EXEC $FOFN_IN $CMP_OUT_SORTED_byread -byread -metrics $METRICS > $OUTDIR/tmp.log

  $ h5diff $CMP_OUT_SORTED_bymetric $CMP_STDOUT_SORTED 
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d+ differences found (re)
  [1]

#Test loadPulses for a zero-alignment cmp.h5 file.
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN_NOALN=$DATDIR/noaln_lp.cmp.h5
  $ $EXEC $FOFN_IN $CMP_IN_NOALN -byread -metrics $METRICS
  WARNING, there is no alignment in the cmp file.


#Test loadPulses -byMetric with a 'large' bas.h5 file of which the dataset size is greater than maxElements. 
  $ FOFN_IN=$DATDIR/ecoli_lp.fofn
  $ CMP_IN=$DATDIR/ecoli_lp_tiny_sorted.cmp.h5
  $ CMP_OUT=$OUTDIR/ecoli_lp_maxEle.cmp.h5
  $ CMP_STDOUT=$STDDIR/ecoli_lp_maxEle.cmp.h5
  $ METRICS=QualityValue,MergeQV,InsertionQV,DeletionQV,DeletionTag,PulseWidth,SubstitutionQV,SubstitutionTag
  $ MAX_ELEMENTS=140000000

  $ rm -rf $CMP_OUT
  $ cp $CMP_IN $CMP_OUT
  $ $EXEC $FOFN_IN $CMP_OUT -bymetric -metrics $METRICS -maxElements $MAX_ELEMENTS
  Loading pulses from .+ by read. (re)
  loading 2 alignments for movie 1
  loading 2 alignments for movie 2

  $ h5diff $CMP_OUT $CMP_STDOUT
  dataset: </FileLog/CommandLine> and </FileLog/CommandLine>
  \d+ differences found (re)
  dataset: </FileLog/Timestamp> and </FileLog/Timestamp>
  \d differences found (re)
  [1]




