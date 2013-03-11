#!/usr/bin/env perl

package AlignmentParser;


sub ParseCompareSequencesLine {
		my ($line, $alignment) = @_;
		@items = split(/\s+/, $line);
#		x1_y40_m091023_093702_Uni_p2_b20 352 3 355 +  gi|9626243|ref|NC_001416.1| 293 5263811 5264104 + -1048 278 7 8 67 TAAGAGTGTTAATGGTCAATT-GTCGAAACTAAAAAAGCTG-CTGTTTTTTTTTGAAAGTTTTTTTTTTTTTAATTTGAAAAAATTTTGAATTTTTGGCCAAAAATCTTTGGCACAATTGGCTCAAAAATAAATCAACAAAAGCAACTTTGGTTTGGCTCTGGCAAGGAAAAAGGTAAAAAATAGC-AGTCGTTAAGC-TTTTTGAATAATATGTTTTTTCTGGGTGG-GTAAATTTCTAATTGGTTAAAATAAAATGTTTGCCTTTTTTTTTTAAATAAATTTTTTTTTTTGTGAAGAAGC-TTTCATTTGCAAGCCCACAAC-T-CAGGTTTTCCACCGACATGGGTTCCAGGGATTT ||||||||||||||||||*||*|||||****||||||||||*|||*******|||**||*******||||**||||||||||||*||||||*|||||**|**|||||*||||||*||||*||||**||||||||||||||||||||||||||||**|||||*|||||||||||*||**||||||*|*|||||||||||*||||||*||*|||||***|||||***|||*||||||||||||||*||**||||||||||*||||****||||||*|||||||********|||*|||||||||*|||||*||||*||||*|||||*|*||||*|||||||||||||**|||*||||||||| TAAGAGTGTTAATGGTCATTTGGTCGA----AAAAAAGCTGCCTG-------TTG--AG-------TTTTAGAATTTGAAAAAA-TTTGAAATTTTG--C--AAATC-TTGGCA-AATT-GCTC--AAATAAATCAACAAAAGCAACTTTGGTT--GCTCT-GCAAGGAAAAA-GT--AAAATA-CAAGTCGTTAAGCTTTTTTG-AT-ATATG---TTTCT---TGGTGTAAATTTCTAATT-GT--AAATAAAATG-TTGC----TTTTTTAAAATAAA--------TTT-TGAAGAAGCTTTTCAGTTGC-AGCCGACAACTTCCAGG-TTTCCACCGACAT--GTT-CAGGGATTT

		${$alignment}{"qTitle"}       = @items[0];
		${$alignment}{"qAlignLength"} = @items[1];
		${$alignment}{"qAlignStart"}  = @items[2];
		${$alignment}{"qAlignEnd"}    = @items[3];
		${$alignment}{"qStrand"}      = @items[4];
		${$alignment}{"tTitle"}       = @items[5];
		${$alignment}{"tAlignLength"} = @items[6];
		${$alignment}{"tAlignStart"}  = @items[7];
		${$alignment}{"tAlignEnd"}    = @items[8];
		${$alignment}{"tStrand"}      = @items[9];
		${$alignment}{"score"}        = @items[10];
		${$alignment}{"nMatch"}       = @items[11];
		${$alignment}{"nMismatch"}    = @items[12];
		${$alignment}{"nIns"}         = @items[13];
		${$alignment}{"nDel"}         = @items[14];
		${$alignment}{"qString"}      = @items[15];
		${$alignment}{"tString"}      = @items[17];
}


return 1;


