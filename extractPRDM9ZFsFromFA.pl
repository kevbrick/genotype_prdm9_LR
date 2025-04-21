#/usr/bin/perl
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Align::PairwiseStatistics;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Math::Round qw (round);
use List::UtilsBy qw(max_by);
use List::Util qw(max);

use Getopt::Long;

GetOptions ('fa=s'   => \(my $seq_input_fa),
            'id=s'   => \(my $sample_id),
            'tmp=s'  => \(my $tmpdir),
            'type=s' => \(my $type),
            'n=i'    => \(my $n_seqs_required = 30),
            'ulim=i' => \(my $upper_limit = 2000),
            'aln=s'  => \(my $aligner = 'blast'),
            'ZFs=s'  => \(my $known_ZFs = '/usr/local/genotypePRDM9/Alleva_et_al_2021_File_S2_PRDM9_ZF_details.txt'),
            'se+'    => \(my $detsToStdErr),
            'so+'    => \(my $simple_outputs),
            'v+'     => \(my $verbose));

unless ($tmpdir){
	if ($ENV{'SLURM_JOBID'} && -e '/lscratch/'.$ENV{'SLURM_JOBID'}){
		$tmpdir = '/lscratch/'.$ENV{"SLURM_JOBID"}.'/';
	}else{
		$tmpdir = $ENV{TMPDIR}?$ENV{TMPDIR}:'';
	}
}

die("FASTA ($seq_input_fa) does not exist                             !") unless (-e $seq_input_fa);
die("Known ZFs file ($known_ZFs) does not exist                       !") unless (-e $known_ZFs);
die("Sample name required (--id)!")                                       unless ($sample_id);
die("Temp folder not found (specify \$TMPDIR or --tmp)")                  unless ($tmpdir);

#count sequences in fasta
my $seq_count = 0;
open(my $fh, '<', $seq_input_fa) or die "Cannot open $seq_input_fa: $!";
while(<$fh>) {
    $seq_count++ if /^>/;
}
close($fh);

if ($seq_count < $n_seqs_required){
  die("**ERROR** - FASTA INPUT FILE - $n_seqs_required sequences required for ZF inference, but only $seq_count found in $seq_input_fa\n");
}

## Make random name stem
my $rand_temp_id = int(rand()*100000000000);
my $tmpinit = $tmpdir;
$tmpdir = $tmpdir."/$rand_temp_id";

my $verbose_outputs = $simple_outputs?0:1;

system("mkdir -p $tmpdir");

my $rand_name_base                         = "zfMatch_".int(rand()*1000000000000000);

## FOR MATCHER
my ($gap_open_penalty,$gap_extend_penalty) = (20,2);

## make BioPerl object to load raw ZF array $fasta
my $obj_fa = Bio::SeqIO->new(-format => 'fasta',
                             -file   => $seq_input_fa);

my $base_name = $tmpinit."/$sample_id.PrZFA".($type?".$type":"");

## For reporting purposes
$type = $type?$type:"longreads";

if ($detsToStdErr){
  open OUTDETS, '>&', STDERR;
}else{
  open OUTDETS, '>', $base_name.".details.txt" if ($verbose_outputs);
}

open READ2ZFASTAT, '>', $base_name.".stats.txt" if ($verbose_outputs);

## Keep stats of how we're calling ZFAs vs total reads
print READ2ZFASTAT join("\t",'total_ZFAs','usable_ZFAs','seq_type','zfa_size','count','pc')."\n" if ($verbose_outputs);

my ($totSeqs, $noZFSeqs, $goodSeqs, $fasta_seq_count);
my (%zf_array_count_by_size, $zf_array_count_total, %two_alleles);
my (%zf_array_sizes_to_use, %zf_aln, $zf_array_sizes_known);

my $alignment_output_count = 0;

## This loop will find the ZFs in each sequence, check that the flanks are found, and build a
## diploid consensus sequence for each. We loop over each sequence in the fasta
while ( my $obj_seq = $obj_fa->next_seq ) {

	## Unnecessary precaution, but I'll leave it here.
	last if ($goodSeqs > $upper_limit);

	$totSeqs++;

	## Triage reads with no ZFs
	my ($oneZF, $rc, $seq2use) = is_there_at_least_one_ZF($obj_seq->seq, $obj_seq->revcom->seq);

	unless ($oneZF){$noZFSeqs++; next}

	## Calls the get_zf_array function. This will find sequences that have a
	## complete zf array and return details of the array
	my ($number_of_zfs, $low_scoring_ZF, $seq_of_entire_zf_array, @zfs) = get_ZF_array($seq2use, $obj_seq->id);

    ## if we have no ZFs, then skip
	next unless ($number_of_zfs);

	## We're going to loop through the ZFs, so set counter to zero
	my $current_zf_pos = 0;

	## This keeps track of the number of zf arrays by size
	$zf_array_count_by_size{$number_of_zfs}++;
    print OUTDETS "ZF COUNT: $number_of_zfs\n" if ($verbose_outputs);

    $zf_array_count_total++;
    $goodSeqs++;
  

    for my $n (sort {$a <=> $b} keys(%zf_array_count_by_size)){
		print READ2ZFASTAT join("\t",$totSeqs,$goodSeqs,$type,$n,$zf_array_count_by_size{$n},sprintf("%4.2f",$zf_array_count_by_size{$n}/$goodSeqs))."\n" if ($verbose_outputs);
	}

	## Try to determine the likely sizes of the two alleles.
	$zf_array_sizes_known = check_zf_array_sizes($zf_array_count_total, \%zf_array_count_by_size, \%two_alleles, \%zf_array_sizes_to_use);

	## If we know the zf sizes, then stop processing zf arrays that are the wrong size
	next unless ($zf_array_sizes_to_use{$number_of_zfs} || not ($zf_array_sizes_known));

	## Next, what we do is keep track of the ZF sequences found at each position in the array
	## We will continue until we reach a reasonable $consensus

	## Create bioperl alignment if necessary
	$zf_aln{$number_of_zfs} = Bio::SimpleAlign->new() unless ($zf_aln{$number_of_zfs});

	my $seq = Bio::LocatableSeq->new(-seq   => $seq_of_entire_zf_array,
	                                 -id    => $obj_seq->primary_id.":::zfa:".$number_of_zfs.":".$zf_array_count_total,
	                                 -start => 1,
	                                 -end   => length($seq_of_entire_zf_array));

	$zf_aln{$number_of_zfs}->add_seq($seq);

  if ($zf_array_sizes_known){
  	print OUTDETS "Array Size Estimate: ZF array sizes are OK [".join(",",keys %zf_array_sizes_to_use)."]... \n" if ($verbose_outputs);
    if ($goodSeqs >= $n_seqs_required && !($goodSeqs % 10)){
      for my $sz(keys %zf_array_sizes_to_use){
		$alignment_output_count++;
        output_fa_from_aln(\$zf_aln{$sz},"$base_name.$sz\ZFs.fa");
      }
    }
  }else{
  	print OUTDETS "Array Size Estimate: ZF Array size ambiguity ... continuing ... \n" if ($verbose_outputs);
  }
}

close OUTDETS if ($verbose_outputs);
close READ2ZFASTAT if ($verbose_outputs);

if ($goodSeqs >= $n_seqs_required && $zf_array_sizes_known){
  for my $sz(keys %zf_array_sizes_to_use){
	$alignment_output_count++;
    output_fa_from_aln(\$zf_aln{$sz},"$base_name.$sz\ZFs.fa");
  }
}

if ($alignment_output_count == 0){
  die("No ZF arrays found in $seq_input_fa\n");
}

system("rm -rf $tmpdir");

###############################################################################
sub is_there_at_least_one_ZF{
	my ($inSeq, $inSeqRC) = @_;

	my $tmpBase  = $tmpdir."/onechk_$rand_name_base";
	makeFASTAs("$tmpBase");

  ## Check both FWD and REVCOMP sequences
	## Keep the one with the best ZF match
	my $inFasta  = "$tmpBase.input.fa";
	open FA, '>', $inFasta;
	print FA ">prdm9\n$inSeq";
	close FA;

	my ($sSeq,$aScore,$aFrom,$aTo)         = get_position_on_read("$tmpBase.ZFA.fa",$inFasta);

	my $inFastaRC  = "$tmpBase.RC.input.fa";
	open FARC, '>', $inFastaRC;
	print FARC ">prdm9\n$inSeqRC";
	close FARC;

	my ($sSeqRC,$aScoreRC,$aFromRC,$aToRC) = get_position_on_read("$tmpBase.ZFA.fa",$inFastaRC);

  ## If neither aligns with a score > 60 (eyeballed this) then there is no ZF
	## Otherwise, keep the better of the two
	if ($aScore > $aScoreRC){
		return($aScore,0,$inSeq) if ($aScore > 60 && ($aTo-$aFrom) > 65);
	}else{
		return($aScoreRC,1,$inSeqRC) if ($aScoreRC > 60 && ($aToRC-$aFromRC) > 65);
	}

	return (0,0,"");
}

###############################################################################
sub get_ZF_array{
	my ($inSeq,$rname) = @_;

  ## Make tep file names
	my $tmpBase  = $tmpdir."/$rand_name_base";
	my $outWater = "$tmpBase.water";
  my $outBlast = "$tmpBase.blast";
	my $inFasta  = "$tmpBase.input.fa";
	my $inFastaZ = "$tmpBase.justZFs.fa";
  my $inFasta2 = "$tmpBase.ZFs.fa";

	open FA, '>', $inFasta;
	print FA ">prdm9\n$inSeq";
	close FA;

	## This simply makes the LHS, RHS and ZF sequences
  my @ZFFastas = makeFASTAs("$tmpBase");

	##Get the location of the LHS and RHS
	my ($lSeq,$lScore,$lFrom,$lTo) = get_position_on_read("$tmpBase.LHS.fa",$inFasta);
	my ($rSeq,$rScore,$rFrom,$rTo) = get_position_on_read("$tmpBase.RHS.fa",$inFasta);

	##Skip this array unless we find both the left and right flanking sequences
	return unless($lTo && $rFrom && ($lTo < $rFrom));

  ## Snip out the likely ZF array and make a fasta file
	open FA2, '>', $inFastaZ;
	print FA2 ">prdm9_zfs\n".substr($inSeq,$lTo-10,($rFrom-$lTo)+10)."\n";
	close FA2;

	## Guess number of ZFs
	my $estimate_ZF_count = round(($rFrom-$lTo)/84);

  my (%zfs,%zfDets);

  if ($aligner eq 'matcher'){
    ## Slow & Accurate

    for my $fastaZF(@ZFFastas){
  	  ## Align ZF sequence to ZF array (allow a max of the expected number of matches)
      system("matcher -alt $estimate_ZF_count  -gapopen $gap_open_penalty -gapextend $gap_extend_penalty -asequence $fastaZF -bsequence $inFastaZ -outfile $outWater -aformat pair 2>/dev/null");
      getZFsFromMatcherAlignment($outWater,\%zfs,\%zfDets);
    }
  }

  if ($aligner eq 'blast'){
    ## FAST & Accurate (10-20x faster)
    system("makeblastdb -dbtype nucl -in $inFastaZ >/dev/null 2>/dev/null");
    system("blastn -query $tmpBase.ZFPUB.fa -db $inFastaZ -word_size 7 -max_hsps 200 -num_alignments 20000 -evalue 1 -culling_limit 20000 -outfmt 0 >$outBlast");

    getZFsFromBlastAlignment($outBlast,\%zfs,\%zfDets);
  }

	my (@outDetail, %zfOK,%numnum,@ZFnumz,$zfs_found, $zfarray_seq, $num_zfs);
	my ($kPrev,$previous_zf_end,$previous_zf_num,$contiguous_ZFs,$contiguous_nums) = (0,0,0,"","");
	my ($lowScoreZF,@zfSequences);

	## Now loop through ZFs in positional order
	foreach my $k (sort { $zfs{$a}->{from} <=> $zfs{$b}->{from} || $zfs{$b}->{score} <=> $zfs{$a}->{score} } keys %zfs) {

		## Guess which ZF this is (not used @ the moment, except for output info)
		my $nZF    = $previous_zf_num?$previous_zf_num + round(($zfs{$k}->{to}-($previous_zf_end))/84):1;
		my $zf_gap = round(($zfs{$k}->{from}-($previous_zf_end)));

		## Skip this ZF if there's already a better quality quality alignment to it
		## Should rarely happen @ this point because we only allow the number of
		## matches to be equal to the number of ZFs
		my $skipZF;
		$skipZF++ if ($zfOK{$nZF} && $zfOK{$nZF}->{score} >= $zfs{$k}->{score});

		## Build new hash of ZFs in order
		unless ($skipZF){
      if ($zfOK{$nZF} ){
        my $xx = 1;
      }
			$zfOK{$nZF}->{seq}     = $zfs{$k}->{seq};
			$zfOK{$nZF}->{zf}      = $nZF;
			$zfOK{$nZF}->{zfNoGap} = $zfs_found++;
			$zfOK{$nZF}->{gapsize} = $zf_gap;
			$zfOK{$nZF}->{zfaFrom} = $zfs{$k}->{from}-$zfDets{minPos};
			$zfOK{$nZF}->{zfaTo}   = $zfs{$k}->{to}-$zfDets{minPos};
			$zfOK{$nZF}->{from}    = $zfs{$k}->{from};
			$zfOK{$nZF}->{to}      = $zfs{$k}->{to};
			$zfOK{$nZF}->{score}   = $zfs{$k}->{score};
			$zfOK{$nZF}->{skip}    = $zfs{$k}->{skip};
			$zfOK{$nZF}->{lOK}     = $lTo?1:0;
			$zfOK{$nZF}->{rOK}     = $rFrom?1:0;

      $previous_zf_end = $zfOK{$nZF}->{to}+1;
      $previous_zf_num = $nZF;
    }
  }

  ## Now loop through ZFs in positional order
  foreach my $nZF (sort { $a <=> $b } keys %zfOK) {
		push @ZFnumz, $nZF unless ($numnum{$nZF}++);

		$zfarray_seq .= $zfs{$nZF}->{seq};

		print STDERR join("\t",$zfOK{$nZF}->{zf},
		                      $zfOK{$nZF}->{zfNoGap},
    											$zfOK{$nZF}->{gapsize},
		      								$zfOK{$nZF}->{from},
					      					$zfOK{$nZF}->{to},
								      		$zfOK{$nZF}->{score},
										      $zfOK{$nZF}->{seq},
    											$lTo,
		      								$rFrom,
					      					sprintf("%4.2f",($rFrom-$lTo)/84))."\n" if ($verbose);

		push @outDetail, join("\t",$fasta_seq_count,
		                           $zfOK{$nZF}->{seq},
		                           $nZF,
		                           $zfOK{$nZF}->{zfaFrom},
		                           $zfOK{$nZF}->{zfaTo},
		                           $lTo." -- ".$rFrom,
		                           $zfOK{$nZF}->{from},
		                           $zfOK{$nZF}->{to},
		                           $zfOK{$nZF}->{score},
		                           $zfOK{$nZF}->{skip}?"SKIP":($zfOK{$nZF}->{score}<200?"LOSCORE":"OK"));

		$lowScoreZF++ if ($zfOK{$nZF}->{score}<200 && $nZF > 1);

		push @zfSequences, $zfOK{$nZF}->{seq};

		$num_zfs++;
	}

  print ("-" x 60)."\n\n" if ($verbose);

  if ($estimate_ZF_count == $num_zfs){
    print OUTDETS "GOOD ESTIMATE: $estimate_ZF_count = $num_zfs\n" if ($verbose_outputs);
  }else{
    print OUTDETS "MEH ESTIMATE: $estimate_ZF_count = $num_zfs\n" if ($verbose_outputs);
  }

	for my $o(@outDetail){print OUTDETS $o."\t".join("\t",$totSeqs,$noZFSeqs,$goodSeqs,"SO")."\n" if ($verbose_outputs)} ;# if ($allZFchk < 999 && $rChk)}

  ## Make return sequence ...
	## 2bp left flank -- ZFs -- 2bp right flank
	## 2bp overhangs are to allow pentamer inference later
	#$zfarray_seq = substr($inSeq,$lTo-2,2).join("",@zfSequences).substr($inSeq,$rFrom,2);

  $zfarray_seq = join("",@zfSequences);

  ## Swap gaps for Xs
	$zfarray_seq =~ s/\-/X/g;

	return($num_zfs, $lowScoreZF, $zfarray_seq, @zfSequences);

}

###############################################################################
sub getZFsFromMatcherAlignment {
  my ($waterAlignmentFile,$h_zfs,$h_zfDets) = @_;

  my ($outSeq,$skip,$padleft,$padright,$nAln);
  my ($arrLOK, $arrROK, $zfarray_seq, $num_zfs);
  my ($minPos,$maxPos) = (1e9,-1e9);

  my $f;

  if (keys(%{$h_zfs})){
    $f = max(keys(%{$h_zfs}));
  }else{
    $f = 0;
  }

  #print STDERR "Water-v-$fastaZF ...\n";
  my @okArray;
  # Import ZF alignment to bioperl obj
  my $in = Bio::AlignIO->new(-format => 'emboss',
                             -file   => $waterAlignmentFile);

  ## Loop through the individual ZF alignments
  while ( my $aln = $in->next_aln ) {

    my (@qry,@tgt,%zf);
    my ($skip,$padleft,$padright,$nAln) = (0,"","",0);

    $f++;

    $zf{$f}->{skip}  = 0;

    my $alnStat = Bio::Align::PairwiseStatistics->new();;

    # Get start, end, and sequence of each zf in query
    my $target = $aln->get_seq_by_pos(1);
    my $query  = $aln->get_seq_by_pos(2);

    next if ($target->seq !~ /^[GATCNX\-]+$/);
    next if ($query->seq !~ /^[GATCNX\-]+$/);
    next if ($query->length < 50);
    next unless ($alnStat->score_nuc($aln));

    $minPos = ($minPos < $f)?$minPos:$f;
    $maxPos = ($maxPos > $f)?$maxPos:$f;

    my @qry = split("",uc($query->seq));
    my @tgt = split("",uc($target->seq));

    ## PAD ZF with N's if we don't quite match at the ends
    if ($target->start != 1){
      $padleft  = "N" x ($target->start-1) ;
    }
    if ($target->end != 84){
      $padright = "N" x (84-$target->end) ;
    }

    ## Build zf-array sequence with left/right padding
    my $mySeq = $padleft;

    for my $i(0..$#qry){
      $qry[$i] = lc($qry[$i]) if ($tgt[$i+1] eq "-");
      $qry[$i] = lc($qry[$i]) if ($tgt[$i-1] eq "-");
      $qry[$i] = lc($qry[$i]) if ($tgt[$i] eq "-");
      $mySeq .= $qry[$i] if ($tgt[$i] ne "-");
    }

    $mySeq = $mySeq.$padright;

    ## PAD N's around gaps
    # $mySeq =~ s/[GATCN\-]{2}([gatcn\-]+)[GATCN\-]{2}/NN$1NN/g;
    # $mySeq =~ s/[GATCN\-]{2}([gatcn\-]+)/NN$1/g;
    # $mySeq =~ s/([gatcn\-]+)[GATCNgatcn\-]{2}/$1NN/g;
    # $mySeq =~ s/[GATCNgatcn\-]([gatcn\-]+)/N$1/g;
    # $mySeq =~ s/([gatcn\-]+)[GATCNgatcn\-]/$1N/g;
    $mySeq =~ s/[n\-]/N/gi;

    ## Account for some parser error
    if ($mySeq =~ /\//){
      $zf{$f}->{skip}  = 1;
    }

    next if ($zf{$f}->{skip} );
    $$h_zfs{$f}->{from}  = $query->start;
    $$h_zfs{$f}->{to}    = $query->end;
    $$h_zfs{$f}->{score} = $alnStat->score_nuc($aln);
    $$h_zfs{$f}->{zfstart} = $target->start;
    $$h_zfs{$f}->{zfend}   = $target->end;
    $$h_zfs{$f}->{seq}   = $mySeq;
    $$h_zfs{$f}->{skip}  = $skip;
    $$h_zfs{$f}->{query} = $query->seq;
    $$h_zfs{$f}->{target} = $target->seq;
    $$h_zfDets{minPos}   = $minPos;
    $$h_zfDets{maxPos}   = $maxPos;
  }
}

###############################################################################
sub getZFsFromBlastAlignment {
  my ($blastAlignmentFile,$h_zfs,$h_zfDets) = @_;

  my ($outSeq,$skip,$padleft,$padright,$nAln);
  my ($arrLOK, $arrROK, $zfarray_seq, $num_zfs);
  my ($minPos,$maxPos) = (1e9,-1e9);

  my $f;
  if (keys(%{$h_zfs})){
    $f = max(keys(%{$h_zfs}));
  }else{
    $f = 0;
  }

  #print STDERR "Water-v-$fastaZF ...\n";
  my @okArray;
  # Import ZF alignment to bioperl obj
  my $in = Bio::SearchIO->new(-format => 'blast',
                             -file   => $blastAlignmentFile);

  ## Loop through the individual ZF alignments
  while ( my $result = $in->next_result ) {
    while ( my $hit = $result->next_hit() ) {
      while( my $hsp = $hit->next_hsp() ) {
          #print join("\t", $hit->name, $hsp->{HIT_START}, $hsp->{HIT_END}, $hsp->bits) . "\n";

          my (@qry,@tgt,%zf);
          ($skip,$padleft,$padright,$nAln) = (0,"","",0);

          $f++;

          my $zfSkip  = 0;

          next if ($hsp->{HIT_SEQ} !~ /^[GATCNX\-]+$/);
          next if ($hsp->{QUERY_SEQ} !~ /^[GATCNX\-]+$/);
          next if ($hsp->{HIT_LENGTH} < 50);

          next unless ($hsp->bits);

          $minPos = ($minPos < $f)?$minPos:$f;
          $maxPos = ($maxPos > $f)?$maxPos:$f;

          my @hit = split("",uc($hsp->{HIT_SEQ}));
          my @zf  = split("",uc($hsp->{QUERY_SEQ}));

          ## PAD ZF with N's if we don't quite match at the ends
          if ($hsp->{QUERY_START} != 1){
            $padleft  = "N" x ($hsp->{QUERY_START}-1) ;
          }

          if ($hsp->{QUERY_END} != 84){
            $padright = "N" x (84-$hsp->{QUERY_END}) ;
          }

          ## Build zf-array sequence with left/right padding
          my $mySeq = $padleft;

          for my $i(0..$#hit){
            $hit[$i] = lc($hit[$i]) if ($zf[$i+1] eq "-");
            $hit[$i] = lc($hit[$i]) if ($zf[$i-1] eq "-");
            $hit[$i] = lc($hit[$i]) if ($zf[$i] eq "-");
            $mySeq .= $hit[$i] if ($zf[$i] ne "-");
          }

          $mySeq = $mySeq.$padright;

          ## PAD N's around gaps
          # $mySeq =~ s/[GATCN\-]{2}([gatcn\-]+)[GATCN\-]{2}/NN$1NN/g;
          # $mySeq =~ s/[GATCN\-]{2}([gatcn\-]+)/NN$1/g;
          # $mySeq =~ s/([gatcn\-]+)[GATCNgatcn\-]{2}/$1NN/g;
          # $mySeq =~ s/[GATCNgatcn\-]([gatcn\-]+)/N$1/g;
          # $mySeq =~ s/([gatcn\-]+)[GATCNgatcn\-]/$1N/g;
          $mySeq =~ s/[n\-]/N/gi;

          ## Account for some parser error
          if ($mySeq =~ /\//){
            $zfSkip  = 1;
          }

          next if ($zfSkip);
          $$h_zfs{$f}->{from}  = $hsp->{HIT_START};
          $$h_zfs{$f}->{to}    = $hsp->{HIT_END};
          $$h_zfs{$f}->{score} = $hsp->bits;
          $$h_zfs{$f}->{zfstart} = $hsp->{QUERY_START};
          $$h_zfs{$f}->{zfend}   = $hsp->{QUERY_START};
          $$h_zfs{$f}->{seq}   = $mySeq;
          $$h_zfs{$f}->{skip}  = $skip;
          $$h_zfs{$f}->{query} = $hsp->{HIT_SEQ};
          $$h_zfs{$f}->{target} = $hsp->{QUERY_SEQ};
          $$h_zfDets{minPos}   = $minPos;
          $$h_zfDets{maxPos}   = $maxPos;
        }
      }
    }
  }

###############################################################################
sub to84{
	my $inLine = shift;
	my $lnRet;
	while ($inLine =~ s/^(.{84})//){
		$lnRet .= $1."\n";
	}
	$lnRet .= $inLine."\n" if ($inLine);
	return($lnRet);
}

###############################################################################
sub get_position_on_read{
	my ($checkFA,$iFA) = @_;

	my $tf = $tmpdir."/$rand_name_base.out";

	## Do Smith-Waterman (LALIGN) alignment
	system("matcher -alt 30 -gapopen $gap_open_penalty -gapextend $gap_extend_penalty -asequence $checkFA -bsequence $iFA -outfile $tf -aformat pair 2>/dev/null");

	## Import aligment to bioperl object
	my $in = Bio::AlignIO->new(-format => 'emboss',
                             -file   => $tf);

	my ($iSeq,$iScore,$iFrom,$iTo);

	while ( my $aln = $in->next_aln ) {
		foreach my $seq ($aln->each_seq) {
			if ($seq->id eq 'prdm9'){
				$iFrom = $seq->start;
				$iTo = $seq->end;
			}
			return(length($seq->seq),$aln->average_percentage_identity,$iFrom,$iTo) if ($iFrom);
		}
	}

	return(0,0,0,0);

}

###############################################################################
sub check_zf_array_sizes{
	my ($tot,$sz,$two_X_alleles,$ok_sizes,$zf_known) = @_;

  my $consensus_threshold = 0.7;

	return 0 if ($tot < $n_seqs_required);

	my @zf_array_rank = sort { $$sz{$b} <=> $$sz{$a} } keys(%{$sz});

  ## Only one possible allele size so far
  if ($#zf_array_rank == 0){
    $$ok_sizes{$zf_array_rank[0]}++;
    $$two_X_alleles{$zf_array_rank[0]}++;
    return 1
  }
  my $ratio12 = ($$sz{$zf_array_rank[0]}/$$sz{$zf_array_rank[1]});
	my $ratio13 = ($$sz{$zf_array_rank[0]}/($zf_array_rank[2]?$$sz{$zf_array_rank[2]}:1));
	my $ratio23 = ($$sz{$zf_array_rank[1]}/($zf_array_rank[2]?$$sz{$zf_array_rank[2]}:1));

	## If the top two ZFs represent N% of all tested ccs
	if (($$sz{$zf_array_rank[0]} + $$sz{$zf_array_rank[1]}) > $tot*$consensus_threshold){
		if ($ratio12 < 2){
			$$ok_sizes{$zf_array_rank[0]}++;
			$$ok_sizes{$zf_array_rank[1]}++;
			return 1;
		}

    if ($ratio12 < 3 && $ratio23 > 2){
      $$ok_sizes{$zf_array_rank[0]}++;
			$$ok_sizes{$zf_array_rank[1]}++;
			return 1;
		}

		if ($ratio12 > 3){
			$$ok_sizes{$zf_array_rank[0]}++;
			$$two_X_alleles{$zf_array_rank[0]}++;
			return 1;
		}

	}

	if ($$sz{$zf_array_rank[0]} > $tot*$consensus_threshold){
		$$ok_sizes{$zf_array_rank[0]}++;
		$$two_X_alleles{$zf_array_rank[0]}++;
		return 1;
	}

	## If we got here, we still don't know
	return 0;
}

###############################################################################
sub makeFASTAs{
	my $stub = shift;

	open  LHSFA, '>', $stub.".LHS.fa";
	print LHSFA ">lhs\n";
	print LHSFA "CACAGCCGTAATGACAAAACCAAAGGTCAAGAGATCAAAGAAAGGTCCAAACTCTTGAATAAAAGGACATGGCAGAGGGAGATT\n";
	print LHSFA "TCAAGGGCCTTTTCTAGCCCACCCAAAGGACAAATGGGGAGCTGTAGAGTGGGAAAAAGAATAATGGAAGAAGAGTCCAGAACA\n";
	print LHSFA "GGCCAGAAAGTGAATCCAGGGAACACAGGCAAATTATTTGTGGGGGTAGGAATCTCAAGAATTGCAAAAGTCAAGTATGGAGAG\n";
	close LHSFA;

	open  RHSFA, '>', $stub.".RHS.fa";
	print RHSFA ">rhs\n";
	print RHSFA "GATGAGTAAGTCATTAGTAATAAAACCTCATCTCAATAGCCACAAAAAGACAAATGTGGTCACCACACACTTGCACACCCCAGC\n";
	print RHSFA "TGTGAGGTGGCTTCAGCGGAAGTCTGCTGACCCCTTATATTCCCCGAGAGTATAAAGAGATCGGAAATAACTGATTAAACAAAT\n";
	print RHSFA "CCGCCACTTTCATGACTAGAGATGAGGAAGAACAAGGGATAGTTCTGTAAGTGTTCGGGGGACATCAGCATGTGTGGTTCTTTC\n";
	close RHSFA;

	open ZFPUB, '>', $stub.".ZFPUB.fa";
  open IN, $known_ZFs;
  my ($ZFCount, @publishedZFFA);
  while (<IN>){
    chomp;
		next if ($_ =~ /^(#|code)/);
    my @F = split(/\t/,$_);
    my $ZFFastaFile = $stub.".".(++$ZFCount).".pub.fa";
    open ZFOUT, '>', $ZFFastaFile;
    print ZFOUT ">".join("_",$F[0..2])."\n$F[3]\n";
    print ZFPUB ">".join("_",$F[0..2])."\n$F[3]\n";
    close ZFOUT;
    push @publishedZFFA, $ZFFastaFile;
  }
  close IN;
  close ZFPUB;

	open  ZFAFA, '>', $stub.".ZFA.fa";
	print ZFAFA ">zfa\n";
	print ZFAFA "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAG\n";
	close ZFAFA;

	open  ZFFA, '>', $stub.".ZF.fa";
	print ZFFA ">zf\n";
	print ZFFA "TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG\n";
	close ZFFA;

  return(@publishedZFFA);
}

###############################################################################
sub output_fa_from_aln {
	my ($objAln, $output_fa_name) = @_;

	my $seq_count = 0;
	open MA_OUT, '>', $output_fa_name;
	for my $s($$objAln->each_seq){
		$seq_count++;
		print MA_OUT ">".($s->display_id?$s->display_id.":ZFAcount=$seq_count":$seq_count)."\n".to84($s->seq);
	}
	close MA_OUT;
}

##################################################################################
sub sort_hash_by_values{
		my ($h,$e) = @_;

		my %ret_h;
		my $hCount = 0;

    ## Get the total number of calls at this locus
		for my $hv(keys (%{$h})){
			$ret_h{tot} += $$h{$hv};
		}

		#Order calls by frequency
		for my $hv(sort {$h->{$b} <=> $h->{$a}} keys (%{$h})){
			## Treat "N" different ... we don't really want Ns to be ranked
			if ($hv eq 'N'){
				$ret_h{$hv}->{key} = $hv;
				$ret_h{$hv}->{N}   = $$h{$hv};
				$ret_h{$hv}->{pc}  = $$h{$hv}/$ret_h{tot};
			}else{
				$ret_h{$hCount}->{key} = $hv;
				$ret_h{$hCount}->{N}   = $$h{$hv};
				$ret_h{$hCount}->{pc}  = $$h{$hv}/($ret_h{tot}*($e?$e:1));
				$hCount++;
			}
		}

    # Store hom/het/neither status
		$ret_h{status} = 0;
		$ret_h{homhet} = 'unk';
		# 40x coverage is a minimum
		if ($ret_h{tot} > 50){
			## Are we looking at pentamers? If so, we have slightly different criteria
			if (length($ret_h{0}->{key}) == 5){

				my $n_diff = ( $ret_h{0}->{key} ^ $ret_h{1}->{key} ) =~ tr/\0//c;

				## OK. So, if the pentamers differ only by the middle base, it's a
				## possible het. Treat it differently.
				if ($n_diff == 1 && (substr($ret_h{0}->{key},2,1) ne substr($ret_h{1}->{key},2,1))) {
					if ($ret_h{0}->{pc} <= $ret_h{1}->{pc}*1.6 ){
							$ret_h{status} = 2;
							$ret_h{homhet} = 'het';
							#$ret_h{keyval} = $ret_h{0}->{key};
							#$ret_h{pcval}  = $ret_h{0}->{pc};
							$ret_h{keyval} = join("/",$ret_h{0}->{key},$ret_h{1}->{key});
							$ret_h{pcval} = join("/",$ret_h{0}->{pc},$ret_h{1}->{pc});
					}else{
						if ($ret_h{tot} > 80 &&
						   ($ret_h{1}->{pc} >= $ret_h{2}->{pc}*4) &&
						    $ret_h{0}->{pc} <= $ret_h{1}->{pc}*3) {
							$ret_h{status} = 2;
							$ret_h{homhet} = 'het';
							#$ret_h{keyval} = $ret_h{0}->{key};
							#$ret_h{pcval}  = $ret_h{0}->{pc};
							$ret_h{keyval} = join("/",$ret_h{0}->{key},$ret_h{1}->{key});
							$ret_h{pcval} = join("/",$ret_h{0}->{pc},$ret_h{1}->{pc});
						}else{
							$ret_h{status} = 1;
					    $ret_h{homhet} = 'hom';
							$ret_h{keyval} = $ret_h{0}->{key};
							$ret_h{pcval}  = $ret_h{0}->{pc};
						}
					}
				}else{
					if ($ret_h{0}->{pc} >= $ret_h{1}->{pc}*1.4 ){
							$ret_h{status} = 1;
							$ret_h{homhet} = 'hom';
							$ret_h{keyval} = $ret_h{0}->{key};
							$ret_h{pcval}  = $ret_h{0}->{pc};
					}
				}
			}else{
				## NOT PENTAMERS

				if ($ret_h{0}->{pc} >= 0.5){
					if ($ret_h{0}->{pc} > $ret_h{1}->{pc}*2 ){
						$ret_h{status} = 1;
					  $ret_h{homhet} = 'hom';
						$ret_h{keyval} = $ret_h{0}->{key};
						$ret_h{pcval}  = $ret_h{0}->{pc};
					}
				}

				if ($ret_h{tot} > 150 &&
					 ($ret_h{1}->{pc} >= $ret_h{2}->{pc}*4) &&
						$ret_h{0}->{pc} <= $ret_h{1}->{pc}*2 ) {
					$ret_h{status} = 2;
					$ret_h{homhet} = 'het';
					#$ret_h{keyval} = $ret_h{0}->{key};
					#$ret_h{pcval}  = $ret_h{0}->{pc};
					$ret_h{keyval} = join("/",$ret_h{0}->{key},$ret_h{1}->{key});
					$ret_h{pcval} = join("/",$ret_h{0}->{pc},$ret_h{1}->{pc});
				}
			}
		}

    # if ($hCount > 1){
		# 	$ret_h{ratio12} = abs(log($ret_h{0}->{pc}/$ret_h{1}->{pc})/log(2));
		#
		# 	if ($ret_h{tot} > 150){
		# 		if ($ret_h{0}->{pc} >= 0.25){
		# 			if ($ret_h{ratio12} < (log(1.2)/log(2))) {
		# 				$ret_h{status} = 2;
		# 				$ret_h{keyval} = join("/",$ret_h{0}->{key},$ret_h{1}->{key});
		# 				$ret_h{pcval} = join("/",$ret_h{0}->{pc},$ret_h{1}->{pc});
		# 			}
		# 		}
		# 	}
		# }

		sub orNA{
			my $v = shift;
			return ($v?$v:"NA");
		}

		my $RRRRRR = 1;

		$ret_h{println} =  join("\t",orNA($ret_h{keyval}),
																 orNA($ret_h{tot}),
																 orNA($ret_h{status}),
																 orNA(sprintf('%4.2f',$ret_h{pcval})),
																 orNA(sprintf('%4.2f',$ret_h{0}->{pc})),
																 orNA(sprintf('%4.2f',$ret_h{1}->{pc})),
																 orNA(sprintf('%4.2f',$ret_h{2}->{pc})));

		$h = \%ret_h;
}
