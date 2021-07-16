# Get an EMBOSS factory
use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Math::Round qw (round);
use List::UtilsBy qw(max_by);
use List::Util qw(max);

use Getopt::Long;

GetOptions ('fa=s'     => \(my $input_fa),
	    			'c=s'      => \(my $consensus_threshold = 0.7),
						'onehap+'  => \(my $one_haplotype),
	    			'n=i'      => \(my $n_seqs_required = 50),
						'ep+'			 => \(my $useErrorProfiling),
						'so+'			 => \(my $simple_outputs),
						'pz=s'	   => \(my $published_prdm9_zfs = '/usr/local/genotypePRDM9/Alleva_et_al_2021_File_S2_PRDM9_ZF_details.txt'),
						'pa=s'	   => \(my $published_prdm9_alleles = '/usr/local/genotypePRDM9/Alleva_et_al_2021_File_S3_PRDM9_allele_details.txt'),
						'countN+'	 => \(my $count_N_values),
						'id=s'   	 => \(my $sample_id),
						'tmp=s'		 => \(my $tmpdir));

unless ($tmpdir){
	if ($ENV{'SLURM_JOBID'} && -e '/lscratch/'.$ENV{'SLURM_JOBID'}){
		$tmpdir = '/lscratch/'.$ENV{"SLURM_JOBID"}.'/';
	}else{
		$tmpdir = $ENV{TMPDIR}?$ENV{TMPDIR}:'';
	}
}

die("FASTA ($input_fa) does not exist!")     unless (-e $input_fa);
die("Sample ID required!")                   unless ($sample_id);
die("Invalid consensus threshold (0 - 1)")   if ($consensus_threshold < 0 || $consensus_threshold >1);
die("Cannot find R (load module)!")          unless (`which R`);
die("Temp folder not found (specify \$TMPDIR or --tmp)")      unless ($tmpdir);

my $rand_temp_id = int(rand()*100000000000);
my $tmpinit = $tmpdir;

$tmpdir = $tmpdir."/$rand_temp_id";
system("mkdir -p $tmpdir");

my $verbose_outputs = $simple_outputs?0:1;

my ($gap_open_penalty,$gap_extend_penalty) = (20,2);
my $rand_name_base                         = "zfMatch_".int(rand()*1000000000000000);
my $base_name = $tmpinit."/".$sample_id;

## Allow a lower consensus threshold if calling an obvious hom
$consensus_threshold = $one_haplotype?0.6:0.7;

## Create bioperl alignment if necessary
my $zfa_seqs = Bio::SeqIO->new(-format => 'Fasta', -file => $input_fa);

my $zf_aln = Bio::SimpleAlign->new();

my ($count, $array_length);

while ( my $s = $zfa_seqs->next_seq() ) {
	my $seq = Bio::LocatableSeq->new(-seq   => $s->seq,
																 	-id    => "zfa:".++$count,
																 	-start => 1,
																 	-end   => length($s->seq));

	$zf_aln->add_seq($seq);
	$array_length = length($s->seq)/84;
}

if ($verbose_outputs){
	open OUTZFS,  '>', $base_name.".$array_length.consensus_zf_seqs.raw";
	open OUTDETS, '>', $base_name.".$array_length.genotypeCalls.details.txt";
}

open OUTHAP,  '>', $base_name.".$array_length.haplotypes.txt";

my ($genotyping_problem,@genotypes) = try_to_guess_genotypes(\$zf_aln);

if ($genotyping_problem){
	if ($verbose_outputs){
		print OUTDETS "********************************************************\n";
		print OUTDETS "* Genotyping problem: $genotyping_problem\n";
		print OUTDETS "********************************************************\n";
	}
}else{
	for my $gt(@genotypes){
		my ($zf_array_sequence,$totZF);
		for my $nzf(sort {$a <=> $b} keys(%{$gt})) {
			print OUTZFS join("\t",$nzf,$gt->{$nzf})."\n" if ($verbose_outputs);
			$zf_array_sequence .= $gt->{$nzf};
			$totZF++;
		}
		my ($prdm9_allele_name, $zf_code) = compare2Published($gt);
		print OUTHAP join("\t",$sample_id,$prdm9_allele_name, $zf_code,$totZF,$zf_array_sequence)."\n";
	}
}
close OUTZFS;
close OUTHAP;
close OUTDETS;

system("rm -rf $tmpdir");

###############################################################################
sub to84{
  my $inLine = shift;
  my $lnRet;
  while ($inLine =~ s/^(.{84})//){
    $lnRet .= $1."\n";
  }
  $lnRet .= $1."\n";
  return($lnRet);
}

###############################################################################
sub getPosOnRead{
  my ($checkFA,$iFA) = @_;

  my $tf = $tmpdir."/$rand_name_base.out";
  ## ALIGN TO GET "A" ALLELE FIRST
  system("matcher -alt 30 -gapopen $gap_open_penalty -gapextend $gap_extend_penalty -asequence $checkFA -bsequence $iFA -outfile $tf -aformat pair 2>/dev/null");

  # Now you might want to get the alignment
  my $in = Bio::AlignIO->new(-format => 'emboss',
                             -file   => $tf);

  my ($iSeq,$iScore,$iFrom,$iTo);

  while ( my $aln = $in->next_aln ) {
      foreach my $seq ($aln->each_seq) {
        if ($seq->id eq 'prdm9'){
          $iFrom = $seq->start;
          $iTo = $seq->end;
        }
        return("",999,$iFrom,$iTo) if ($iFrom);
      }
  }
  return(0,0,0,0);
}

###############################################################################
sub do_we_have_a_consensus_diploid_sequence{
  my ($h_consensus, $h_sizes) = @_;

  return 0 unless (%{$h_consensus});

  my ($num_ok, $num_required);

  for my $sz(keys(%{$h_sizes})){
    my @zf_pos_with_consensus = keys(%{$$h_consensus{$sz}});
    $num_ok++ if (($#zf_pos_with_consensus+1) == $sz);
    $num_required++;
  }

  return ($num_ok && ($num_ok == $num_required))?1:0;
}

###############################################################################
sub try_to_guess_genotypes{

  my $aln = shift;

  my ($problems, @returnGenotypes);

  my $minimum_alleles_required = 1;
	my $maximum_alleles_required = $one_haplotype?1:2;

  my (@seq_array, @locSeqArray);

  for my $s($$aln->each_seq){
    push @seq_array, Bio::Seq->new(-seq => $s->seq,-id => $s->id);
    push @locSeqArray, $s;
  }

	my @heterozygosity_pos;
  my $pos = 0;

	my $array_size = round(length($seq_array[0]->seq)/84);

	## FYI: MSA=Multiple Sequence Alignment
	my $hom_msa_filename     = join("_",$base_name,"hom",$array_size,"msa.fa");
	my $hom_err_filename     = join("_",$base_name,"hom",$array_size,"errorprofile.txt");
	my $hom_scores_filename  = join("_",$base_name,"hom",$array_size,"consensus_scores.txt");

	my $het1_msa_filename    = join("_",$base_name,"hetallele1",$array_size,"msa.fa");
	my $het1_scores_filename = join("_",$base_name,"hetallele1",$array_size,"consensus_scores.txt");
	my $het2_msa_filename    = join("_",$base_name,"hetallele2",$array_size,"msa.fa");
	my $het2_scores_filename = join("_",$base_name,"hetallele2",$array_size,"consensus_scores.txt");

	system("rm -rf $hom_msa_filename $het1_msa_filename $het2_msa_filename");
	system("rm -rf $hom_scores_filename $het1_scores_filename $het2_scores_filename");
	system("rm -rf $hom_err_filename");

  ## Custom consensus function ... ignores Ns
  #for my $pc_id($aln->consensus_conservation()){
	my @zf_errors = build_zf_error_profile($aln,$hom_err_filename);
	output_fa_from_aln($aln,$hom_msa_filename) if ($verbose_outputs);

	my @consensus_hom = get_consensus($aln,\@zf_errors,'score');
	output_consensus_scores(\@consensus_hom,$hom_scores_filename) if ($verbose_outputs);

	for my $pc_id(@consensus_hom){
    if ($pc_id < $consensus_threshold*100){
      #print $pos."\n" ;
      push @heterozygosity_pos, $pos;
    }
    $pos++;
  }

  ## If there's no reason to think its a het ... return consensus.
  if (@heterozygosity_pos){
    my $short_aln = Bio::SimpleAlign->new();
    my (@seqArrayShort, @ss);
    my $s = 0;

		my %altGenotypes;
    for my $sq(@seq_array){
      my $selectSeq;
      $selectSeq .= $sq->subseq($_+1,$_+1) foreach(@heterozygosity_pos);
      push @ss, $selectSeq;

			$altGenotypes{$selectSeq}++;

      my $seq = Bio::LocatableSeq->new(-seq => "$selectSeq",
                     -id  => "selseq".$s++,
                     -alphabet => 'dna',
                     -start => 1,
                     -end   => length($selectSeq));

      $short_aln->add_seq($seq);

      my $sSeq = Bio::Seq->new(-seq => "$selectSeq",
                               -id  => "seq".$s++);
      push @seqArrayShort, $sSeq;
    }

    my $matrixFile = 'tst_'.int(rand()*1000000000000).'.csv';
    clusterMe($matrixFile,@ss);
		my $RS = makeRscript();
    system("Rscript $RS --distmatrix $matrixFile 2>/dev/null");

    my $clusterFile = $matrixFile.'.clusters.csv';
    open IN, $clusterFile;
    open my $CF, '<', $clusterFile;
    chomp(my @clusters = <$CF>);
    close $CF;

    my %clustAln;
    $clustAln{1} = Bio::SimpleAlign->new();
    $clustAln{2} = Bio::SimpleAlign->new();

    my $nSeq;
    for my $i(@clusters){
				## only consider top two clusters
			  if ($i == 1 || $i == 2){
        	$clustAln{$i}->add_seq($locSeqArray[$nSeq]);
				}
				$nSeq++;
    }

    my (%consensus, %cons1, %cons2, $ambiguity1, $ambiguity2);
		my $het_name = join("_",$sample_id,"hetallele1",$array_size);

    $ambiguity1 = splitToZFs(\$clustAln{1}, \%cons1, \@zf_errors, $het1_scores_filename);
    $ambiguity2 = splitToZFs(\$clustAln{2}, \%cons2, \@zf_errors, $het2_scores_filename);

    $problems .= ',ambiguity_in_haplotype1'  if ($ambiguity1);
    $problems .= ',ambiguity_in_haplotype2'  if ($ambiguity2);

		output_fa_from_aln(\$clustAln{1},$het1_msa_filename);
		output_fa_from_aln(\$clustAln{2},$het2_msa_filename);

    push @returnGenotypes, \%cons1;
    push @returnGenotypes, \%cons2;
  }else{
    my %consensus;
    splitToZFs($aln, \%consensus);

    push @returnGenotypes, \%consensus;
  }

	my $problems;
	$problems .= ',too_many_genotypes('.($#returnGenotypes+1).')' if ($#returnGenotypes > 1); ## Zero-based
	$problems .= ',too_few_genotypes('.($#returnGenotypes+1).')'  if (($#returnGenotypes+1) < $minimum_alleles_required );

	return($problems,@returnGenotypes);
}

###############################################################################
sub clusterMe{
  my ($matrix_csv,@s) = @_;
  open OUT, '>', $matrix_csv;
  my $matrix;
  for my $ii(0..$#s){
    for my $jj(0..$#s){
      my $dist = 0;
      if ($s[$ii] eq $s[$jj]){
        $matrix .= "0,";
      }else{
        my @seq1 = split(//,$s[$ii]);
        my @seq2 = split(//,$s[$jj]);
        for my $kk(0..$#seq1){
          $dist++ if ($seq1[$kk] ne $seq2[$kk]);
        }

        $matrix .= "$dist,";
      }
    }
    $matrix =~ s/[,\s]+$//;
    print OUT $matrix."\n";
    $matrix="";
  }

  close OUT;

  my $xsada=0
}

###############################################################################
sub splitToZFs{
  my ($oAln, $hRet, $zf_err_arr, $out_scores_file) = @_;
  my ($zfStart,$current_pos,$zfNum,$uncertain,$some_ambiguity) = (0,1,1,0,0);

  ## Custom consensus function ... ignores Ns
  #for my $pc_id($oAln->consensus_conservation()){
	my @consensus_arr = get_consensus($oAln,$zf_err_arr,'score');
	output_consensus_scores(\@consensus_arr,$out_scores_file);

  my $consensus_threshold_as_pc = $consensus_threshold*100;
	for my $pc_id(@consensus_arr){
      if ($pc_id < $consensus_threshold_as_pc){
        $uncertain = 1;
        $some_ambiguity = 1;
      }
      ## Each ZF
      unless (++$current_pos % 84){
				## NOTE: here we set the consensus threshold for the string as 0 because we've already decided that we meet the
				## necessary criteria for acceptance of this position. This is necessary for cases where there are many "N"s / "Xs"
				## that we expect to get from the error profile
				$$hRet{$zfNum++} = $uncertain?"U"x84:substr(get_consensus($oAln, $zf_err_arr,'seq'),$zfStart,84);
        ($zfStart,$uncertain) = ($current_pos,0);
      }
  }

  return($some_ambiguity);
  #$$hRet{$zfNum++} = $uncertain?"U"x84:substr($oAln->consensus_string(80),$zfStart,84);
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

	open  ZFAFA, '>', $stub.".ZFA.fa";
        print ZFAFA ">zfa\n";
        print ZFAFA "TGTGGACAAGGTTTCAGTGTTAAATCAGATGTTATTACACACCAAAGGACACATACAGGGGAGAAGCTCTACGTCTGCAGGGAG\n";
	close ZFAFA;

	open  ZFFA, '>', $stub.".ZF.fa";
        print ZFFA ">zf\n";
        print ZFFA "TGTGGGCGGGGCTTTAGCCGGCAGTCAGTCCTCCTCACTCACCAGAGGAGACACACAGGGGAGAAGCCCTATGTCTGCAGGGAG\n";
	close ZFFA;

}

###############################################################################
sub get_consensus{
	my ($o_aln, $error_profile, $retType) = @_;

	## Can either return a consensus score OR a consensus string
	$retType = 'score' unless ($retType);

	my (%scores, @ret_score, @ret_string);
	my $pos = 0;

	foreach my $o_seq ($$o_aln->each_seq) {
		$pos = 0;
		for my $nt(split(//,$o_seq->seq)){
			$scores{$pos++}->{$nt}++;
		}
	}

	for my $p(0..($pos-1)){
			my $tot;

			if ($count_N_values){ ## This is unnecessary, but keeping anyway
				$tot = $scores{$p}->{'A'} + $scores{$p}->{'C'} + $scores{$p}->{'G'} + $scores{$p}->{'T'} + $scores{$p}->{'N'} +
						 	 $scores{$p}->{'a'} + $scores{$p}->{'c'} + $scores{$p}->{'g'} + $scores{$p}->{'t'};
			}else{
				$tot = $scores{$p}->{'A'} + $scores{$p}->{'C'} + $scores{$p}->{'G'} + $scores{$p}->{'T'} +
						   $scores{$p}->{'a'} + $scores{$p}->{'c'} + $scores{$p}->{'g'} + $scores{$p}->{'t'};
			}

			## Get ZF position
			my $zf_position = ($p % 84);

			if ($zf_position == 59){
				my $debug = 1;
			}

			if ($p == 429){
				my $debug = 1;
			}

			## Blunt correction for anticipated @zf_errors
			## Based on ZF error profile
			$tot *= (1-$$error_profile[$zf_position]) if ($useErrorProfiling);

			my ($sc, $nt);
			if ($tot){
				if ($p == 1008 || $p == 1009){
					my $rr = 1;
				}
				$sc = max($scores{$p}->{'A'}/$tot,
									$scores{$p}->{'C'}/$tot,
									$scores{$p}->{'G'}/$tot,
									$scores{$p}->{'T'}/$tot,
									$scores{$p}->{'a'}/$tot,
									$scores{$p}->{'c'}/$tot,
									$scores{$p}->{'g'}/$tot,
									$scores{$p}->{'t'}/$tot)*100;

				for my $n('A','C','G','T','a','c','g','t'){
					$nt = $n if ($scores{$p}->{$n}/$tot*100 == $sc);
				}

				$sc = ($scores{$p}->{uc($nt)} + $scores{$p}->{lc($nt)})/$tot*100;
			}else{
				$sc = 0;
			}
			push @ret_score, $sc;
			push @ret_string, $nt;
	}

	return (@ret_score)         if ($retType eq 'score');
	return join("",@ret_string) if ($retType eq 'seq');
}

###############################################################################
sub output_fa_from_aln {
	my ($objAln, $output_fa_name) = @_;

	my $c = 0;
	open MA_OUT, '>', $output_fa_name;
	for my $s($$objAln->each_seq){
		print MA_OUT ">".(++$c)."\n".$s->seq."\n";
	}
	close MA_OUT;
}

###############################################################################
sub build_zf_error_profile{
	my ($objAln,$error_log_file) = @_;

	if ($verbose_outputs){
		open ZFERR, '>', $error_log_file ;
		print ZFERR "zf\tpos\terror\n";
	}

	my (%scores, @zf_error_profile);
	my $pos = 0;

	foreach my $o_seq ($$objAln->each_seq) {
		$pos = 0;
		for my $nt(split(//,$o_seq->seq)){
			$scores{$pos++}->{$nt}++;
		}
	}

	my $zf_pos = 0;
	my $num_zf = 1;

	for my $p(0..($pos-1)){

			my $tot = $scores{$p}->{'A'} + $scores{$p}->{'C'} + $scores{$p}->{'G'} + $scores{$p}->{'T'} + $scores{$p}->{'X'} + $scores{$p}->{'N'};
			my $err = $scores{$p}->{'X'} + $scores{$p}->{'N'};

			print ZFERR join("\t",$num_zf,$zf_pos,($err/$tot))."\n" if ($verbose_outputs);
			$zf_error_profile[$zf_pos] += $err/$tot;

			$zf_pos++;
			if ($zf_pos == 84){
				$zf_pos = 0 ;
				$num_zf++ ;
			}
	}

	for my $i(0..$#zf_error_profile){
		$zf_error_profile[$i]/=$num_zf;
	}

	close ZFERR if ($verbose_outputs);
	return (@zf_error_profile);
}

###############################################################################
sub output_consensus_scores{
	my ($cons,$log_file) = @_;

	open CONS, '>', $log_file;
	my $zf = 1 ;
	my $zfpos = 0;

	print CONS join("\t","zf","zfpos","pos","score")."\n";

	for my $ipos(0..$#{$cons}){
		print CONS join("\t",$zf,$zfpos,$ipos,$$cons[$ipos])."\n";
		if ($zfpos++ == 84){
			$zfpos = 0;
			$zf++;
		}
	}

	close CONS;
}

###############################################################################
sub compare2Published{

	my $zf_array = shift;

	## Parse ZFs first
	#open ZF, 'humanPRDM9ZFcodes.BergJeffreys.txt';
	open ZF, $published_prdm9_zfs;

	my (%zfs,%zfCounts);

	while (<ZF>){
		chomp;
		next if ($_ =~ /^(#|code)/);
		my @F = split(/\t/,$_);
		$zfs{$F[0]} = $F[3];
		$zfs{$F[3]} = $F[0];
	}

	close ZF;

	## Parse alleles second
	#open PRDM9, 'humanPRDM9alleles.BergJeffreys.txt';
	open PRDM9, $published_prdm9_alleles;

	my (%prdm9, %prdm9Counts);

	while (<PRDM9>){
		chomp;
		next if ($_ =~ /^(#|ID\s)/);
		my @F = split(/\t/,$_);
		$prdm9{$F[2]} = $F[0];
	}

	close PRDM9;

	## Now ...
	my $allele;
	for my $nZF (sort {$a <=> $b} keys(%{$zf_array})) {

		my $zfSeq = $$zf_array{$nZF};
		my $zfname = $zfs{$zfSeq}?$zfs{$zfSeq}:(($zfSeq =~ /U/)?"UU":"NU");

		$allele .= $zfname;
	}

	my $name =  $prdm9{$allele}?$prdm9{$allele}:(($allele =~ /UU/)?"Unk":"new");

	return($name, $allele);

}

###############################################################################
sub makeRscript{

	my $scriptName = $rand_name_base.".R";
	open RS, '>', $scriptName;
	print RS '
	#!/usr/bin/env Rscript
	library("optparse")
	library("tidyverse")

	option_list = list(
	  make_option(c("--distmatrix"),
	              type="character",
	              default=NULL,
	              help="Dissimilarity matrix csv",
	              metavar="character")
	)

	opt_parser = OptionParser(option_list=option_list);
	opt = parse_args(opt_parser);

	c_counts <- data.frame(numclusters=999,clust=999,Freq=999);

	#################################################################################################
	getOptimalClusterNumber <- function(hc,dm,max_clusters=20){
	  ## Vector to store cluster mean distances
	  means <- rep(999,max_clusters)

	  ## check each number of clusters and find optimum
	  ## Optimum = minimum within cluster mean distance
	  ##           (only consider two largest clusters; 1&2)
	  for (num_clusters in 1:max_clusters){
	    cdata                <- cutree(hc,num_clusters)

	    ## per cluster mean distances
	    means[num_clusters]  <- mean(dm[which(cdata == 1),which(cdata == 1)])

	    if (num_clusters > 1){
	      means[num_clusters] <- means[num_clusters] + mean(dm[which(cdata == 2),which(cdata == 2)])
	    }
	  }

	  optimal_clusters <- min(which(means==min(means)))
	  cdata            <- cutree(hc,optimal_clusters)

	  clusterOrder <- as.data.frame(table(cdata)) %>% arrange(desc(Freq));

	  badClusters <- as.numeric(clusterOrder$cdata[3:length(clusterOrder$cdata)])

	  ## Only keep the top two largest clusters
		repClusters <- rep(9,length(cdata))
    repClusters[cdata == clusterOrder$cdata[1]] <- 1
    repClusters[cdata == clusterOrder$cdata[2]] <- 2

	  #cdata[cdata %in% badCluster]         <- 9
	  #cdata[cdata == clusterOrder$cdata[1]] <- 1
	  #cdata[cdata == clusterOrder$cdata[2]] <- 2

	  return(repClusters)
	}

	##################################################################
	full_dist_matrix <- as.matrix(read.csv(opt$distmatrix,header=FALSE))
	pngOut           <- paste0(opt$distmatrix,".png")
	clustOut         <- paste0(opt$distmatrix,".clusters.csv")

	dist_matrix      <- full_dist_matrix
	dist_matrix[upper.tri(dist_matrix,diag=TRUE)] <- NA
	diag(dist_matrix) <- 0

	dist_matrix <- as.dist(dist_matrix, diag = TRUE)

	h <-hclust(dist_matrix)
	png(filename = pngOut)
	print(plot(h))
	dev.off()

	clust <- getOptimalClusterNumber(h,full_dist_matrix,20)

	write.table(clust,
	            file=clustOut,
	            quote = FALSE,
	            row.names = FALSE,
	            col.names = FALSE,
	            sep = ",")';

	close RS;

	return $scriptName;
}
