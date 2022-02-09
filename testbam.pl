use Bio::DB::Sam;
# use Bio::DB;
use POSIX;
my $sam = Bio::DB::Sam->new(-expand_flags  => 1,
                            -bam  =>"P6_1_12178_ND608_N616_2.sort.markdup.bam");


print strftime("First time read bamFile begin: %H:%M:%S\n", localtime);
my @alignments = $sam->features();
for my $a (@alignments) {
  #  my $seqid  = $a->seq_id;
  # my $seg = $a->isize; #insert_size(can be less than 0)
  #  my $start  = $a->start;
  #  my $end    = $a->end;
  # #  my ($first_mate,$second_mate) = $a->get_SeqFeatures;
  # #  my $fs = $first_mate;
  # #  my $strand = $a->strand;
  # #  my $ref_dna= $a->dna;
 
   my $data  = $a->data;
  # #  my $query_end    = $a->query->end;
  # #  my $query_strand = $a->query->strand;
  # #  my $query_dna    = $a->query->dna;
  #  my $flag =$a->flag;
  #  my $cigar     = $a->cigar_str;
  #  my $scores    = $a->_qscore;     # per-base quality scores
  #  my $match_qual= $a->qual;       # mapping quality
  #  my $paired = $a->get_tag_values('PAIRED'); #get flag_tag 
  #  my $sa = $a->get_tag_values('SA'); #get flag_tag 
  #  my @tag = $a->get_all_tags; #get all tag 
  #  my $length  = $a->length; #read length
  #  my $read = $a->display_name; # read_id
  #  my $type = $a->type; # match or other
  #  my $strand = $a->strand; # 1 or -1
  # my $contig = $a->seq_id; #which contig
  print "$data\n";
   }
print strftime("First time read bamFile begin: %H:%M:%S\n", localtime);
