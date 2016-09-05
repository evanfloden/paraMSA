#!/usr/bin/env perl

use strict;
use warnings;
use File::Find::Rule;

# Take as input a directory containing results with the following filename format:
# CLUSTALW/tips16_2.0_001.0400/nodeSupport/nodeSupportFor_tips16_2.0_001.0400_random_100_Tree.result


my %da;
my %depth;
my $datasetID;
my ($tips, $sym, $len, $set, $tree);

my $dir = $ARGV[0];

my @files = File::Find::Rule->file()
                            ->name("*.result")
                            ->in($dir);


foreach my $resultFile (@files) {
    #print "Results File $resultFile\n";
    if ($resultFile=~/[\s\S]+\/(CLUSTALW|MAFFT|PRANK)\/[\s\S]+\/nodeSupport\/nodeSupportFor_(tips\d+)_(\d.\d)_(\d+).(\d+)_(\S+)_Tree.result/) {
        my $aligner =$1;
        my $tips    =$2;
        my $sym     =$3;
        my $set     =$4;
        my $len     =$5;
        my $tree    =$6;
       	my $datasetID="$1_$2_$3";
	print "DatasetID $datasetID\n";
        my @bootstrapMethodList;
        my $node;
        
        open (F, $resultFile);
        print "Open $resultFile\n";
        while (<F>) {
            my $line=$_;
            # #REPLICATE:   2 partialBootstrapTreesConcat.nwk   15
	    if ( ($line=~/\#REPLICATE:/)) {
                my $bootstrapMethod;
                # weightedBootstrapTreesConcat.nwk
                #print "LINE=$line\n";
		$line=~/([\S_]+)_concatenatedSupportTrees.nwk/;
                $bootstrapMethod=$1;
                @bootstrapMethodList=(@bootstrapMethodList,$bootstrapMethod);
            }
            elsif ( !($line=~/\#/)) {
                $node++;
                my @v=split ('\s+', $line);
                my $count=3;
                foreach my $bootstrapMethod (@bootstrapMethodList) {

                    #             0            1        2          $count < 3,4,5 etc..>
                    #           <list>      <depth> <correct> <support from replicate X..>
                    # @v = 1100000000000000    2        1     1.0000 0.9333 0.5333 1.0000

                    $depth{$datasetID}{$v[1]}++;
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'list'}           =$v[0];
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'depth'}          =$v[1];
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'treeCorrect'}    =$v[2]; #Is the tree correct: 0 or 1
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'pp'}                   +=$v[2]; #Number nodes correct in tree
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'nt'}++;
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'}          =$v[$count]; #BS support for the node
		    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{bs}                +=$v[$count]; #Sum BS supports in tree    
		    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{nnode}++;            
		    #       print "Dataset $datasetID\nTree $tree\nBootstrap $bootstrapMethod\nSet $set\nNode $node\n\n";
                    $count++;
	        
                    # Avg Number of Correct Nodes
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}=$da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'pp'}/$da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'tree'}{'nt'};

                    # Avg Accuracy Score of Node 
                    # If Tree is 0 (incorrect) and BS Method is 0.1, score is 0.9.  
                    # If Tree is 1 (correct) and BS Method is 0.8, score is 0.8.              
                    my $score;
                    if ($da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'treeCorrect'}) {
                        $score = $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'};
                    }
                    else {
                        $score = 1 - $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'bs'};
                    }
                              
                    $da{$datasetID}{$tree}{$bootstrapMethod}{$set}{'node'}{$node}{'score'} = $score;               

                }
            }
        }
        close(F);   
    }
}


my @list_exp               = sort (keys (%da));
my @list_trees             = sort (keys (%{$da{$list_exp[0]}}));
my @list_bootstrapMethods  = sort (keys (%{$da{$list_exp[0]}{$list_trees[0]}}));
my @list_sets              = sort (keys (%{$da{$list_exp[0]}{$list_trees[0]}{$list_bootstrapMethods[0]}}));
my @list_tips              = ('tips16');
my @list_syms              = ('2.0');
my @list_aligners          = ('CLUSTALW', 'MAFFT');

print "@list_exp\n";
print "@list_trees\n";
print "@list_bootstrapMethods\n";
print "@list_sets\n";


my %ma;
my %T;
my %bacc;
my $nexp;
foreach my $exp (sort (keys (%da))) {
    $nexp++;
    print "exp = |$exp|\n";
    foreach my $tree (sort (keys (%{$da{$exp}}))) {
        my $j=1;
        my $y;
        my $bootstrapMethod = $list_bootstrapMethods[0];
	foreach my $set (sort (keys %{$da{$exp}{$tree}{$bootstrapMethod}})) {
            my @ll=keys(%{$da{$exp}{$tree}{$bootstrapMethod}});
            my $nset=$#ll+1;
            $y->[$j][1]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'};
            $y->[$j][2]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'};
            $j++;

            $ma{$exp}{$tree}{$bootstrapMethod}{'acc'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}/$nset;
            $ma{$exp}{$tree}{$bootstrapMethod}{'bs'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'}/$nset;
            
        }
    }
}

print "Table 1.1: Average Topological Accuracy\n";
foreach my $header (@list_exp) {
  print "$header\t";
}
print "\n";
foreach my $tree (@list_trees) {
  printf "$tree\t";
  foreach my $exp (@list_exp) {
    my $bootstrapMethod = $list_bootstrapMethods[0];
    printf "%.3f\t", $ma{$exp}{$tree}{$bootstrapMethod}{acc};
  }    
  print "\n";
}
print "\n";


## ROC Data


my $filename = "ROC.data";
open my $fh, '>', $filename or die;

foreach my $exp (@list_exp) { 
  my @tmpArray = split(/_/,$exp);
  my $tree_for_comparison = $tmpArray[0];
  $tree_for_comparison = lc($tree_for_comparison);
  foreach my $bs_method(@list_bootstrapMethods) {
    print $fh "${exp}_${bs_method}.labels,";
    foreach my $set (sort (keys %{$da{$exp}{$tree_for_comparison}{$bs_method}})) {
      foreach my $node (keys %{$da{$exp}{$tree_for_comparison}{$bs_method}{$set}{node}}) {
        print $fh "$da{$exp}{$tree_for_comparison}{$bs_method}{$set}{'node'}{$node}{'treeCorrect'},";
      }
    }
  }
print $fh "\b \b\n";
}


foreach my $exp (@list_exp) {
  my @tmpArray = split(/_/,$exp);
  my $tree_for_comparison = $tmpArray[0];
  $tree_for_comparison = lc($tree_for_comparison);
  print "T4C =|$tree_for_comparison|\n";
  foreach my $bs_method(@list_bootstrapMethods) {
    print $fh "${exp}_${bs_method}.labels,";
    foreach my $set (sort (keys %{$da{$exp}{$tree_for_comparison}{$bs_method}})) {
      foreach my $node (keys %{$da{$exp}{$tree_for_comparison}{$bs_method}{$set}{node}}) {
        print $fh "$da{$exp}{$tree_for_comparison}{$bs_method}{$set}{'node'}{$node}{'treeCorrect'},";
      }
    }
  } 
print $fh "\b \b\n";
}



sub list2bmcc
  {
     my ($x)=@_;
     my $n= scalar(@{$x});
     my ($bmcc,$pbs,$bbs);
     my ($tp,$tn,$fp,$fn);
     my @array;
     for (my $i=0; $i<$n; $i++)
       {
         $array[$i][0]=0+$x->[$i][0];
         $array[$i][1]=0+$x->[$i][1];
         $tp+=0+$array[$i][0];
         $fp+=1-$array[$i][0];
       }
     @array=sort { $a->[1] <=> $b->[1]} @array;


     $bmcc=-999;
     $pbs=-1;
     $bbs=-1;
     for (my $i=0; $i<$n; $i++)
       {
         my $bs=$array[$i][1];
         my $pp=$array[$i][0];
         my $pn=1-$pp;
         if ($pbs!=-1 && $bs!=$pbs)
           {
             my $matN=($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn);
             my $mcc=($matN==0)?-999:(($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
             if ($mcc >$bmcc)
               {
                 $bmcc=$mcc; $bbs=$pbs;
               }
           }
         $tp-=$pp;
         $tn+=$pn;
         $fp-=$pn;
         $fn+=$pp;
         $pbs=$bs;
       }
     my $matN=($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn);
     my $mcc=($matN==0)?-999:(($tp*$tn)-($fp*$fn))/sqrt(($tp+$fp)*($tp+$fn)*($tn+$fp)*($tn+$fn));
     if ($mcc >$bmcc){$bmcc=$mcc; $bbs=$pbs;}
     return $bmcc;
   }

