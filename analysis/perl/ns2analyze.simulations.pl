#!/usr/bin/env perl

use strict;
use warnings;

# Take as input a directory containing nodesupport reslts with the following filename format:
# tips16_0.5_0400_001_baseTree.result

my %da;
my %depth;
my $datasetID;
my ($tips, $sym, $len, $set, $tree);

my @fileArray;
my $dir=$ARGV[0];
opendir (DIR, $dir) or die $!;
while (my $file = readdir(DIR)) {
  push @fileArray, $file;
}
close(DIR);

foreach my $resultFile (@fileArray) {
    #print "Results File $resultFile\n";
    if ($resultFile=~/(tips\d+)_(\d.\d)_(\d+)_(\d+)_(Base|Concatenated)Tree.result/) {
        $datasetID="$1_$2";
	#print "DatasetID $datasetID\n";
        my $tips = $1;
	my $sym  = $2;
	my $len  = $3;
	my $set  = $4;
	my $tree = $5;

        my @bootstrapMethodList;
        my $node;
        
        open (F, "$dir/$resultFile");
        print "Open $resultFile\n";
        while (<F>) {
            my $line=$_;
            # #REPLICATE:   2 partialBootstrapTreesConcat.nwk   15
	    if ( ($line=~/\#REPLICATE:/)) {
                my $bootstrapMethod;
                # weightedBootstrapTreesConcat.nwk
                #print "LINE=$line\n";
		$line=~/(partial|weighted|base|super)BootstrapTreesConcat.nwk/;
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
                    print "Dataset $datasetID\nTree $tree\nBootstrap $bootstrapMethod\nSet $set\nNode $node\n\n";
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
my @list_trees             = sort (keys (%{$da{@list_exp[0]}}));
my @list_bootstrapMethods  = sort (keys (%{$da{@list_exp[0]}{$list_trees[0]}}));
my @list_sets              = sort (keys (%{$da{@list_exp[0]}{$list_trees[0]}{$list_bootstrapMethods[0]}}));




my %ma;
my %T;
my %bacc;
my $nexp;
foreach my $exp (sort (keys (%da))) {
    $nexp++;
    foreach my $tree (sort (keys (%{$da{$exp}}))) {
        my $j=1;
        my $y;
        my @bootstrapMethods = (sort (keys (%{$da{$exp}{$tree}})));
        my $bootstrapMethod = shift(@bootstrapMethods);
        foreach my $set (sort (keys %{$da{$exp}{$tree}{$bootstrapMethod}})) {
            my @ll=keys(%{$da{$exp}{$tree}{$bootstrapMethod}});
            my $nset=$#ll+1;
            $y->[$j][1]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'};
            $y->[$j][2]=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'};
            $j++;

            $ma{$exp}{$tree}{$bootstrapMethod}{'acc'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}/$nset;
            $ma{$exp}{$tree}{$bootstrapMethod}{'bs'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'}/$nset;
            
            $ma{all}{$tree}{$bootstrapMethod}{'acc'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'acc'}/$nset;
            $ma{all}{$tree}{$bootstrapMethod}{'bs'}+=$da{$exp}{$tree}{$bootstrapMethod}{$set}{'tree'}{'bs'}/$nset;
            
        }
    }
}


my @treeList=('Base','Concatenated');
my @methodsList=('base','partial','super','weighted');
my @tipsList=('tips16', 'tips32', 'tips64');
my @symList=('0.5','1.0','2.0');

print "Table 1.1: Average Topological Accuracy\n";
foreach my $tree (@treeList) {
    print "$tree ";
    foreach my $tips (@tipsList) {
        foreach my $sym (@symList) {
            my $exp="$tips\_$sym";
            my @bootstrapMethods = (keys (%{$da{$exp}{$tree}}));
            my $bootstrapMethod = shift(@bootstrapMethods);
            printf "%.3f ", $ma{$exp}{$tree}{$bootstrapMethod}{acc};
         }
    }
    printf "%.3f ",$ma{all}{$tree}{'base'}{acc}/$nexp;
    print "\n";
}



## ROC Data


my $filename = "ROC.data";
open my $fh, '>', $filename or die;

foreach my $method (@methodsList) {
  print $fh "$method.labels,";
  foreach my $tips (@tipsList) {
    foreach my $sym (@symList) {
      my $exp="$tips\_$sym";
      foreach my $tree (@treeList) {
          foreach my $set (sort (keys %{$da{$exp}{$tree}{$method}})) {
              foreach my $node (keys %{$da{$exp}{$tree}{$method}{$set}{node}}) {
                  print $fh "$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'treeCorrect'},";
          }
        }
      }
    }
  }
  print $fh "\b \b\n";
}


foreach my $method (@methodsList){
  print $fh "$method.predictions,";
  foreach my $tips (@tipsList) {
    foreach my $sym (@symList) {
      my $exp="$tips\_$sym";
      foreach my $tree (@treeList) {
          foreach my $set (sort (keys %{$da{$exp}{$tree}{$method}})) {
              foreach my $node (keys %{$da{$exp}{$tree}{$method}{$set}{node}}) {
                  print $fh "$da{$exp}{$tree}{$method}{$set}{'node'}{$node}{'bs'},";
          }
        }
      }
    }
  }
  print $fh "\b \b\n";
}


# Best Mathews Correlation Coefficent

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


## Accuracy 
