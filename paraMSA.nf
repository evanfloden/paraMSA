
/* 
 * Main paraMSA pipeline script
 *
 * @authors
 * Jia-Ming Chang <chang.jiaming@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Cedric Notredame <cedric.notredame@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.name          = "Concatenated MSAs and Consensus Phylogentic Trees Analysis"
params.seq           = "$baseDir/tutorial/data/tips16_0.5_001.0400.fa"
params.ref           = "$baseDir/tutorial/data/asymmetric_0.5.unroot.tree"
params.output        = "$baseDir/tutorial/results"
params.straps        = 10
params.alignments    = 15
params.aligner       = "CLUSTALW"


log.info "p a r a M S A ~  version 0.1"
log.info "====================================="
log.info "name                      : ${params.name}"
log.info "input sequence (FASTA)    : ${params.seq}"
log.info "reference tree (NEWICK)   : ${params.ref}"
log.info "output (DIRECTORY)        : ${params.output}"
log.info "straps                    : ${params.straps}"
log.info "alignments                : ${params.alignments}"
log.info "aligner                   : ${params.aligner}"
log.info "\n"


/*
 * Input parameters validation
 */

alignments_num = params.alignments as int
straps_num     = params.straps as int

/*
 * Create a channel for input sequence files 
 */
 
Channel
    .fromPath( params.input )
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .map { file -> tuple( file.baseName, file ) }
    .set { datasets }

Channel
    .fromPath ( params.ref )
     .ifEmpty { error "Cannot find the reference tree file matching: ${params.ref}" }
    .set { referenceTrees }


/*
 * Generate the alternative alignments using guidance2
 */

process alignments {
    tag "guidance2: $datasetID"

    publishDir "${params.output}/${params.aligner}/${datasetID}/alignments", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file(datasetFile) from datasets

    output:
    set val(datasetID), file ("alternativeMSA") into alignments_1, alignments_2
    set val(datasetID), file ("base_alignment.aln") into baseAlignments, baseAlignmentsForBootstrap   
  
    script:
    //
    // Guidance2: Generating alternative alignments
    //
    
    """
    perl \$GUIDANCE2DIR/www/Guidance/guidance.pl \
        --seqFile ${datasetFile} \
        --msaProgram ${params.aligner} \
        --seqType aa \
        --outDir \$PWD \
        --mafft mafft \
        --clustalw clustalw2 \
        --prank prank

    cp MSA.${params.aligner}.aln.With_Names base_alignment.aln

    mkdir \$PWD/alternativeMSA
    tar -zxf MSA.${params.aligner}.Guidance2_AlternativeMSA.tar.gz -C alternativeMSA

    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    files=(\$(ls "alternativeMSA/"|sort -R |tail --lines=${alignments_num}))

    ls "alternativeMSA/" |while read file; 
    do
    if containsElement "\$file" "\${files[@]}"
        then
        echo "Keeping \$file";
        else
            echo "Removing \$file";
            rm "alternativeMSA/\$file"
        fi
    done

    """
}

process create_strap_alignments {
    tag "strap alignments: $datasetID"
    publishDir "${params.output}/${params.aligner}/${datasetID}/strap_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alignments_2
    set val(datasetID), file(baseAlignment) from baseAlignments


    output:
    set val (datasetID), file("paramastrap_phylips") into paramastrapPhylipsDir
    set val (datasetID), file("bootstrap.phylip") into bootstrapPhylips

    script:
    //
    // This task generates strap alignments in phylip format.
    // For each alternative alignment, ${straps_num} strap alignments are generated.
    // The base (default) alignment is also bootstrapped 
    //

    """
    seed=4533
    alternativeMSAsFASTA=( ${alternativeAlignmentsDir}/*".fasta" )
    alternativeMSAsPHYLIP=("\${alternativeMSAsFASTA[@]/.fasta/.phylip}")
    straps=\$((${straps_num} + 1))

    x=0
    mkdir -p paramastrap_phylips/alternativeMSA
    for i in "\${alternativeMSAsFASTA[@]}"
    do
        outfileBase=\${i%.fasta}
        esl-reformat phylip \$i > \$outfileBase.phylip
        echo -e "\$outfileBase.phylip\nR\n\${straps}\nY\n\$seed\n" | seqboot
        mv outfile paramastrap_phylips/\${outfileBase}_strap.phylip
        mv paramastrap_phylips/\${outfileBase}_strap.phylip paramastrap_phylips/.
    done
    rm -rf paramastrap_phylips/alternativeMSA

    esl-reformat phylip ${baseAlignment} > base.phylip
    echo -e "base.phylip\nR\n\${straps}\nY\n\$seed\n" | seqboot
    mv outfile bootstrap.phylip

    """
}


paramastrapPhylipsDir.flatMap {  id, dir -> dir.listFiles().collect { [id, it] } }
                     .view()
                     .set { paramastrapPhylips }

/*
 *  Split the strap alignment phylips into chunks of one alignment:
 *
 *   For every item in paramastrapPhylips (contains a set of val(datasetID), file(multi_phylip_file)) :
 *       get the alignmentID from multi_phylip_file.baseName
 *       split file(multi_phylip_file) with splitPhylip
 *       create a new channel that emits a set containing val(datasetID), val(alignmentID), val(strapID), file(single_phylip)
 */

def splitPhylip(file) {
    def chunks = Channel.create() 
    def buffer = new StringBuffer()
    println file
    file.eachLine { line ->
      if( line ==~ /^ +\d+\s+\d+/ ) {
         if( buffer.length() ) 
             chunks << buffer.toString()
         buffer.length = 0
      }
      buffer.append(line).append('\n')

    }

   return chunks
}

/*
 * Split each directory in the channel paramastrapPhylipsDir into 
 */


paramastrapPhylips
    .flatMap { set ->
        def datasetID = set[0]
        def alignmentID = set[1].baseName
        def file =  set[1]
        def strapID = 0 
        splitPhylip(file).collect{ phylip -> tuple(datasetID, alignmentID, strapID++, phylip) } 
     }
     .set { splitPhylips  }  

bootstrapPhylips
    .flatMap { set ->
        def datasetID = set[0]
        def alignmentID = set[1].baseName
        def file =  set[1]
        def strapID = 0
        splitPhylip(file).collect{ phylip -> tuple(datasetID, alignmentID, strapID++, phylip) }
     }
     .set { splitBasePhylips  }


process create_strap_trees {
    tag "bootstrap samples: $datasetID"
    publishDir "${params.output}/${params.aligner}/$datasetID/strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(alignmentID), val(strapID), val(phylip) from splitPhylips
    
    output:
    set val (datasetID), file("${alignmentID}_${strapID}.nwk") into paramastrapTrees
    set val (datasetID), file("${alignmentID}_${strapID}.phylip") into paramastrapSplitPhylips

    script:
    //
    // Generate Strap Trees: Generate strap trees in Newick format
    //

    """
    echo "${phylip}" | tee ${alignmentID}_${strapID}.phylip  
    FastTree ${alignmentID}_${strapID}.phylip > ${alignmentID}_${strapID}.nwk

    """
}

process create_default_strap_trees {
    tag "bootstrap samples: $datasetID"
    publishDir "${params.output}/${params.aligner}/$datasetID/strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(alignmentID), val(strapID), val(phylip) from splitBasePhylips

    output:
    set val (datasetID), file("${alignmentID}_${strapID}.nwk") into paramastrapBaseTrees
    set val (datasetID), file("${alignmentID}_${strapID}.phylip") into paramastrapBaseSplitPhylips

    script:
    //
    // Generate Strap Trees: Generate strap trees in Newick format
    //

    """
    echo "${phylip}" | tee ${alignmentID}_${strapID}.phylip
    FastTree ${alignmentID}_${strapID}.phylip > ${alignmentID}_${strapID}.nwk

    """
}

process tree_lists {
    tag "tree_list: Alignments-${x} Bootstraps-${y} {dataset_ID}"
    publishDir "$results_path/CLUSTALW/$datasetID/supportSelection", mode: 'copy', overwrite: 'true'

    input:
    each x from [1,2,5,10,25,50,100,400]
    each y from [1,2,5,10,25,50,100]
    set val (datasetID), file (alternativeAlignmentsDir) from alternativeAlignmentsDirectories
    set val (datasetID), file (supportTreesDir) from supportTreesDirectories

    output:
    set val(datasetID), file("${x}_${y}_supportFileList.txt") into supportFileLists
    set val(datasetID), file("${x}_${y}_concatenatedSupportTrees.nwk") into concatenatedSupportTrees

    script:
    """
    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    perl $baseDir/bin/hhsearch_cluster.pl ${alternativeAlignmentsDir} a3m $x ${x}_${y}_alignmentFileList.txt

    while IFS=\\= read var_${x}_${y} value_${x}_${y}; do
        vars_${x}_${y}+=(\$var_${x}_${y})
        values_${x}_${y}+=(\value_${x}_${y})
    done < ${x}_${y}_alignmentFileList.txt

    $alternativeAlignmentsDir | while read file;
    do
    if containsElement \$file \${vars_${x}_${y}[@]}
        then
            echo "Selected: \$file for ${y} bootstraps"
            name=\${file%.fasta}
            ls -1a ${supportTreesDir}/\${name}_strap_{0..${y}}.nwk >> ${x}_${y}_supportFileList.txt
        else
            echo "Skipping: \$file"
        fi
    done

    while read p; do
      cat \$p >> ${x}_${y}_concatenatedSupportTrees.nwk
    done <${x}_${y}_supportFileList.txt

    """
}

process node_support {
    
    tag "node_support: treeID-${tree_ID} ${dataset_ID}"
    publishDir "$results_path/CLUSTALW/$datasetID/nodeSupport", mode: 'copy', overwrite: 'true'

    input: 
    set val(datasetID), file (
    file (referenceTree) from referenceTrees
    set val(datasetID), val(treeID), file (phylogeneticTree) from phlyogeneticTrees

    output:
    set val(datasetID), file("${x}_${y}_filelist.txt") into fileLists
    set val(datasetID), file("${x}_${y}_concatenatedSupportTrees.nwk") into concatenatedSupportTrees

    script:
    """
    t_coffee -other_pg seq_reformat \
             -in $phylogeneticTree \
             -in2 supportList.txt \
             -action +tree2ns $referenceTree \
             > nodeSupportFor_${treeID}_Tree.result
    """
}




}
  if [[ $datasetID =~ ^tips([0-9]+)_([0-9].[0-9])_[0-9]+.[0-9]+\$ ]]
    then
        reference_tree=${baseDir}/tutorial/data/ref_trees/tips\${BASH_REMATCH[1]}/asymmetric_\${BASH_REMATCH[2]}.unroot.tree
    else
        echo "unable to parse string ${datasetID}"
    fi

    base_tree=${resultsDir}/strap_trees/bootstrap_0.nwk;

    for file in ${resultsDir}/strap_trees/bootstrap_*.nwk;
    do
        cat \$file >> base_10_TreesConcat.nwk
    done

    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    perl $baseDir/bin/hhsearch_cluster.pl ${resultsDir}/alignments/alternativeMSA a3m 1 1_10_alignmentList.txt
    perl $baseDir/bin/hhsearch_cluster.pl ${resultsDir}/alignments/alternativeMSA a3m 10 10_1_alignmentList.txt

    while IFS=\\= read var1x10 value1x10; do
        vars1x10+=(\$var1x10)
        values1x10+=(\$value1x10)
    done < 1_10_alignmentList.txt

    while IFS=\\= read var10x1 value10x1; do
        vars10x1+=(\$var10x1)
        values10x1+=(\$value10x1)
    done < 10_1_alignmentList.txt

    ls ${resultsDir}/alignments/alternativeMSA | while read file;
    do
    if containsElement \$file \${vars1x10[@]}
        then
            echo "Selected: \$file for 10 bootstraps"
            name=\${file%.fasta}
            ls -1a ${resultsDir}/strap_trees/\${name}_strap_{0..9}.nwk >> 1x10_filelist.txt
        else
            echo "Skipping: \$file"
        fi
    done

    while read p; do
      cat \$p >> 1x10_TreesConcat.nwk
    done <1x10_filelist.txt

    echo ">1x10_TreesConcat.nwk" > supportList.txt

    ls ${resultsDir}/alignments/alternativeMSA | while read file;
    do
    if containsElement \$file \${vars10x1[@]}
        then
            echo "Selected: \$file for 1 bootstraps"
            name=\${file%.fasta}
            ls -1a ${resultsDir}/strap_trees/\${name}_strap_0.nwk >> 10x1_filelist.txt
        else
            echo "Skipping: \$file"
        fi
    done

    while read p; do
      cat \$p >> 10x1_TreesConcat.nwk
    done <10x1_filelist.txt

    echo ">10x1_TreesConcat.nwk" >> supportList.txt


    t_coffee -other_pg seq_reformat \
              -in \$base_tree \
              -in2 supportList.txt \
              -action +tree2ns \$reference_tree \
              > nodeSupportForBaseTree.result

    """
}

