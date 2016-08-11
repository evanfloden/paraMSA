
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
    .fromPath( params.seq )
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.seq}" }
    .map { file -> tuple( file.baseName, file ) }
    .set { datasets }

/*
 * Generate the alternative alignments using guidance2
 * Outputs a directory containing the alternative alignments and the default alignment
 */

process alignments {
    tag "guidance2: $datasetID"

    publishDir "${params.output}/${params.aligner}/${datasetID}/alignments", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file(datasetFile) from datasets

    output:
    set val(datasetID), file ("alternativeMSA") into alternativeAlignmentDirectories_A, alternativeAlignmentDirectories_B,  alternativeAlignmentDirectories_C
    set val(datasetID), file ("default_alignment.aln") into defaultAlignments
  
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

    cp MSA.${params.aligner}.aln.With_Names default_alignment.aln

    mkdir \$PWD/alternativeMSA
    tar -zxf MSA.${params.aligner}.Guidance2_AlternativeMSA.tar.gz -C alternativeMSA

    """
}


/*
 * Uses seqboot to create the bootstrap alignments in phylip format
 * Outputs a directory containing the bootstrapped alternative alignment replicates 
 * as well as the bootstrapped default alignment replicates
 */


process create_strap_alignments {
    tag "strap alignments: $datasetID"
    publishDir "${params.output}/${params.aligner}/${datasetID}/strap_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alternativeAlignmentDirectories_A
    set val(datasetID), file(defaultAlignment) from defaultAlignments


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

    esl-reformat phylip ${defaultAlignment} > base.phylip
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


/*
 * Create bootstrap trees from each boottrap alignment replicate
 * Output both the NEWICK tree and the bootstrap alignment used to generate it
 */


process strap_trees {
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

process default_strap_trees {
    tag "default_strap_trees: $datasetID"
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

process alternative_alignment_concatenate {
    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alternativeAlignmentDirectories_B
    each z from (1, 2, 5, 10, 25, 50, 100, 400)

    output:
    set val(datasetID), val(z), file("${z}_concatenated_alignments.phylip") into concatenatedAlignmentPhylips

    script:
    //
    // Concatenate Process: Generate the concatenated alignments in PHYLIP format
    // 

    """
    shopt -s nullglob
    alternativeMSAsFASTA=( ${alternativeAlignmentsDir}/*".fasta" )
    alternativeMSAsPHYLIP=("\${alternativeMSAsFASTA[@]/.fasta/.phylip}")

    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    perl $baseDir/bin/hhsearch_cluster.pl ${alternativeAlignmentsDir} a3m ${z} ${z}_alignmentFileList.txt

    while IFS=\\= read var_${z} value_${z}; do
    vars_${z}+=(\$var_${z})
    values_${z}+=(\$value_${z})
    done < ${z}_alignmentFileList.txt

    declare -a selectedAlignmentsArray

    $alternativeAlignmentsDir | while read file;
    do
    if containsElement \$file \${vars_${z}[@]}
        then
            echo "Concatenating: \$file"
            name=\${file%.fasta}
            selectedAlignmentsArray+=(\$file)
        else
            echo "Skipping: \$file"
        fi
    done

    concatenate.pl --aln \${selectedAlignmentsArray[@]} --out ${z}_concatenated_alignments.aln
    esl-reformat phylip ${z}_concatenated_alignments.aln > ${z}_concatenated_alignments.phylip


    """
}


process phylogenetic_trees {

    publishDir "$results_path/CLUSTALW/$datasetID/phylogeneticTrees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(z), file(concatenatedAlns) from concatenatedAlignmentPhylips

    output:
    set val(datasetID), val(z), file("phylogeneticTree.nwk") into phylogeneticTrees

    script:
    """
    
    FastTree ${concatenatedAlns} > phylogeneticTree.nwk

    """
}



/*
 * Collect all files from paramastrapTrees from the same dataset and put them in a directory
 *
 *  paramastrapTrees <- set val (datasetID), file("${alignmentID}_${strapID}.nwk")
 *
 */

paramastrapTrees
    .collectFile() { datasetID, file ->
        def dir = file("$baseDir/work/tmp/datasetID") 
        dir.mkdirs()
        [ dir / file.name, file ]
    }
    .set {supportTreesDirectories}


process support_tree_lists {
    tag "tree_list: Alignments-${x} Bootstraps-${y} {dataset_ID}"
    publishDir "$results_path/CLUSTALW/$datasetID/supportSelection", mode: 'copy', overwrite: 'true'

    input:
    each x from (1,2,5,10,25,50,100,400)
    each y from (1,2,5,10,25,50,100)
    set val (datasetID), file (alternativeAlignmentsDir) from alternativeAlignmentDirectories_C
    set val (datasetID), file (supportTreesDir) from supportTreesDirectories.first()

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
        values_${x}_${y}+=(\$value_${x}_${y})
    done < ${x}_${y}_alignmentFileList.txt

    $alternativeAlignmentsDir | while read file;
    do
    if containsElement \$file \${vars_${x}_${y}[@]}
        then
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
    set val(datasetID), file(concatenatedSupportTrees) from concatenatedSupportTrees
    set val(datasetID), val(treeID), file(phylogeneticTree) from phylogeneticTrees

    output:
    set val(datasetID), val(treeID), file("nodeSupportFor_${treeID}_Tree.result") into nodeSupports

    script:
    """
    t_coffee -other_pg seq_reformat \
             -in $phylogeneticTree \
             -in2 supportList.txt \
             -action +tree2ns ${params.ref} \
             > nodeSupportFor_${treeID}_Tree.result
    """
}

