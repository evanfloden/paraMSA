
/* 
 * Main paraMSA pipeline script
 *
 * @authors
 * Jia-Ming Chang <chang.jiaming@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Cedric Notredame <cedric.notredame@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */

params.name             = "tutorial_data"
params.seq              = "$baseDir/tutorial/data/tips16_0.5_00{1,2}.0400.fa"
params.ref              = "$baseDir/tutorial/data/asymmetric_0.5.unroot.tree"
params.output           = "$baseDir/tutorial/results/"
params.straps           = 16
params.straps_list      = '1,2,4,8,16'
params.alignments       = 16
params.alignments_list  = '1,2,4,8,16'
params.aligner          = "CLUSTALW"

log.info "p a r a M S A ~  version 0.1"
log.info "====================================="
log.info "name                            : ${params.name}"
log.info "input sequence (FASTA)          : ${params.seq}"
log.info "reference tree (NEWICK)         : ${params.ref}"
log.info "output (DIRECTORY)              : ${params.output}"
log.info "straps                          : ${params.straps}"
log.info "straps list                     : ${params.straps_list}"
log.info "alignments [multiple of 4]      : ${params.alignments}"
log.info "alignments list                 : ${params.alignments_list}"
log.info "aligner [CLUSTALW|MAFFT|PRANK]  : ${params.aligner}"
log.info "\n"

/*
 * Input parameters validation
 */

alignments_num     = params.alignments as int
straps_num         = params.straps as int

strapsList         = params.straps_list.tokenize(',')
alignmentsList     = params.alignments_list.tokenize(',')

def guidanceBootraps = alignments_num / 4 

/*
 * Create a channel for input sequence files 
 */



Channel
    .fromPath( params.seq )
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.seq}" }
    .map { file -> tuple( file.baseName, file ) }
    .into { datasetsA; datasetsB }



/*
 * Generate the alternative alignments using guidance2
 * Outputs a directory containing the alternative alignments and the default alignment
 */

process alignments {
    tag "guidance2: $datasetID"

    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/alignments", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file(datasetFile) from datasetsA

    output:
    set val(datasetID), file ("alternativeMSA") into alternativeAlignmentDirectories_A, alternativeAlignmentDirectories_B,  alternativeAlignmentDirectories_C
  
    script:
    //
    // Guidance2: Generating alternative alignments
    //
    
    """
    perl ${params.guidance2dir}/www/Guidance/guidance.pl \
        --seqFile ${datasetFile} \
        --msaProgram ${params.aligner} \
        --seqType aa \
        --outDir \$PWD \
        --bootstraps ${guidanceBootraps} \
        --mafft mafft \
        --clustalw clustalw2 \
        --prank prank

    mkdir \$PWD/alternativeMSA
    tar -zxf MSA.${params.aligner}.Guidance2_AlternativeMSA.tar.gz -C alternativeMSA

    """
}

/*
 *  Create the default alignments using CLUSTALW2, MAFFT, PRANK and T-COFFEE
 */

process default_alignments {
    tag "default alignments: $datasetID"

    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/default_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(datasetFile) from datasetsB

    output:
    set val(datasetID), val("clustalw"), file ("${datasetID}_default_clustal_alignment.fa") into defaultClustalAlignments
    set val(datasetID), val("mafft"), file ("${datasetID}_default_mafft_alignment.fa") into defaultMafftAlignments
    set val(datasetID), val("prank"), file ("${datasetID}_default_prank_alignment.fa") into defaultPrankAlignments
    set val(datasetID), val("tcoffee"), file ("${datasetID}_default_tcoffee_alignment.fa") into defaultTcoffeeAlignments
   

    script:
    //
    // Default Alignments
    //

    """
    clustalw2 ${datasetFile} -outfile=${datasetID}_default_clustal_alignment.aln
    esl-reformat afa ${datasetID}_default_clustal_alignment.aln > ${datasetID}_default_clustal_alignment.fa

    mafft ${datasetFile} > ${datasetID}_default_mafft_alignment.aln
    esl-reformat afa ${datasetID}_default_mafft_alignment.aln > ${datasetID}_default_mafft_alignment.fa

    prank -d=${datasetFile}
    mv output.best.fas ${datasetID}_default_prank_alignment.aln 
    esl-reformat afa ${datasetID}_default_prank_alignment.aln > ${datasetID}_default_prank_alignment.fa

    t_coffee -in ${datasetFile} -outfile ${datasetID}_default_tcoffee_alignment.aln
    esl-reformat afa ${datasetID}_default_prank_alignment.aln >${datasetID}_default_tcoffee_alignment.fa

    """
}


defaultClustalAlignments
    .concat(defaultMafftAlignments, defaultPrankAlignments, defaultTcoffeeAlignments)
    .set { defaultAlignments }
   
defaultAlignments.into { defaultAlignmentsA; defaultAlignmentsB }

process default_phylogenetic_trees {
    tag "default phylogenetic trees: $datasetID - $aligner"

    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/default_phylogenetic_trees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(aligner), file (alignment) from defaultAlignmentsA

    output:
    set val(datasetID), val(aligner), file("${datasetID}_${aligner}_phylogeneticTree.nwk") into defaultPhylogeneticTrees
 
    script:
    //
    // Phylogenetic Trees from Default Alignment
    //

    """
    FastTree ${alignment} > ${datasetID}_${aligner}_phylogeneticTree.nwk

    """
}

/*
 * Uses seqboot to create the bootstrap alignments in phylip format
 * Outputs a directory containing the bootstrapped alternative alignment replicates 
 * as well as the bootstrapped default alignment replicates
 */

process create_strap_alignments {
    tag "strap alignments: $datasetID"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/strap_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alternativeAlignmentDirectories_A

    output:
    set val (datasetID), file("paramastrap_phylips") into paramastrapPhylipsDir

    script:
    //
    // This task generates strap alignments in phylip format.
    // For each alternative alignment, ${straps_num} strap alignments are generated.
    //

    """
    seed=4533
    alternativeMSAsFASTA=( ${alternativeAlignmentsDir}/*".fasta" )
    alternativeMSAsPHYLIP=("\${alternativeMSAsFASTA[@]/.fasta/.phylip}")
    straps=\$((${straps_num} + 1))

    x=0
    for i in "\${alternativeMSAsFASTA[@]}"
    do
        outfileBase=\${i%.fasta}
        outfileName=\${outfileBase##*/}

        esl-reformat phylip \$i > \${outfileName}.phylip

        echo -e "\$outfileName.phylip\nR\n\${straps}\nY\n\$seed\n" | seqboot

        mv outfile \${outfileName}_strap.phylip

    done

    mkdir paramastrap_phylips

    mv *_strap.phylip paramastrap_phylips/.
    """
}


paramastrapPhylipsDir.flatMap {  id, dir -> dir.listFiles().collect { [id, it] } }
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

/*
 * Create bootstrap trees from each boottrap alignment replicate
 * Output both the NEWICK tree and the bootstrap alignment used to generate it
 */


process strap_trees {
    tag "bootstrap samples: ${datasetID} - ${alignmentID} - ${strapID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(alignmentID), val(strapID), val(phylip) from splitPhylips
    
    output:
    set val (datasetID), file("${alignmentID}_${strapID}.nwk") into paramastrapTrees
    set val (datasetID), file("${datasetID}_${alignmentID}_${strapID}.phylip") into paramastrapSplitPhylips

    script:
    //
    // Generate Strap Trees: Generate strap trees in Newick format
    //

    """
    echo "${phylip}" | tee ${datasetID}_${alignmentID}_${strapID}.phylip  
    FastTree ${datasetID}_${alignmentID}_${strapID}.phylip > ${alignmentID}_${strapID}.nwk

    """
}


process create_default_alignment_straps {

    tag "default phylogenetic trees: $datasetID - $aligner"

    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/default_phylogenetic_trees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(aligner), file (alignment) from defaultAlignmentsB

    output:
    set val(datasetID), val(aligner), file("${datasetID}_${aligner}_default_bootstraps.phylip") into defaultBootstrapPhylips

    script:
    //
    // This task generates strap alignments in phylip format.
    // For each deafault alignment, ${straps_num} strap alignments are generated.
    //

    """
    seed=4533
    straps=\$((${straps_num} + 1))
    esl-reformat phylip ${alignment} > ${aligner}_default.phylip
    echo -e "${aligner}_default.phylip\nR\n\${straps}\nY\n\$seed\n" | seqboot
    mv outfile ${datasetID}_${aligner}_default_bootstraps.phylip

    """
}

defaultBootstrapPhylips
    .flatMap { set ->
        def datasetID = set[0]
        def aligner = set[1]
        def file =  set[2]
        def strapID = 0
        splitPhylip(file).collect{ phylip -> tuple(datasetID, aligner, strapID++, phylip) }
     }
     .set { defaultBootstrapPhylipsSplit }

process create_default_strap_trees {
    tag "default_strap_trees: $datasetID"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(aligner), val(strapID), val(phylip) from defaultBootstrapPhylipsSplit

    output:
    set val (datasetID), val(aligner), val(strapID), file("${datasetID}_${aligner}_${strapID}.nwk") into defaultAlignmentsBootstrapTrees
    set val (datasetID), val(aligner), val(strapID), file("${datasetID}_${aligner}_${strapID}.phylip") into defaultAlignmentsBootstrapPhylips

    script:
    //
    // Generate Strap Trees: Generate strap trees in Newick format for the default alignment bootstraps
    //

    """
    echo "${phylip}" | tee ${datasetID}_${aligner}_${strapID}.phylip
    FastTree ${datasetID}_${aligner}_${strapID}.phylip > ${datasetID}_${aligner}_${strapID}.nwk

    """
}

process alternative_alignment_concatenate {
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/concatenated_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alternativeAlignmentDirectories_B
    each treeID from alignmentsList

    output:
    set val(datasetID), val(treeID), file("distant_${treeID}_concatenated_alignments.phylip") into distantConcatenatedAlignmentPhylips
    set val(datasetID), val(treeID), file("distant_${treeID}_concatenated_alignments.aln") into distantConcatenatedAlignmentFastas

    set val(datasetID), val(treeID), file("random_${treeID}_concatenated_alignments.phylip") into randomConcatenatedAlignmentPhylips
    set val(datasetID), val(treeID), file("random_${treeID}_concatenated_alignments.aln") into randomConcatenatedAlignmentFastas

    //
    // Concatenate Process: Generate the concatenated alignments in PHYLIP format
    // 

    """
    shopt -s nullglob
    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }

    ####
    # Cluster the alignments for the most distant 
    ####
    perl $baseDir/bin/hhsearch_cluster.pl ${alternativeAlignmentsDir} a3m ${treeID} ${treeID}_alignmentFileList.txt

    while IFS=\\= read var_${treeID}; do
        vars_${treeID}+=(\$var_${treeID})
    done < ${treeID}_alignmentFileList.txt

    set selectedAlignmentsArray
    declare -a selectedAlignmentsArray=('')

    while read file; do
        if containsElement \$file "\${vars_${treeID}[@]}"
            then
                echo "Adding \$file to selectedAlignmentsArray";
                selectedAlignmentsArray=(\${selectedAlignmentsArray[@]} ${alternativeAlignmentsDir}/\${file})
            else
                echo "\$file not in vars_${treeID}";
        fi
    done <<<"\$(ls -1 ${alternativeAlignmentsDir})"

    echo "selected = |\${selectedAlignmentsArray[@]}|"

    concatenate.pl --aln \${selectedAlignmentsArray[@]} --out distant_${treeID}_concatenated_alignments.aln
    esl-reformat phylip distant_${treeID}_concatenated_alignments.aln > distant_${treeID}_concatenated_alignments.phylip
    ####


    ####
    # Select alignments randomly
    ####
    set completeAlignmentsArray
    declare -a completeAlignmentsArray=('')
    while read file; 
    do
      echo "Adding \$file to completeAlignmentsArray";
      completeAlignmentsArray=(\${completeAlignmentsArray[@]} ${alternativeAlignmentsDir}/\${file})
    done <<<"\$(ls -1 ${alternativeAlignmentsDir})"

    set randomAlignmentsArray
    declare -a randomAlignmentsArray=('')
    set delete
    declare -a delete=('')

    for ((i=1; i<=$treeID; i++));
    do
      selectedFile=\${completeAlignmentsArray[\$RANDOM % \${#completeAlignmentsArray[@]}]}
      randomAlignmentsArray=(\${randomAlignmentsArray[@]} \${selectedFile})
      delete=(\$selectedFile)
      completeAlignmentsArray=completeAlignmentsArray[@]\\\$delete 
    done

    concatenate.pl --aln \${randomAlignmentsArray[@]} --out random_${treeID}_concatenated_alignments.aln
    esl-reformat phylip random_${treeID}_concatenated_alignments.aln > random_${treeID}_concatenated_alignments.phylip

    """
}

process concatenated_phylogenetic_trees {
 
    tag "concatenated phylogenetic_trees: ${datasetID} - ${treeID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/phylogeneticTrees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(treeID), file(distantConcatenatedAlns) from distantConcatenatedAlignmentFastas
    set val(datasetID), val(treeID), file(randomConcatenatedAlns) from randomConcatenatedAlignmentFastas

    output:
    set val(datasetID), val("distant_${treeID}"), file("${datasetID}_distant-${treeID}_phylogeneticTree.nwk") into distantPhylogeneticTrees
    set val(datasetID), val("random_${treeID}"), file("${datasetID}_random-${treeID}_phylogeneticTree.nwk") into randomPhylogeneticTrees

    script:
    """
    FastTree ${distantConcatenatedAlns} > ${datasetID}_distant-${treeID}_phylogeneticTree.nwk
    FastTree ${randomConcatenatedAlns} > ${datasetID}_random-${treeID}_phylogeneticTree.nwk
    """
}


/*
 * Combine phlyogenetic trees to be grouped by ${datasetID}
 */

distantPhylogeneticTrees
    .concat (randomPhylogeneticTrees, defaultPhylogeneticTrees)
    .set { phylogeneticTrees }

/*
 *  Collect all files from paramastrapTrees from the same dataset 
 *  groupTuple -> emits [ val ($datasetID) [ file("${alignmentID}_${strapID}.nwk"), file("${alignmentID}_${strapID}.nwk") ... ] ] 
 */

paramastrapTrees
    .groupTuple()	
    .set { supportTreesFiles }



process support_tree_lists {
    tag "tree_list: Alignments-${x} Bootstraps-${y} ${datasetID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/supportSelection", mode: 'copy', overwrite: 'true'

    input:
    each x from alignmentsList
    each y from strapsList
    set val (datasetID), file (alternativeAlignmentsDir) from alternativeAlignmentDirectories_C
    set val (datasetID), file (supportTrees) from supportTreesFiles

    output:
    set val(datasetID), file("${x}_${y}_supportFileList.txt") into supportFileLists
    set val(datasetID), file("${x}_${y}_concatenatedSupportTrees.nwk") into concatenatedSupportTrees
    set val (datasetID), val("${x}_${y}_concatenatedSupportTrees.nwk") into supportCombinationsTested

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
    done < ${x}_${y}_alignmentFileList.txt

    declare -i w
    w=\$((${y} - 1))
    while read file; do
        if containsElement \$file "\${vars_${x}_${y}[@]}"
            then
                name=\${file%.fasta}
                if (( w > 0 ))
                  then 
                    eval ls -1a \${name}_strap_{0..\${w}}.nwk >> ${x}_${y}_supportFileList.txt
                  else
                    ls -1a \${name}_strap_0.nwk >> ${x}_${y}_supportFileList.txt
                fi
                echo "Adding \$file to supportFileList.txt";
            else
                echo "\$file not in vars_${x}_${y}";
        fi
    done <<<"\$(ls -1 ${alternativeAlignmentsDir})"

    while read p; do
      cat \$p >> ${x}_${y}_concatenatedSupportTrees.nwk
    done <${x}_${y}_supportFileList.txt

    """
}

concatenatedSupportTrees
    .groupTuple()
    .set{concatenatedSupportTreesGrouped}


/*
 *  a) Create file with list of support combinations to be tested
 *  b) map with a datasetID and the file
 *  c) phase to the concatenatedSupportTrees based on datasetID -> Now emits a list containing [ [datasetID, file(support_list_file)] , [ datasetID, [list of support tree files] ] ]  
 *  d) map this to datasetID -> emits `set val(datasetID), file(support_list_file), list(support tree files)`
 *  e) cross this with phylogeneticTrees -> emits [ 
 *                                                  [ val(datasetID), file(support_list_file), list(support tree files) ],
 *                                                  [ val(datasetID), val(treeID), file("${datasetID}_${treeID}_phylogeneticTree.nwk") ] 
 *                                                ] 
 *  f) map back to datasetID to emit [ val(datasetID), val(treeID), file("${datasetID}_${treeID}_phylogeneticTree.nwk"), file(support_list_file), file(list of support tree files)] 
 */


supportCombinationsTested
    .collectFile() { item -> 
        ["${item[0]}", ">" + item[1] + '\n' ] 
    }
    .map { item -> [ item.name, item ] } 
    .phase(concatenatedSupportTreesGrouped)
    .map { item -> [ item[0][0], item[0][1], item[1][1] ] }
    .cross (phylogeneticTrees)
    .map { item -> [ item[0][0], item[1][1], item[1][2], item[0][1], item[0][2] ] }
    .set { supportCombinationsTestedGrouped }


process node_support {
    
    tag "node_support: ${datasetID} - ${treeID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/nodeSupport", mode: 'copy', overwrite: 'true'

    input: 
    set val(datasetID), val(treeID), file(phylogeneticTree), file(supportFileList), file(listOfSupportTrees) from supportCombinationsTestedGrouped

    output:
    set val(datasetID), val(treeID), file("nodeSupportFor_${datasetID}_${treeID}_Tree.result") into nodeSupports

    script:
    """
    t_coffee -other_pg seq_reformat \
             -in $phylogeneticTree \
             -in2 ${supportFileList} \
             -action +tree2ns ${params.ref} \
             > nodeSupportFor_${datasetID}_${treeID}_Tree.result
    """
}

