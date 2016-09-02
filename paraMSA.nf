
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
params.straps_list      = '1,4,16'
params.alignments       = 16
params.alignments_list  = '1,4,16'
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

Channel
    .fromPath ( params.ref )
    .set { refTree }



/**************************
 *
 *   C R E A T E   A L T E R N A T I V E   A L I G N M E N T S
 *
 *   Generate the alternative alignments using guidance2
 *   Outputs a directory containing the alternative alignments and the default alignment
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
 *
 **************************/


/**************************
 *
 *   C R E A T E   D E F A U L T   A L I G N M E N T S
 *
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

    mafft ${datasetFile} > ${datasetID}_default_mafft_alignment.aln

    prank -d=${datasetFile}
    mv output.best.fas ${datasetID}_default_prank_alignment.aln

    t_coffee -in ${datasetFile} -outfile ${datasetID}_default_tcoffee_alignment.aln.tmp
    t_coffee -other_pg seq_reformat -in ${datasetID}_default_tcoffee_alignment.aln.tmp -out ${datasetID}_default_tcoffee_alignment.aln
    
    esl-reformat afa ${datasetID}_default_clustal_alignment.aln > ${datasetID}_default_clustal_alignment.fa
    esl-reformat afa ${datasetID}_default_mafft_alignment.aln > ${datasetID}_default_mafft_alignment.fa
    esl-reformat afa ${datasetID}_default_prank_alignment.aln > ${datasetID}_default_prank_alignment.fa
    esl-reformat afa ${datasetID}_default_tcoffee_alignment.aln > ${datasetID}_default_tcoffee_alignment.fa

    """
}

/*
 *
 **************************/



/**************************
 *
 *   D E F A U L T   A L I G N M E N T   P H L Y O G E N E T I C   T R E E S
 *
 */

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
 *
 **************************/

/**************************
 *
 *         A L T E R N A T I V E   A L I G N M E N T   S E L E C T I O N
 *
 *                                   A N D
 *
 *    C O N C A T E N A T I O N   F O R   P H Y L O G E N E T I C   T R E E S
 */

process alternative_alignment_selection_and_concatenatation {
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/concatenated_alignments", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), file(alternativeAlignmentsDir) from alternativeAlignmentDirectories_B
    each numberOfAlignments from alignmentsList

    output:
    set val(datasetID), val(numberOfAlignments), file("distant_${numberOfAlignments}_concatenated_alignments.phylip") into distantConcatenatedAlignmentPhylips
    set val(datasetID), val(numberOfAlignments), file("distant_${numberOfAlignments}_concatenated_alignments.aln") into distantConcatenatedAlignmentFastas

    set val(datasetID), val(numberOfAlignments), file("random_${numberOfAlignments}_concatenated_alignments.phylip") into randomConcatenatedAlignmentPhylips
    set val(datasetID), val(numberOfAlignments), file("random_${numberOfAlignments}_concatenated_alignments.aln") into randomConcatenatedAlignmentFastas

    set val(datasetID), val(numberOfAlignments), val('random'), file("random_${numberOfAlignments}_alignmentFileList.txt") into randomAlignmentLists
    set val(datasetID), val(numberOfAlignments), val('distant'), file("distant_${numberOfAlignments}_alignmentFileList.txt") into distantAlignmentLists

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
    perl $baseDir/bin/hhsearch_cluster.pl ${alternativeAlignmentsDir} a3m ${numberOfAlignments} distant_alignmentFileList.txt

    while read line; do
      selectedFileMinusPath=\${line##*/}
      selectedFileMinusExtension=\${selectedFileMinusPath%.fasta}
      echo "\$selectedFileMinusExtension" >> distant_${numberOfAlignments}_alignmentFileList.txt
    done < distant_alignmentFileList.txt

    while IFS=\\= read var_${numberOfAlignments}; do
        vars_${numberOfAlignments}+=(\$var_${numberOfAlignments})
    done < distant_alignmentFileList.txt

    set selectedAlignmentsArray
    declare -a selectedAlignmentsArray=('')

    while read file; do
        if containsElement \$file "\${vars_${numberOfAlignments}[@]}"
            then
                selectedAlignmentsArray=(\${selectedAlignmentsArray[@]} ${alternativeAlignmentsDir}/\${file})
        fi
    done <<<"\$(ls -1 ${alternativeAlignmentsDir})"

    concatenate.pl --aln \${selectedAlignmentsArray[@]} --out distant_${numberOfAlignments}_concatenated_alignments.aln
    esl-reformat phylip distant_${numberOfAlignments}_concatenated_alignments.aln > distant_${numberOfAlignments}_concatenated_alignments.phylip
    #######

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

    for ((i=1; i<=$numberOfAlignments; i++));
    do
      selectedFile=\${completeAlignmentsArray[\$RANDOM % \${#completeAlignmentsArray[@]}]}
      randomAlignmentsArray=(\${randomAlignmentsArray[@]} \${selectedFile})
      delete=(\$selectedFile)
      selectedFileMinusPath=\${selectedFile##*/}
      selectedFileMinusExtension=\${selectedFileMinusPath%.fasta}
      echo "\$selectedFileMinusExtension" >> random_${numberOfAlignments}_alignmentFileList.txt
      completeAlignmentsArray=(\${completeAlignmentsArray[@]/\$delete})
    done

    concatenate.pl --aln \${randomAlignmentsArray[@]} --out random_${numberOfAlignments}_concatenated_alignments.aln
    esl-reformat phylip random_${numberOfAlignments}_concatenated_alignments.aln > random_${numberOfAlignments}_concatenated_alignments.phylip

    """
}

/*
 *
 **************************/


/**************************
 *
 *         A L T E R N A T I V E   A L I G N M E N T   P H Y L O G E N E T I C   T R E E S
 */

process concatenated_phylogenetic_trees {

    tag "concatenated phylogenetic_trees: ${datasetID} - ${numberOfAlignments}"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/phylogeneticTrees", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(numberOfAlignments), file(distantConcatenatedAlns) from distantConcatenatedAlignmentFastas
    set val(datasetID), val(numberOfAlignments), file(randomConcatenatedAlns) from randomConcatenatedAlignmentFastas

    output:
    set val(datasetID), val("distant_${numberOfAlignments}"), file("distant_${numberOfAlignments}_phylogeneticTree.nwk") into distantPhylogeneticTrees
    set val(datasetID), val("random_${numberOfAlignments}"), file("random_${numberOfAlignments}_phylogeneticTree.nwk") into randomPhylogeneticTrees

    script:
    """
    FastTree ${distantConcatenatedAlns} > distant_${numberOfAlignments}_phylogeneticTree.nwk
    FastTree ${randomConcatenatedAlns} > random_${numberOfAlignments}_phylogeneticTree.nwk
    """
}

/*
 *
 **************************/



/**************************
 *
 *   D E F A U L T   A L I G N M E N T   B O O T S T R A P   R E P L I C A T E S
 *
 *   Uses seqboot to create the bootstrap alignments in phylip format
 */

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

/*
 *
 **************************/


/**************************
 *
 *   D E F A U L T   A L I G N M E N T   B O O T S T R A P   T R E E S
 *
 */


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
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/default_strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(aligner), val(strapID), val(phylip) from defaultBootstrapPhylipsSplit

    output:
    set val (datasetID), val(aligner), val(strapID), file("${datasetID}_${aligner}_${strapID}.nwk") into defaultAlignmentsBootstrapTrees
    set val (datasetID), val(aligner), val(strapID), file("${datasetID}_${aligner}_${strapID}.phylip") into defaultAlignmentsBootstrapPhylips


    """
    echo "${phylip}" | tee ${datasetID}_${aligner}_${strapID}.phylip
    FastTree ${datasetID}_${aligner}_${strapID}.phylip > ${datasetID}_${aligner}_${strapID}.nwk

    """
}

/*
 *
 **************************/




/**************************
 *
 *   A L T E R N A T I V E   A L I G N M E N T   B O O T S T R A P   R E P L I C A T E S
 *
 *   Uses seqboot to create the bootstrap alignments in phylip format
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
        rm \${outfileName}.phylip
        mv outfile \${outfileName}.phylip
    done

    mkdir paramastrap_phylips
    mv *.phylip paramastrap_phylips/.
    """
}

paramastrapPhylipsDir
    .flatMap {  id, dir -> dir.listFiles().collect { [id, it] } }
    .set { paramastrapPhylips }

/*
 *  Split the paramastrap alignment phylips into chunks of one alignment:
 *       For every item in paramastrapPhylips (contains a set of val(datasetID), file(multi_phylip_file)) :
 *           *) get the alignmentID from multi_phylip_file.baseName
 *           *) split file(multi_phylip_file) with splitPhylip
 *           *) create a new channel that emits a set containing val(datasetID), val(alignmentID), val(strapID), file(single_phylip)
 */

paramastrapPhylips
    .flatMap { set ->
        def datasetID = set[0]
        def alignmentID = set[1].baseName
        def file =  set[1]
        def strapID = 1
        splitPhylip(file).collect{ phylip -> [[datasetID, alignmentID, strapID],datasetID, alignmentID, strapID++, phylip] }
        
     }
     .map{ item -> 
         String stringStrap = item[0][2]; 
         [ tuple(item[0][0],item[0][1], stringStrap), item[1], item[2], item[3], item[4]]
     }  
     .set { splitPhylips  }

/*
 *
 **************************/

/**************************
 *
 *   C R E A T E   F I L E   O F   R E Q U I R E D   P A R A M A S T A R P   R E P L I C A T E S
 *
 */

randomAlignmentLists
    .concat(distantAlignmentLists)
    .set {concatAlignmentLists}

concatAlignmentLists.into { concatAlignmentLists1; concatAlignmentLists2}

process strap_trees_required {
    tag "strap trees required: ${datasetID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/strap_trees_required", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(numberOfAlignments), val(type), file(alignments_required) from concatAlignmentLists1

    output:
    set val(datasetID), val(numberOfAlignments), val(type), file("${datasetID}_${numberOfAlignments}_${type}_treesRequired.txt") into requiredStrapTrees, requiredStrapTrees2

    script:
    """
    while read p; do
      strapsNum=\$((${straps_num} / ${numberOfAlignments}))
      for ((i=1; i<=\$strapsNum; i++)); do
        echo "\$p \$i" >> ${datasetID}_${numberOfAlignments}_${type}_treesRequired.txt
       done
    done <${alignments_required}

    """
}

/*
 *
 **************************/

/**************************
 *
 *   C R E A T E   C H A N N E L   O F   R E Q U I R E D   P A R A M A S T R A P   R E P L I C A T E S
 *
 */

def splitListFile(file) {
    file.readLines().findAll {it}.collect {line -> line.tokenize(' ')}
}

requiredStrapTrees
  .map { set ->
    def datasetID = set[0]
    def file = set[3]
    splitListFile(file).collect { item -> tuple(datasetID, item[0], item[1]) }
  }
  .unique()
  .flatMap {item ->  item }
  .map { item -> [ tuple(item[0],item[1], item[2]), item[0], item[1], item[2] ] }
  .phase (splitPhylips) 
  .map {item -> [item[0][1], item[0][2], item[0][3], item[1][4]]}
  .set {requiredSplitPhylips}

/*
 *
 **************************/


/**************************
 *
 *   C R E A T E   P A R A M A S T R A P   T R E E   R E P L I C A T E S
 *
 */

process strap_trees {
    tag "bootstrap samples: ${datasetID} - ${alignmentID} - ${strapID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/strap_trees", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(alignmentID), val(strapID), val(phylip) from requiredSplitPhylips

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

/*
 *
 **************************/

/**************************
 *
 *   C R E A T E   L I S T   O F   D E F A U L T   A L I G N M E N T   S U P P O R T   V A L U E S
 *
 */

defaultAlignmentsBootstrapTrees.into { defaultAlignmentsBootstrapTreesA; defaultAlignmentsBootstrapTreesB }

defaultAlignmentsBootstrapTreesA
    .map {item -> [ item[0], item[3] ] }
    .groupTuple()
    .set { defaultAlignmentBootstrapFiles }

process default_support_tree_lists {
    tag "default_tree_list: Aligner ${y} ${datasetID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/supportSelection", mode: 'copy', overwrite: 'true'

    input:
    each y from ('clustalw', 'mafft', 'prank', 'tcoffee')
    set val (datasetID), file (defaultSupportTreesFiles) from defaultAlignmentBootstrapFiles

    output:
    set val (datasetID), file("${y}_supportFileList.txt") into defaultSupportFileLists
    set val (datasetID), file("${y}_concatenatedSupportTrees.nwk") into defaultConcatenatedSupportTrees
    set val (datasetID), val("${y}_concatenatedSupportTrees.nwk") into defaultSupportCombinationsTested

    script:
    """
    for ((i=0; i<$straps_num; i++)); do
       ls -1a ${datasetID}_${y}_\${i}.nwk >> ${y}_supportFileList.txt
    done

    while read p; do
      cat \$p >> ${y}_concatenatedSupportTrees.nwk
    done <${y}_supportFileList.txt
    """

}

/*
 *
 **************************/


/**************************
 *
 *   C R E A T E   L I S T   O F   A L T E R N A T I V E   A L I G N M E N T   S U P P O R T   V A L U E S
 *
 */


paramastrapTrees
    .groupTuple()
    .set { supportTreesFiles }

supportTreesFiles
    .cross(requiredStrapTrees2) 
    .view()
    .map {item -> [item[1][0], item[1][1], item[1][2], item[1][3], item[0][1]]}
    .set {supports}

process support_tree_lists {
    tag "distant and random tree_list: ${datasetID} - Alignments:${numberOfAlignments}"
    publishDir "${params.output}/${params.name}/${params.aligner}/$datasetID/supportSelection", mode: 'copy', overwrite: 'true'

    input:
    set val (datasetID), val(numberOfAlignments), val(type), file(treesRequired), file(fullListOfSupportTrees) from supports

    output:
    set val(datasetID), file("${type}_${numberOfAlignments}_supportFileList.txt") into SupportFileLists
    set val(datasetID), file("${type}_${numberOfAlignments}_concatenatedSupportTrees.nwk") into concatenatedSupportTrees
    set val(datasetID), val("${type}_${numberOfAlignments}_concatenatedSupportTrees.nwk") into supportCombinationsTested

    script:
    """
    while read p; do
      stringarray=(\$p)
      alingmentID=\${stringarray[0]}
      strapID=\${stringarray[1]}
      echo "\${alingmentID}_\${strapID}.nwk" >> ${type}_${numberOfAlignments}_supportFileList.txt
    done <${treesRequired}

    while read r; do
      cat \$r >> ${type}_${numberOfAlignments}_concatenatedSupportTrees.nwk
    done <${type}_${numberOfAlignments}_supportFileList.txt

    """
}

/*
 *
 **************************/



/**************************
 *
 *   G R O U P   P H Y L O G E N E T I C   T R E E S   B Y   D A T A S E T
 *
 *   Combine phlyogenetic trees to be grouped by ${datasetID}
 */

distantPhylogeneticTrees
    .concat (randomPhylogeneticTrees, defaultPhylogeneticTrees)
    .set { phylogeneticTrees }

/*
 *
 **************************/


/**************************
 *
 *   S Y N C  A N D   C O M B I N E   D A T A  B Y   D A T A S E T
 *
 */

concatenatedSupportTrees
    .concat(defaultConcatenatedSupportTrees)
    .groupTuple()
    .set{fullConcatenatedSupportTrees}

/*
 *  supportCombinationsTested <- set val(datasetID), val("${type}_${numberOfAlignments}_concatenatedSupportTrees.nwk")
 *  .collectFile emits -> datasetID.txt containing "> suppport1" etc 
 *  .map emits -> dataset, file(list_of_supports)
 *  
 *  fullConcatenatedSupportTrees <- set val(datasetID), file("${type}_${numberOfAlignments}_concatenatedSupportTrees.nwk")
 *
 *
 */

supportCombinationsTested
    .concat(defaultSupportCombinationsTested)
    .collectFile() { item ->
        ["${item[0]}", ">" + item[1] + '\n' ]
    }
    .map { item -> [ item.name, item ] }
    .phase(fullConcatenatedSupportTrees)
    .map { item -> [ item[0][0], item[0][1], item[1][1] ] }
    .cross (phylogeneticTrees)
    .map { item -> [ item[0][0], item[1][1], item[1][2], item[0][1], item[0][2] ] }
    .set { supportCombinationsTestedGrouped }

/*
 *
 **************************/


/**************************
 *
 *   E V A L U A T E   T R E E S   A N D   S U P P O R T S
 *
 */

process node_support {

    tag "node_support: ${datasetID} - ${treeID}"
    publishDir "${params.output}/${params.name}/${params.aligner}/${datasetID}/nodeSupport", mode: 'copy', overwrite: 'true'

    input:
    set val(datasetID), val(treeID), file(phylogeneticTree), file(supportFileList), file(listOfSupportTrees) from supportCombinationsTestedGrouped
    file(referenceTree) from refTree.first()

    output:
    set val(datasetID), val(treeID), file("nodeSupportFor_${datasetID}_${treeID}_Tree.result") into nodeSupports

    script:
    """
    t_coffee -other_pg seq_reformat \
             -in $phylogeneticTree \
             -in2 ${supportFileList} \
             -action +tree2ns ${referenceTree} \
             > nodeSupportFor_${datasetID}_${treeID}_Tree.result
    """
}

/*
 *
 **************************/


/**************************
 *
 *   F U N C T I O N S
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
 *
 **************************/

