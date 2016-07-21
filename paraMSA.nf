
/* 
 * Main paraMSA-NF pipeline script
 *
 * @authors
 * Jia-Ming Chang <chang.jiaming@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Cedric Notredame <cedric.notredame@gmail.com>
 * Evan Floden <evanfloden@gmail.com> 
 */


params.name          = "Concatenated MSAs and Consensus Phylogentic Trees Analysis"
params.input         = "$baseDir/tutorial/data/*.fa"
params.output        = "$baseDir/tutorial/results"
params.straps        = 10
params.alignments    = 10
params.aligner       = "CLUSTALW"


log.info "p a r a M S A - N F  ~  version 0.1"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "input                  : ${params.input}"
log.info "output                 : ${params.output}"
log.info "straps                 : ${params.straps}"
log.info "alignments             : ${params.alignments}"
log.info "aligner                : ${params.aligner}"
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

    x=0
    mkdir -p paramastrap_phylips/alternativeMSA
    for i in "\${alternativeMSAsFASTA[@]}"
    do
        outfileBase=\${i%.fasta}
        esl-reformat phylip \$i > \$outfileBase.phylip
        echo -e "\$outfileBase.phylip\nR\n${straps_num}\nY\n\$seed\n" | seqboot
        mv outfile paramastrap_phylips/\${outfileBase}_strap.phylip
        mv paramastrap_phylips/\${outfileBase}_strap.phylip paramastrap_phylips/.
    done
    rm -rf paramastrap_phylips/alternativeMSA

    esl-reformat phylip ${baseAlignment} > base.phylip
    echo -e "base.phylip\nR\n${straps_num}\nY\n\$seed\n" | seqboot
    mv outfile bootstrap.phylip

    """
}


paramastrapPhylipsDir.flatMap {  id, dir -> dir.listFiles().collect { [id, it] } }
                     .view()
                     .set { paramastrapPhylips }

/*
 *  Split the strap alignment phylips into chunks of one alignment
 */

/* For every item in paramastrapPhylips (contains a set of val(datasetID), file(multi_phylip_file)) :
 *   get the alignmentID from multi_phylip_file.baseName
 *   split file(multi_phylip_file) with splitPhylip
 *   create a new channel that emits a set containing val(datasetID), val(alignmentID), val(strapID), file(single_phylip)
 *
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

splitPhylips = Channel.create()
paramastrapPhylips
    .map { set ->
        def datasetID = set[0]
        def alignmentID = set[1].baseName
        def file =  set[1]
        splitPhylip(file).eachWithIndex{ phylip, strapID -> splitPhylips << tuple(datasetID, alignmentID, strapID, phylip) } 
     }


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

    raxmlHPC-SSE3 -m PROTGAMMAJTT -p 4533 -T 1 -s ${alignmentID}_${strapID}.phylip -n tmp
    cp RAxML_bestTree.tmp ${alignmentID}_${strapID}.nwk
    rm *.tmp

    """
}


