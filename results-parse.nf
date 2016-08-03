params.input         = "$baseDir/tutorial/results/CLUSTALW/tips16_0.5_010.0400/"
params.output        = "$baseDir/tutorial/results/nodeSupport"

/*
 * Input parameters validation
 */

input_folder            = file(params.input)
results_path            = file(params.output)

/*
 * Create a channel for input file 
 */
 
datasets = Channel
    .fromPath( params.input )
    .ifEmpty { error "Cannot find any input folder matching: ${params.input}" }
    .map { file -> tuple( file.baseName, file ) }
 
process tree_lists {
    publishDir "$results_path/CLUSTALW/$datasetID/", mode: 'copy', overwrite: 'true'
  
    input:
    set val(datasetID), file (resultsDir) from datasets

    output:
    set val(datasetID), file('base_10_TreesConcat.nwk') into baseConcatenatedTrees
    file('1x10_TreesConcat.nwk') into ConcatenatedTrees_1x10
    file('10x1_TreesConcat.nwk') into ConcatenatedTrees_10x1


    script:
    """
    for file in ${resultsDir}/strap_trees/bootstrap_*.nwk;
    do 
        cat \$file >> base_10_TreesConcat.nwk
    done 
   
    containsElement () {
        local e
        for e in "\${@:2}"; do [[ "\$e" == "\$1" ]] && return 0; done
        return 1
    }


    perl $baseDir/bin/hhsearch_cluster.pl ${resultsDir}/alignments/alternativeMSA a3m 1 1_10_treesList.txt

    perl $baseDir/bin/hhsearch_cluster.pl ${resultsDir}/alignments/alternativeMSA a3m 10 10_1_treesList.txt

    while IFS=\\= read var1x10 value1x10; do
        vars1x10+=(\$var1x10)
        values1x10+=(\$value1x10)
    done < 1_10_treesList.txt

    while IFS=\\= read var10x1 value10x1; do
        vars10x1+=(\$var10x1)
        values10x1+=(\$value10x1)
    done < 10_1_treesList.txt

    ls ${resultsDir}/alignments/alternativeMSA | while read file; 
    do
    if containsElement \$file \${vars1x10[@]}
        then
            echo "Selected: \$file for 10 bootstraps"
            name=\${file%.fasta}
            cat ${resultsDir}/strap_trees/\${name}_strap_{0..9}.nwk >> 1x10_TreesConcat.nwk
        else
            echo "Skipping: \$file"
        fi
    done

    ls ${resultsDir}/alignments/alternativeMSA | while read file;
    do
    if containsElement \$file \${vars10x1[@]}
        then
            echo "Selected: \$file for 1 bootstraps"
            name=\${file%.fasta}
            cat ${resultsDir}/strap_trees/\${name}_strap_0.nwk >> 10x1_TreesConcat.nwk
        else
            echo "Skipping: \$file"
        fi
    done 


    """
}

