

process prepfold {
    label 'prepfold'
    container  "${params.presto_singularity_image}"
    memory '5 GB'
    cpus 2
    errorStrategy 'ignore'
    executor 'condor'
    cache 'lenient'



    input:
    tuple val(fold_command_id), val(prepfold_cmd), val(prefix), val(filterbank_file), val(start_frac), val(end_frac), val(additional_flags)

    output:
    stdout

    script:
    """
    #!/bin/bash
    mkdir ${params.out_dir}/$fold_command_id -p
    cd ${params.out_dir}/$fold_command_id
    fil_prefix=\$(basename $filterbank_file | sed -e 's/_01.fil//g'| sed -e 's/J0514-4002A_//g')
    final_prefix="\${fil_prefix}_${prefix}_\$(echo $additional_flags | sed -e s/-//g | tr -s ' ' '_')_${start_frac}_${end_frac}"
    $prepfold_cmd $additional_flags -start $start_frac -end $end_frac -o \$final_prefix $filterbank_file
    """

}


workflow  {
    def fold_candidate_channel = Channel.fromPath("${params.csv_file}")
                            .splitCsv(header: true, sep:',')
                            .map{ row ->
                                def fold_candidate_id = row.foldcand_hex
                                def prepfold_cmd = row.prepfold_cmd
                                def filterbank_file = row.fil_file
                                def prefix = String.format("%.02f", Double.parseDouble(row.dm))

                                return tuple(fold_candidate_id, prepfold_cmd, prefix, filterbank_file) 
                            }
    


    def planObj = new groovy.json.JsonSlurper().parseText(params.segmentation_plan)
    def segmentParamsList = planObj.segmentation.collectMany { segName, segData ->
            segData.chunks.collect { chunkIndex ->
                def start = segData.offset + chunkIndex * segData.fraction
                def end = segData.offset + (chunkIndex + 1) * segData.fraction
                if(end > 1.0) {
                    return null
                }
                return tuple(start, end)
            }
        }.findAll{it != null}

   def segment_params_init_ch = Channel.from(segmentParamsList)

   def flags_ch = ["-nosearch", "-fine", " "]
   //segment_params_init_ch.view()

   def combined_channel = fold_candidate_channel.combine(segment_params_init_ch).combine(flags_ch)
   combined_channel.count().view { n -> "number of elements: ${n}" }


    prepfold(combined_channel)
    
    
}