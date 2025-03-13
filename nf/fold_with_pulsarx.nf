nextflow.enable.dsl=2


process refold_with_pulsarx {
    label 'refold_with_pulsarx'
    container  "${params.psrfold_singular_image}"
    memory '6 GB'
    cpus 4
    executor 'condor'
    cache 'lenient'
    maxForks 50
    errorStrategy 'ignore'



    input:
    tuple val(beam), val(dm), val(psrfold_cmd), val(start_frac), val(end_frac), val(additional_flags)

    output:
    stdout

    script:
    def beam_dm_str = String.format("%s_%.02f",beam, Double.parseDouble(dm))
    """
    #!/bin/bash
    mkdir ${params.psrfold_out_dir}/${beam_dm_str} -p
    cd ${params.psrfold_out_dir}/${beam_dm_str}

    L=\$(echo $start_frac $end_frac | awk '{print int(7200*(\$2-\$1)/64)}')

    beam_num=\$(echo $beam | awk -F_ '{print \$NF}')

    final_prefix="${beam_dm_str}_\$(echo $additional_flags | sed -e s/-//g | tr -s ' ' '_' | sed -e s/search//g)_${start_frac}_${end_frac}"
    $psrfold_cmd $additional_flags --frac ${start_frac} ${end_frac} -o \$final_prefix -L \$L -i \${beam_num} -t ${task.cpus}
    """

}

workflow  {
    def fold_candidate_channel = Channel.fromPath(params.psrfold_csv_file)
                            .splitCsv(header: true, sep:',')
                            .map{ row ->
                                def beam = row.beam
                                def dm = row.dm
                                def psrfold_cmd = row.psrfold_cmd

                                return tuple(beam, dm, psrfold_cmd) 
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

   def flags_ch = ["--nosearch", "--nof1search --nof0search", ""]
   //def flags_ch = ["--nosearch"]
   //segment_params_init_ch.view()

   def combined_channel = fold_candidate_channel.combine(segment_params_init_ch).combine(flags_ch)
   combined_channel.count().view { n -> "number of elements: ${n}" }
    refold_with_pulsarx(combined_channel)
    
    
}
   
 
