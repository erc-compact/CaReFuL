#!/usr/bin/env nextflow
nextflow.enable.dsl=2



process filtool {
    label 'CaReFuL_filtool'
    container "${params.apptainer_images.pulsarx}"

    input:
    tuple path(root_dir), path(candFile), path(dmFile), val(filterbanksList), val(idx)

    output:
    tuple path(root_dir), path(candFile), path(dmFile), path("*.fil")

    beforeScript "eval \$(ssh-agent -s) && ssh-add ~/.ssh/hercules && if [ ! -f  /tmp/${idx}/TRAPUM_ARCHIVE/link.exists ]; then mkdir -p /tmp/${idx}; sshfs archive:/p/MFR/MEERKAT/TRAPUM /tmp/${idx}/ -o reconnect -o ro; fi"


    script:
    filterbankFiles = filterbanksList.join(" ")
    filterbankFiles = filterbankFiles.replaceAll("UUID", "${idx}")
    """
    #!/bin/bash
    workdir=\$(pwd)
    
    publish_dir="${root_dir}/01_FILTOOLED/"
    mkdir -p \${publish_dir}
    cd \${publish_dir}
    beam_name=${root_dir}
    filtool_prefix="tmp_\${beam_name}"
    filtool_tmp_output="tmp_\${beam_name}_01.fil"
    filtool_output="\${beam_name}_01.fil"
    if [ ! -f \${filtool_output} ]; then
        echo "Running filtool on \${filterbankFiles}"
        filtool -t ${task.cpus} --td ${params.filtool.time_decimation_factor} --fd ${params.filtool.freq_decimation_factor} --telescope ${params.filtool.telescope} -z ${params.filtool.rfi_filter} -o \$filtool_prefix ${params.filtool.extra_args} -f ${filterbankFiles} 
        mv \${filtool_tmp_output} \${filtool_output}
    else
        echo "Filtool output already exists, skipping filtool step"
    fi
    cd \${workdir}
    ln -s \${publish_dir}/\${filtool_output} \${filtool_output}
    fusermount -u /tmp/${idx}
    echo "Unmounted /tmp/${idx}"
    """
}

process get_fft_size {
    label 'CaReFuL_get_fft_size'
    container "${params.apptainer_images.presto}"

    input:
    path(filterbankFile)

    output:
    env(fft_size)

    script:
    """
    #!/bin/bash
    output=\$(readfile ${filterbankFile})
    echo "Output: \$output"

    value=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')

    log2=\$(python -c "import math; print(math.log(\$value, 2))")
    rounded_log2=\$(echo "\$log2" | awk -F"." '{if (\$2 >= 35) print \$1+1; else print \$1}')

    export fft_size=\$((2**\$rounded_log2))
    echo "nsamples: \$value"
    echo "log2: \$log2"
    echo "rounded_log2: \$rounded_log2"
    echo "FFT size: \$fft_size"

    """
}

process peasoup {
    label 'CaReFuL_peasoup'
    container "${params.apptainer_images.peasoup}"
    errorStrategy 'ignore'
    clusterOptions = '--gres=gpu:1'
    

    input:
    tuple path(root_dir), path(candFile), path(dmFile), path(filterbankFile), val(fft_size)

    output:
    tuple path(root_dir), path(candFile), path(filterbankFile), path("overview.xml")

    script:
    """
    #!/bin/bash

    # Accumulate optional arguments
    optional_args=""
    if [ "${params.peasoup.birdie_list}" != "null" ]; then
        optional_args="\${optional_args} -z ${params.peasoup.birdie_list}"
    fi

    if [ "${params.peasoup.chan_mask}" != "null" ]; then
        optional_args="\${optional_args} -k ${params.peasoup.chan_mask}"
    fi

    if [ "${params.peasoup.nsamples}" != "null" ]; then
        optional_args="\${optional_args} --nsamples ${params.peasoup.nsamples}"
    fi

    if [ "${params.peasoup.extra_args}" != "null" ]; then
        optional_args="\${optional_args} ${params.peasoup.extra_args}"
    fi

    max_acc=\$(cat ${candFile} | grep -v '#' | awk '{print \$3}' | tr -s "-" " " | sort -g | tail -n 1)
    acc_range=\$(echo \${max_acc} ${params.peasoup.acc_scale} | awk '{print \$1*\$2}')

    acc_start=\$(echo \${acc_range} -1 | awk '{print \$1*\$2}')
    acc_end=\$(echo \${acc_range} +1 | awk '{print \$1*\$2}')


    #cp /scratch/vkrishna/COMPACT/code/CaReFuL/nf/work/3c/a6b6560b2bce40c1e6fdb0aa3c8382/overview.xml . 
     

    peasoup -i ${filterbankFile} \
    --acc_start \${acc_start} \
    --acc_end \${acc_end} \
    --acc_tol ${params.peasoup.accel_tol} \
    -m ${params.peasoup.min_snr} \
    --ram_limit_gb ${params.peasoup.ram_limit_gb} \
    --nharmonics ${params.peasoup.nh} \
    --limit ${params.peasoup.total_cands} \
    --fft_size ${fft_size} \
    --start_sample ${params.peasoup.start_sample} \
    --cdm ${params.peasoup.coherent_dm} \
    \${optional_args} \
    --dm_file ${dmFile} \
    -o .
   
    """
}

process shortlist_folds {
    label 'CaReFuL_shortlist_folds'

    input:
    tuple path(root_dir), path(candFile), path(filterbankFile), path(xmlFile)

    output:
    tuple path(root_dir), path(candFile), path(filterbankFile), path("${root_dir}/04_FOLDED/output.candfile"), path("${root_dir}/04_FOLDED/pulsarx_cmds.txt")


    script:
    """
    #!/bin/bash
    mkdir -p  ${root_dir}/04_FOLDED/

    python ${baseDir}/scripts/shortlist_cand_create_candfile.py --input_candfile ${candFile} \
        --xml ${xmlFile}  --ncands ${params.shortlist_folds.ncands} --ptol ${params.shortlist_folds.ptol} \
        --dm_tol ${params.shortlist_folds.dm_tol} --template_dir ${baseDir}/include/ --cdm ${params.pulsarx.coherent_dm} --out_dir \$(realpath ${root_dir})/04_FOLDED/ -b \$(basename ${root_dir})

    """

}

process fold_with_pulsarx {
    label 'CaReFuL_pulsarx'
    container "${params.apptainer_images.pulsarx}"


    input:
    tuple path(root_dir), path(inputCandFile), path(filterbankFile), path(outputCandFile), path(commandsFile)
    
    output:
    tuple path("${root_dir}/04_FOLDED/direct_fold/*"), path("${root_dir}/04_FOLDED/search_fold/*"), path("pulsarx_cmds.txt"), path(inputCandFile), path(outputCandFile)


    script:
    """
    #!/bin/bash
    workdir=\$(pwd)
    mkdir -p  ${root_dir}/04_FOLDED/
    cd ${root_dir}/04_FOLDED/

    mkdir -p direct_fold search_fold

    while read -r line; do
    eval "\$line" &
    done < ${commandsFile}

    #get total list of psrfold_fil2 command running 

    njobs=\$(ps -ef | grep "psrfold_fil2" | grep -v grep | wc -l)
    echo "Total number of jobs running: \$njobs"

    #wait for all jobs to finish
    while [ \$njobs -gt 0 ]; do
        sleep 5
        njobs=\$(ps -ef | grep "psrfold_fil2" | grep -v grep |  wc -l)
        echo "Total number of jobs running: \$njobs"
    done
    echo "All jobs finished"

    """

}


workflow {

    Channel.fromPath("/scratch/vkrishna/COMPACT/FOLLOW_UP/3HM_FOLLOW_UP/2021-01-15-15:01:25/cfbf*/", type: 'dir')
            .map {
                dir ->
                def candFile = dir + "/input.candfile"
                def dmFile = dir + "/input.dmfile"
                def filterbankListFile = dir + "/input_fil.list"
                def filterbanksList = filterbankListFile.readLines()
                def UUID = UUID.randomUUID().toString()
                return tuple(
                    dir,
                    candFile,
                    dmFile,
                    filterbanksList,
                    UUID
                )
            }.set { inputFiles }


    def input = inputFiles
    def filtool_output = filtool(input)
    // append the fft size to the filtool output

   def output_filterbanks = filtool_output.map { tuple ->
        return tuple[3]
    }

    def fft_sizes = get_fft_size(output_filterbanks)
    def peasoup_input = filtool_output.merge(fft_sizes)

    def peasoup_output = peasoup(peasoup_input)

    def shortlist_folds_output = shortlist_folds(peasoup_output)

    def fold_output = fold_with_pulsarx(shortlist_folds_output)


   
    



    
   
    
// /scratch/vkrishna/COMPACT/code/CaReFuL/nf/work/3c/a6b6560b2bce40c1e6fdb0aa3c8382


    
}