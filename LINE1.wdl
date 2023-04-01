version 1.0

workflow line1 {
    input {
    	String output_directory
        File bam_file
        File bam_index_file
        #line1 workflow parameters
        Boolean paired_end = true
        String read1_range = 'L1HS:1-20'
        String read2_range = 'L1HS'
        String unpaired_range = 'L1HS:1-150'
        Int max_edit_distance = 2
        Int clip_max = 20
        Float downsample_to = 1.0
        #general parameters
        Int cpu = 8
        String memory = "64G"
        String docker = "mparikhbroad/line1_workflow:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call run_line1_counting {
        input:
            output_dir = output_directory_stripped,
            bam_file = bam_file,
            bam_index_file = bam_index_file,
            paired_end = paired_end,
            read1_range = read1_range,
            read2_range = read2_range,
            unpaired_range = unpaired_range,
            max_edit_distance = max_edit_distance,
            clip_max = clip_max,
            downsample_to = downsample_to,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    output {
        File line1_counts_file = run_line1_counting.counts_file
    }
}

task run_line1_counting {

    input {
        String output_dir
        File bam_file
        File bam_index_file
        Boolean paired_end
        String read1_range
        String read2_range
        String unpaired_range
        Int max_edit_distance
        Int clip_max
        Float downsample_to
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        eval "$(micromamba shell hook --shell=bash)"
        micromamba activate base
        python <<CODE
        from subprocess import check_call
        import os

        paired_end = ~{true="True" false="False" paired_end}
        read1_range = "~{read1_range}"
        read2_range = "~{read2_range}"
        unpaired_range = "~{unpaired_range}"
        max_edit_distance = ~{max_edit_distance}
        clip_max = ~{clip_max}
        downsample_to = ~{downsample_to}

        bam = pysam.AlignmentFile("~{bam_file}",'rb')

        if(paired_end):
            r1_UMI_dict = dict()
            for read in bam.fetch(region=read1_range):
                if random.random() < downsample_to and read.is_read1 and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and not read.is_reverse and read.cigarstring[-1]!='S' and 'N' not in read.cigarstring:
                    if read.cigartuples[0][0] != 4 or (read.cigartuples[0][1] <= clip_max):
                        CB = read.get_tag("CB")
                        UB = read.get_tag("UB")
                        if CB not in r1_UMI_dict:
                            r1_UMI_dict[CB] = set()
                        r1_UMI_dict[CB].add(UB)

            r2_UMI_dict = dict()
            for read in bam.fetch(region=read2_range):
                if read.is_read2 and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and read.is_reverse and 'S' not in read.cigarstring and 'N' not in read.cigarstring:
                    CB = read.get_tag("CB")
                    UB = read.get_tag("UB")
                    if CB not in r2_UMI_dict:
                        r2_UMI_dict[CB] = set()
                    r2_UMI_dict[CB].add(UB)
                    
            UMI_dict = dict()
            for CB in r1_UMI_dict:
                if CB in r2_UMI_dict:
                    UMI_dict[CB] = r1_UMI_dict[CB].intersection(r2_UMI_dict[CB])
        else:
            UMI_dict = dict()
            for read in bam.fetch(region=unpaired_range):
                if random.random() < downsample_to and read.has_tag("nM") and read.get_tag("nM") <= max_edit_distance and read.has_tag("CB") and read.has_tag("UB") and read.is_reverse and 'N' not in read.cigarstring:
                    if (read.cigartuples[0][0] != 4 or (read.cigartuples[0][1] <= clip_max)) and (read.cigartuples[-1][0] != 4 or (read.cigartuples[-1][1] <= clip_max)):
                        CB = read.get_tag("CB")
                        UB = read.get_tag("UB")
                        if CB not in UMI_dict:
                            UMI_dict[CB] = set()
                        UMI_dict[CB].add(UB)

        count_df = pd.DataFrame(index=list(UMI_dict.keys()), columns=['counts'])
        count_df['counts'] = [len(val) for val in UMI_dict.values()]
        count_df.to_csv('LINE1_counts.csv')
        
        CODE

        gsutil -m cp LINE1_counts.csv ~{output_dir}/
    >>>

    output {
        File counts_file = 'LINE1_counts.csv'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(bam_file, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}