from pathlib import Path

# TODO: is this a good mechanism?
if build_primary:
    ruleorder: build_primary > build
else:
    ruleorder: build > build_primary


temp_dir_config = f"--disk-swap {config['tmpdir']}" if 'tmpdir' in config else '',


rule build:
    input: seqs_file_list_path
    output: graph_path
    threads: max_threads 
    params:
        #input_files=lambda wildcards, input: '\n'.join(input),
        k=config['k'],
        tempdir_opt=temp_dir_config,
        mem_cap=get_mem_cap,
        disk_cap=get_disk_cap,
        log_writing=get_log_opt("build")
    resources:
        mem_mb=get_mem('build_memory_mb'),
        disk_mb=get_disk('build_disk_mb'),
    shell:
         # TODO: --swap-dir
        """
        cat {input} | {exec_cmd} build \
        --parallel {threads} \
        -k {params.k} \
        -o {output} \
        --mem-cap-gb {params.mem_cap} \
        --disk-cap-gb {params.disk_cap} \
        {verbose_opt} \
        {params.tempdir_opt} 2>&1 {params.log_writing}
        """



### Build Primary

canonical_graph_path=wdir/f'{graph}_canonical.dbg'

rule build_canonical_graph:
    input: seqs_file_list_path
    output: canonical_graph_path
    threads: max_threads
    params:
        k=config['k'],
    shell:
        """
        cat {input} | {exec_cmd} build --parallel {threads} --mode basic -k {params.k} -o {output}
        """

rule primarize_canonical_graph:
    input: canonical_graph_path
    output: wdir/f'{graph}_primary.fasta.gz'
    threads: max_threads
    shell:
        """
        echo "{input}" | {exec_cmd} transform --to-fasta --primary-kmers --parallel {threads} -o {output}
        """


rule build_primary:
    input: wdir/f'{graph}_primary.fasta.gz'
    output: graph_path
    threads: max_threads
    params:
        k=config['k'],
    shell:
        """
        {exec_cmd} build --parallel {threads} --mode primary -k {params.k} -o {output} {input}
        """