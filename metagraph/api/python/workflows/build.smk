from pathlib import Path

# TODO: is this a good mechanism?
if build_primary:
    ruleorder: build_primary > build
else:
    ruleorder: build > build_primary

mem_cap_factor=0.9
disk_cap_factor=0.9

def get_mem(key):
    return config[key] if key in config else max_memory_mb

def get_disk(key):
    return config[key] if key in config else config['default_disk_mb']

def get_mem_cap(wildcards, resources):
    return max(int(resources.mem_mb * mem_cap_factor/1024), 1)

def get_disk_cap(wilcards, resources):
    return max(int(resources.disk_mb * disk_cap_factor/1024), 1)

temp_dir_config = f"--disk-swap {config['tmpdir']}" if 'tmpdir' in config else '',


rule build:
    input: input_files_list_path
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
    input: input_files_list_path
    output: canonical_graph_path
    threads: max_threads
    params:
        k=config['k'],
    shell:
        """
        cat {input} | {exec_cmd} build --parallel {threads} --canonical -k {params.k} -o {output}
        """

rule primarize_canonical_graph:
    input: canonical_graph_path
    output: wdir/f'{graph}_primary.fasta.gz'
    threads: max_threads
    shell:
        """
        echo "{input}" | {exec_cmd} transform --to-fasta --primary-kmers --parallel {threads} --canonical -o {output}
        """


rule build_primary:
    input: wdir/f'{graph}_primary.fasta.gz'
    output: graph_path
    threads: max_threads
    params:
        k=config['k'],
    shell:
        """
        {exec_cmd} build --parallel {threads} -c -k {params.k} -o {output} {input}
        """