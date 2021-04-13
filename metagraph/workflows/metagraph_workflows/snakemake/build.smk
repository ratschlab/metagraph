from pathlib import Path

# TODO: is this a good mechanism?
if build_primary:
    ruleorder: build_primary > build
else:
    ruleorder: build > build_primary

from metagraph_workflows.resource_management import BuildResources
from metagraph_workflows import cfg_utils, constants

rule_name='build'
rule build:
    input: seqs_file_list_path
    output: graph_path
    threads: max_threads 
    params:
        k=config['k'],
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildResources().get_mem_cap(config),
        disk_cap=cfg_utils.get_rule_specific_config(rule_name, constants.DISK_CAP_MB_KEY, config),
        log_writing=get_log_opt(rule_name)
    resources:
        mem_mb=BuildResources().get_mem(config),
    shell:
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