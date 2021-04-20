from metagraph_workflows.resource_management import BuildResources, ResourceConfig
from metagraph_workflows import cfg_utils, constants, rule_utils

if build_primary:
    ruleorder: build_joint_primary > build
else:
    ruleorder: build > build_joint_primary


rule_name="build"
rule build:
    input: seqs_file_list_path
    output: graph_path
    threads: max_threads
    resources:
        mem_mb=BuildResources().get_mem(config),
    params:
        k=config['k'],
        tempdir_opt=cfg_utils.temp_dir_config(config),
        mem_cap=BuildResources().get_mem_cap(config),
        disk_cap=cfg_utils.get_rule_specific_config(rule_name, constants.DISK_CAP_MB_KEY, config),
        log_writing=get_log_opt(rule_name)
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

canonical_graphs_dir=wdir/f'{graph}_canonical'
canonical_graph_path=wdir/f'{graph}_canonical.dbg'

joint_contigs_path=wdir/f'{graph}_primary.fasta.gz'



sample_ids_spec = False
orig_samples_path=wdir/'orig_samples'

rule_name="stage_samples"
rule stage_samples:
    output: orig_samples_path/"{sample_id}.fasta.gz"
    params:
        staging_script_path=config[constants.SAMPLE_STAGING_SCRIPT_PATH]
    shell:
        """
        bash {params.staging_script_path} {wildcards.sample_id} {output}
        """


rule_name="build_canonical_graph_single_sample"
rule build_canonical_graph_single_sample:
    input: rule_utils.get_build_single_sample_input(config, orig_samples_path, seq_ids_dict)
    output: canonical_graphs_dir/"{sample_id}.dbg"
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(rule_name).get_mem(config),
    params:
        k=config['k'],
    shell:
        """
        echo "{input}" | {exec_cmd} build --parallel {threads} --mode basic -k {params.k} -o {output}
        """


rule_name="primarize_canonical_graph_single_sample"
rule primarize_canonical_graph_single_sample:
    input: canonical_graphs_dir/"{sample_id}.dbg"
    output: contigs_dir/"{sample_id}_primary.fasta.gz"
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(rule_name).get_mem(config),
    shell:
        """
        echo "{input}" | {exec_cmd} transform --to-fasta --primary-kmers --parallel {threads} -o {output}
        """


rule_name="build_joint_graph"
rule build_joint_graph:
    input: rule_utils.get_build_joint_input(config, contigs_dir, seq_ids_dict, seqs_file_list_path)
    output: canonical_graph_path
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(rule_name).get_mem(config),
    params:
        k=config['k'],
        separate_build=str(bool(config[constants.PRIMARIZE_SAMPLES_SEPARATELY])).lower(),
    shell:
        """
        if {params.separate_build}; then
            SEQ_PATHS={wdir}/seqs_paths.txt
            echo "{input}" | tr ' ' '\n' > $SEQ_PATHS
        else
            SEQ_PATHS="{input}"
        fi

        # TODO: canonical or basic?
        cat $SEQ_PATHS | {exec_cmd} build --parallel {threads} --mode basic -k {params.k} -o {output} 
        
        """

rule_name="primarize_joint_graph"
rule primarize_joint_graph:
    input: canonical_graph_path
    output: joint_contigs_path
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(rule_name).get_mem(config),
    shell:
        """
        echo "{input}" | {exec_cmd} transform --to-fasta --primary-kmers --parallel {threads} -o {output}
        """


rule_name="build_joint_primary"
rule build_joint_primary:
    input: joint_contigs_path
    output: graph_path
    threads: max_threads
    resources:
        mem_mb=ResourceConfig(rule_name).get_mem(config),
    params:
        k=config['k'],
    shell:
        """
        {exec_cmd} build --parallel {threads} --mode primary -k {params.k} -o {output} {input}
        """
