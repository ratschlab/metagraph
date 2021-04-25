import sys
import subprocess
from subprocess import PIPE


def get_lookup_label_taxid(path_to_lookup_label_taxid: str) -> {str, str}:
    lookup_label_taxid = {}
    lines = open(path_to_lookup_label_taxid).readlines()
    for line in lines[1:]:
        lookup_label_taxid[line.split("\t")[1].strip()]=line.split("\t")[2].strip()
    return lookup_label_taxid

def get_taxo_parents(path_to_taxo_tree: str) -> ({int, int}, int):
    taxo_parent = {}
    taxo_root = 0
    taxo_lines = open(path_to_taxo_tree).readlines()
    for line in taxo_lines:
        act_node = int(line.split('\t')[0])
        act_parent = int(line.split('\t')[2])
        taxo_parent[act_node] = act_parent
        if act_node == act_parent:
            taxo_root = act_node
    return taxo_parent, taxo_root

def is_descendant(taxo_parent:{int, int}, taxo_root:int, node:str, query:str):
    node = int(node)
    query = int(query)
    while query != taxo_root:
        if not query in taxo_parent: # this should not be here
            print("query", query, "taxo_root", taxo_root)
            return False
        query = taxo_parent[query]
        if query == node:
            return True
    return False

def main(argv: [str]):
    path_to_dbg = argv[1]
    path_to_taxo_db = argv[2]
    path_to_queries = argv[3]
    path_to_taxo_tree = argv[4]
    path_to_lookup_label_taxid = argv[5]
    lca_coverage = 0.66
    print("before lookup")
    lookup_label_taxid = get_lookup_label_taxid(path_to_lookup_label_taxid)
    print("before taxo parents")
    taxo_parent, taxo_root = get_taxo_parents(path_to_taxo_tree)
    print("after taxo parents")

    tax_class_command = '{exe} tax_class -i {dbg} {fasta_queries} --taxonomic-tree {taxoDB} \
                            --lca-coverage-fraction {lca_coverage} -p 10'.format(
        exe="./metagraph",
        dbg=path_to_dbg,
        fasta_queries=path_to_queries,
        taxoDB=path_to_taxo_db,
        lca_coverage=lca_coverage
    )
    res = subprocess.run([tax_class_command], shell=True, stdout=PIPE, stderr=PIPE)
    print("finished tax_class")
    if res.returncode != 0:
        print(res)
        exit(1)
    res_lines = res.stdout.decode().rstrip().split('\n')

    with open('tmp.txt', 'w') as f:
        for line in res_lines:
            print(line, file=f)


    percent_perfect_match = 0
    percent_ancestor_match = 0
    percent_descendant_match = 0
    percent_fail = 0
    percent_root_match = 0
    nonexistent_label = 0
    not_classified = 0

    num_lines = 0

    saved_lines = {}

    for line in res_lines:
        if line == "":
            continue
        print("line->" + line)
        num_lines += 1
        if num_lines % 10000 == 0:
            print("processing", num_lines)
            print(line)
        label = line.split(" ")[1].split("|")[3]

        if label not in lookup_label_taxid:
            nonexistent_label += 1
            continue

        expected_taxid = lookup_label_taxid[label]
        predicted_taxid = line.split(" ")[7].split("'")[1]

        if int(predicted_taxid) == 0:
            not_classified += 1
            continue

        if int(expected_taxid) == int(predicted_taxid):
            percent_perfect_match += 1
            continue
        if int(predicted_taxid) == taxo_root:
            percent_root_match += 1
        if int(predicted_taxid) not in taxo_parent:
            inexistent_label += 1
            continue

        if is_descendant(taxo_parent, taxo_root, expected_taxid, predicted_taxid):
            percent_descendant_match += 1
            continue
        if is_descendant(taxo_parent, taxo_root, predicted_taxid, expected_taxid):
            percent_ancestor_match += 1
            continue
        percent_fail += 1

    percent_perfect_match /= num_lines
    percent_ancestor_match /= num_lines
    percent_descendant_match /= num_lines
    percent_fail /= num_lines
    percent_root_match /= num_lines
    nonexistent_label /= num_lines
    not_classified /= num_lines

    print("percent_perfect_match", "%.5f" %  percent_perfect_match)
    print("percent_root_match", "%.5f" % percent_root_match)
    print("percent_ancestor_match", "%.5f" % percent_ancestor_match)
    print("percent_descendant_match", "%.5f" % percent_descendant_match)
    print("percent_fail", "%.5f" % percent_fail)
    print("nonexistent_label", "%.5f" % nonexistent_label)
    print("not_classified", "%.5f" % not_classified)

if __name__ == '__main__':
    main(sys.argv)
