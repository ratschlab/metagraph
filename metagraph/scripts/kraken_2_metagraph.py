import sys

def get_results_from_metagraph(filepath: str):
    lines = open(filepath).readlines()
    results = {}

    for line in lines:
        if line.strip() == "":
            continue
        if len(line.split(' ')) != 8:
            print("skip1:", line)
            continue
        if line.split(' ')[0] != "Sequence":
            print("skip2:", line)
            continue
        label = line.split(" ")[1]
        prediction = line.split(" ")[7].split("'")[1]
        results[label] = int(prediction)
    return results

def get_results_from_kraken(filepath: str):
    lines = open(filepath).readlines()
    results = {}

    not_classified = 0

    for line in lines:
        if len(line.split('\t')) != 5:
            print("skip1:", line)
            continue
        if line.split('\t')[0] == 'U':
            not_classified += 1
            continue

        if line.split('\t')[0] != 'C':
            print("skip2:", line)
            continue

        label = line.split('\t')[1]
        prediction = line.split('\t')[2]
        results[label] = int(prediction)
    return results, not_classified

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

def get_prediction_topology(curr_taxid: int, expected_taxid: int, taxo_root:int, taxo_parent:{int, int}, dist_to_root: int):
    dist_to_expected = 0
    aux_taxid = curr_taxid
    while aux_taxid != expected_taxid:
        if aux_taxid == taxo_root:
            break
        aux_taxid = taxo_parent[aux_taxid]
        dist_to_expected += 1
    if aux_taxid == expected_taxid:
        return dist_to_expected + dist_to_root

    aux_tax = expected_taxid
    dist_to_expected = 0
    while aux_tax != curr_taxid:
        if aux_tax == taxo_root:
            break
        aux_tax = taxo_parent[aux_tax]
        dist_to_expected += 1
    if curr_taxid == aux_tax:
        return dist_to_root - dist_to_expected
    return 18


def main(argv: [str]):
    path_to_metagraph_output = argv[1]
    path_to_kraken2_output = argv[2]
    path_to_lookup_label_taxid = argv[3]
    path_to_taxo_tree = argv[4]

    metagraph_results = get_results_from_metagraph(path_to_metagraph_output)
    print("before kraken results")
    kraken_results, kraken_fn = get_results_from_kraken(path_to_kraken2_output)

    print("len(metagraph_results)", len(metagraph_results))
    print("len(kraken_results)", len(kraken_results))
    print("kraken_fn", kraken_fn)

    assert len(metagraph_results) == len(kraken_results) + kraken_fn

    print("before lookup")
    lookup_label_taxid = get_lookup_label_taxid(path_to_lookup_label_taxid)

    print("lookup size", len(lookup_label_taxid))

    print("before taxo parents")
    taxo_parent, taxo_root = get_taxo_parents(path_to_taxo_tree)
    print("after taxo parents")

    print("taxo parent size", len(taxo_parent))
    print("taxo root", taxo_root)

    num_unclassified = 0

    import numpy

    real_vs_metagraph = numpy.zeros((20, 20))
    real_vs_kraken    = numpy.zeros((20, 20))

    not_in_metagraph = 0
    not_in_kraken = 0

    for label in metagraph_results:
        # acc_version = label.split("|")[3]
        acc_version = label.split("|")[2].split("-")[0]
        if acc_version not in lookup_label_taxid:
            num_unclassified += 1
            # print("acc vers not in lookup table", acc_version)
            continue
        expected_taxid = int(lookup_label_taxid[acc_version].strip())
        if expected_taxid not in taxo_parent:
            num_unclassified += 1
            # print("taxid not in taxotree:", expected_taxid)
            continue

        dist_to_root = 0
        aux_tax = expected_taxid
        while aux_tax != taxo_root:
            dist_to_root += 1
            aux_tax = taxo_parent[aux_tax]

        if dist_to_root > 15:
            print("bigger dist to root", dist_to_root)
            return 1


        curr_taxid = metagraph_results[label]
        if curr_taxid == 0:
            not_in_metagraph += 1
        elif curr_taxid not in taxo_parent:
            real_vs_metagraph[19][dist_to_root] += 1
        else:
            dist_predict = get_prediction_topology(curr_taxid=curr_taxid, expected_taxid=expected_taxid, taxo_root=taxo_root, taxo_parent=taxo_parent, dist_to_root=dist_to_root)
            real_vs_metagraph[dist_predict][dist_to_root] += 1

        label = label.split("'")[1]
        if label not in kraken_results:
            not_in_kraken += 1
            continue
        curr_taxid = kraken_results[label]
        #print("curr taxid kraken", curr_taxid)
        if curr_taxid not in taxo_parent:
            real_vs_kraken[19][dist_to_root] += 1
        else:
            dist_predict = get_prediction_topology(curr_taxid=curr_taxid, expected_taxid=expected_taxid, taxo_root=taxo_root, taxo_parent=taxo_parent, dist_to_root=dist_to_root)
            real_vs_kraken[dist_predict][dist_to_root] += 1

    total_predictions = len(metagraph_results)

    num_unclassified /= total_predictions
    not_in_metagraph /= total_predictions
    not_in_kraken /= total_predictions

    real_vs_metagraph /= total_predictions
    real_vs_kraken /= total_predictions

    #real_vs_metagraph_d = [[0.0 for x in range(8)] for y in range(10)]
    #real_vs_kraken_d    = [[0.0 for x in range(8)] for y in range(10)]
    #for x in range(9):
    #    for y in range(8):
    #        real_vs_metagraph_d[x][y] = 1.0 * real_vs_metagraph[x][y] / total_predictions
    #        real_vs_kraken_d[x][y] = 1.0 * real_vs_kraken[x][y] / total_predictions

    metagraph_pred_not_in_tree = 0.0
    kraken_pred_not_in_tree = 0.0

    metagraph_wrong_class = 0.0
    kraken_wrong_class = 0.0

    metagraph_tp = 0.0
    metagraph_fn = 0.0
    metagraph_vp = 0.0
    metagraph_fp = 0.0

    kraken_tp = 0.0
    kraken_fn = 0.0
    kraken_vp = 0.0
    kraken_fp = 0.0

    metagraph_dist_pred = numpy.zeros(40)
    print("real_vs_metagraph")

    print("Metagr", end="\t")
    for j in range(real_vs_metagraph.shape[1]):
        print(j, end="\t")
    print(" ")
    for i in range(real_vs_metagraph.shape[0]):
        print(i, end="\t")
        for j in range(real_vs_metagraph.shape[1]):
            print("%.4f" % real_vs_metagraph[i, j], end="\t")
            if i == 18:
                metagraph_wrong_class += real_vs_metagraph[i, j]
                continue
            if i == 19:
                metagraph_pred_not_in_tree += real_vs_metagraph[i, j]
                continue
            metagraph_dist_pred[i - j + 20] += real_vs_metagraph[i, j]
            if i >= j:
                metagraph_tp += real_vs_metagraph[i, j]
            else:
                metagraph_vp += real_vs_metagraph[i, j]
        print(" ")

    kraken_dist_pred = numpy.zeros(40)
    print("\n\n real_vs_kraken")
    print("kraken", end="\t")
    for j in range(real_vs_kraken.shape[1]):
        print(j, end="\t")
    print(" ")
    for i in range(real_vs_kraken.shape[0]):
        print(i, end="\t")
        for j in range(real_vs_kraken.shape[1]):
            print("%.4f" % real_vs_kraken[i, j], end="\t")
            if i == 18:
                kraken_wrong_class += real_vs_kraken[i, j]
                continue
            if i == 19:
                kraken_pred_not_in_tree += real_vs_kraken[i, j]
                continue
            kraken_dist_pred[i - j + 20] += real_vs_kraken[i, j]
            if i >= j:
                kraken_tp += real_vs_kraken[i, j]
            else:
                kraken_vp += real_vs_kraken[i, j]
        print(" ")


    print("\n num_unclassified (not in label tax map)", num_unclassified)
    print("metagraph pred not in tree", metagraph_pred_not_in_tree)
    print("kraken pred not in tree   ", kraken_pred_not_in_tree)
    print("not classified by metagraph", not_in_metagraph)
    print("not classified by kraken   ", not_in_kraken)

    print("metagraph wrong class:", metagraph_wrong_class)
    print("kraken wrong class:   ", kraken_wrong_class)

    print("\nmetagraph distances")
    for i in range(40):
        print(i-20, "\t",  "%.6f" % metagraph_dist_pred[i])

    print("\nkraken distances")
    for i in range(40):
        print(i-20, "\t", "%.6f" % kraken_dist_pred[i])

    print("metagraph genus acc:", metagraph_dist_pred[20] + metagraph_dist_pred[19])
    print("kraken genus acc:   ", kraken_dist_pred[20] + kraken_dist_pred[19])

    metagraph_fn = num_unclassified + metagraph_pred_not_in_tree + not_in_metagraph
    metagraph_fp = metagraph_wrong_class

    kraken_fn = num_unclassified + kraken_pred_not_in_tree + not_in_kraken
    kraken_fp = kraken_wrong_class

    print("\ntp = exact match or descendant")
    print("fn = failed to assign any class")
    print("vp = ancestor match")
    print("fp = classification that is not correct")

    print("metagr_tp", metagraph_tp)
    print("kraken_tp", kraken_tp)
    print("metagr_fn", metagraph_fn)
    print("kraken_fn", kraken_fn)
    print("metagr_vp", metagraph_vp)
    print("kraken_vp", kraken_vp)
    print("metagr_fp", metagraph_fp)
    print("kraken_fp", kraken_fp)

if __name__ == '__main__':
    main(sys.argv)
