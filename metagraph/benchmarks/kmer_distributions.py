import matplotlib.pyplot as plt
import importlib
import scipy
# plot histograms

# statistical tests whether they are poisson?
# import the python file from build/
def evaluate_cpp_probability_distribution(dir):
    i = importlib.import_module(dir + "_distribution")
    i.results_exhaustive

def statistical_test(kmer_counts_row):
    result = scipy.stats.kstest(kmer_counts_row, "poisson")
    return result.pvalue, result.sign

def plot_distribution(kmer_counts_row, title):
    n, bins, patches = plt.hist(x=kmer_counts_row, bins='auto', color='#0504aa',
                                alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel("K-mer's occurance")
    plt.ylabel('Frequency')
    plt.title(title)
    pvalue, sign = statistical_test(kmer_counts_row)
    plt.suptitle("Kolmogorov-Smirnov: p-value "+ pvalue + ", sign: " + sign)
    plt.savefig(dir.replace(".", "\\") + "_dist.png")
    plt.show()

print("start")
folder = "C:/Users/myrth/Desktop/coding/metagraph/metagraph/build/"

i = importlib.import_module(folder + "kmer_dist_table.py")
counts_matrix = i.counts_matrix_in
counts_matrix = i.counts_matrix_out

random_rows = [2, 18, 20, 300, 240]# only sample even numbers
for i in random_rows:
    plot_distribution(counts_matrix_in[i], "in labels") # for in_labels
    plot_distribution(counts_matrix_out[i], "out labels") # for in_labels
    plot_distribution(counts_matrix_in[i] + counts_matrix_out[i], "union of in and out labels") # for in_labels



# read in the .fai file, and select the hepatitis.  'C:/Users/myrth/Desktop/coding/diff_assembly_counts/split/virus.shuf1.100.fa.fai'
#     f = open("C:/Users/myrth/Desktop/coding/metagraph/metagraph/build/virus_data/virus.shuf1.100.fa.fai","r")
#     lines = f.readlines()
#     labels_in = []
#     labels_out = []
#
#     for line in lines:
#         line = line.replace(",", "")
#         line_list = line.split("|")
#         #line_list = line_list.split("\t")
#         filename = line_list[3]
#         tag = line_list[4]
#         if tag[:8] == "Hepatitis"[:8]:
#             labels_in.append(filename + ".fa")
#         else:
#             labels_out.append(filename + ".fa")
#     print("\", \"".join(labels_in))
#     print("\", \"".join(labels_out))
#

# Once you have established this, the next question is whether the distribution is homogeneous or not. This means whether all parts of the data are handled by the same poisson distribution or is there a variation in this based on some aspect like time or space. Once you have convinced of these aspects, try the following three tests