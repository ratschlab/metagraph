//
// Created by studenyj on 7/8/19.
//

//#pragma omp parallel
//{
// vector<char> routing_table_array_local;
//#pragma omp for
// for (int64_t node = 1; node <= this->graph.num_nodes(); node++) {
// routing_table_array_local.push_back('#');// to always start a block with #
// if (PathDatabaseDynamicCore<DRT, DIT>::node_is_split(node)) {
// auto &dynamic_table = PathDatabaseDynamicCore<DRT, DIT>::routing_table;
// alt_assert(dynamic_table.size(node));
// int encoded = 0;
// for (int64_t i = 0; i < dynamic_table.size(node); i++) {
// routing_table_array_local.push_back(dynamic_table.get(node, i));
// encoded++;
//}
// alt_assert(encoded);
//}
//}
//
// for (int t = 0; t < omp_get_num_threads(); t++) {
//#pragma omp barrier
// if (t == omp_get_thread_num()) {
// routing_table_array.encode(routing_table_array[i]);(routing_table_array.end(),all(routing_table_array_local));
//}
//}
//
//};

#include <sdsl/int_vector_buffer.hpp>
#include <string>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char *argv[]) {
    srand(time(NULL));
    auto dir = getenv("TMPDIR");
    string tmp_dir = dir ? dir : "/tmp";
    string tmp_file = tmp_dir + "/" + "tmp_local_"s + to_string(rand()) + ".bin";
    size_t size = 10000000;

    // create an int_vector_buffer
    int_vector_buffer<> test(
            tmp_file, // filename
            std::ios::out, // we do not want to open an existing file, but create a new one
            1024 * 1024, // use a buffer of about 1MB
            2, // use 64bit for each integer
            false); // use int_vector format

    // write sequentially random values to disk
    for (uint64_t i = 0; i < size; ++i) {
        for (int j = 0; j < 100; i++, j++) {
            test.push_back(0);
        }
        test.push_back(1);
    }

    memory_monitor::start();

    auto start = timer::now();
    construct(cst, argv[1], 1);
    auto stop = timer::now();
    cout << "construction cst time in seconds: "
         << duration_cast<seconds>(stop - start).count() << endl;

    memory_monitor::stop();

    std::cout << "peak usage = " << memory_monitor::peak() / (1024 * 1024) << " MB"
              << std::endl;

    std::ofstream cstofs("cst-construction.html");
    cout << "writing memory usage visualization to cst-construction.html\n";
    memory_monitor::write_memory_log<HTML_FORMAT>(cstofs);
    cstofs.close();
    util::clear(cst);
    test.close(true); // close buffer and (3) remove the file

    return 0;
}