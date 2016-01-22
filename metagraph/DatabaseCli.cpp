// TODO implement some sort of command line interface?

#include <iostream>
#include <RocksdbImpl.cpp>

/**
 * Basic command line interface to the RocksDB database without having
 * to rely on the RocksDB commandline tools. Currently supports:
 *     - kmer queries.
 */
int main(int argc, char const ** argv) {

    if (argc != 3) {
        std::cout << "usage: database-cli <path/to/rocksdb> <kmer query>" << std::endl;
        return -1;
    }

    std::string path_to_rocksdb = argv[1];
    std::string kmer_query = argv[2];

    dbg_database::IDatabase *therocksdb = new dbg_database::RocksdbImpl(path_to_rocksdb);

    auto annotations = therocksdb->get_annotation(kmer_query);
    for (auto it = annotations.begin(); it != annotations.end(); ++it)
        std::cout << *it << std::endl;

    return 0;
}
