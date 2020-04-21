#include <spdlog/sinks/stdout_color_sinks.h>

#include "common/logger.hpp"
#include "common/algorithms.hpp"
#include "cli/config/config.hpp"
#include "cli/build.hpp"
#include "cli/annotate.hpp"
#include "cli/stats.hpp"
#include "cli/augment.hpp"
#include "cli/clean.hpp"
#include "cli/merge.hpp"
#include "cli/align.hpp"
#include "cli/query.hpp"
#include "cli/assemble.hpp"
#include "cli/server.hpp"
#include "cli/transform_graph.hpp"
#include "cli/transform_annotation.hpp"

using mg::common::logger;


int main(int argc, char *argv[]) {
    auto config = std::make_unique<Config>(argc, argv);

    logger->set_level(utils::get_verbose() ? spdlog::level::trace : spdlog::level::info);
    //logger->set_pattern("%^date %x....%$  %v");
    //spdlog::set_pattern("[%H:%M:%S %z] [%n] [%^---%L---%$] [thread %t] %v");
    //console_sink->set_color(spdlog::level::trace, "\033[37m");
    spdlog::flush_every(std::chrono::seconds(1));

    logger->trace("Metagraph started");

    switch (config->identity) {
        case Config::BUILD:
            return build_graph(config.get());

        case Config::EXTEND:
            return augment_graph(config.get());

        case Config::ANNOTATE:
            return annotate_graph(config.get());

        case Config::ANNOTATE_COORDINATES:
            return annotate_graph_with_genome_coordinates(config.get());

        case Config::MERGE_ANNOTATIONS:
            return merge_annotation(config.get());

        case Config::QUERY:
            return query_graph(config.get());

        case Config::SERVER_QUERY:
            return run_server(config.get());

        case Config::COMPARE:
            return compare(config.get());

        case Config::CONCATENATE:
            return concatenate_graph_chunks(config.get());

        case Config::MERGE:
            return merge_graph(config.get());

        case Config::CLEAN:
            return clean_graph(config.get());

        case Config::STATS:
            return print_stats(config.get());

        case Config::TRANSFORM_ANNOTATION:
            return transform_annotation(config.get());

        case Config::TRANSFORM:
            return transform_graph(config.get());

        case Config::ASSEMBLE:
            return assemble(config.get());

        case Config::RELAX_BRWT:
            return relax_multi_brwt(config.get());

        case Config::ALIGN:
            return align_to_graph(config.get());

        case Config::NO_IDENTITY:
            assert(false);
    }

    return 0;
}
