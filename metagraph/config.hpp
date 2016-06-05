#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

#include <cstring>
#include <string>
#include <vector>

class Config {
    public:
        Config();
        Config(int argc, const char *argv[]);

        ~Config();
       
        void print_usage(std::string prog_name, int identity = -1);
        //void print_call(string prog_name);

        bool verbose;
        bool quiet;
        bool integrate;
        bool print_graph;
        
        unsigned int k;
        unsigned int distance;

        std::vector<std::string> fname;
        std::string outfbase;
        std::string infbase;
        std::string sqlfbase;

        enum identities {noidentity = -1, 
                         build = 1, 
                         merge, 
                         compare, 
                         align,
                         stats};
        int identity;

    private:
        void init();
};
#endif
