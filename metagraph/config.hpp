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
       
        void print_usage(std::string prog_name);
        //void print_call(string prog_name);

        bool verbose;
        bool integrate;
        
        unsigned int k;
        unsigned int distance;

        std::vector<std::string> fname;
        std::string outfbase;
        std::string infbase;
        std::string sqlfbase;
        std::string merge;
        std::string compare;
        std::string align;

    private:
        void init();
};
#endif
