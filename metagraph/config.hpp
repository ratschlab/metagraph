#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

    struct CFG 
    {
        //String<char> > fname;
        std::vector<std::string> fname;
        std::string outfbase;
        std::string infbase;
        std::string sqlfbase;
        std::string merge;
        std::string compare;
        bool verbose;
        bool integrate;
        unsigned int k;
        std::string db_connect_string;

        CFG() :
            verbose(false), k(31), db_connect_string("")
        {}
    };




#endif
