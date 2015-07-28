#ifndef __CONFIG_HPP__
#define __CONFIG_HPP__

    struct CFG 
    {
        //String<char> > fname;
        std::vector<std::string> fname;
        bool verbose;
        unsigned int k;

        CFG() :
            verbose(false), k(31)
        {}
    };




#endif
