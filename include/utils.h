#ifndef UTILS_H
#define UTILS_H

#include <sys/stat.h>
#include <unistd.h>
#include <cstdlib>
#include <string>

#include "hooks.h"

// DataType definitions
using CSROrdinal = unsigned int;
using Value = double;

inline bool FileExist( const std::string& name ) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

class Config {
    public:
        std::string A;
        std::string B;
        std::string C;
        bool TrackIndividualTime;
        size_t InterationsWarmUp;
        size_t InterationsExecution;

        Config(): A(""), B(""), C(""), TrackIndividualTime(false), InterationsWarmUp(1), InterationsExecution(10) {}
};

void printUsage(const int argc, char * const argv[]) {
    std::cout << "usage: "<< std::endl;
    std::cout << "  -h print the help message" << std::endl;
    std::cout << std::endl;
    std::cout << "  -A input file for matrix A" << std::endl;
    std::cout << "  -B input file for matrix B, use A if not given" << std::endl;
    std::cout << "  -C output file for matrix C, don't write to output if not given" << std::endl;
    std::cout << "  -t set true for printing the runtime for each step, false if not set" << std::endl;
    std::cout << "  -w the number of iterations for warm up" << std::endl;
    std::cout << "  -e the number of iterations for execution" << std::endl;
    std::cout << std::endl;
    std::cout << "  Usage example: " << argv[0] << " -A ./exampleData/gre_185_l.mtx -B ./exampleData/gre_185_r.mtx -t -w 1 -e 2" << std::endl;
}

Config parseCommandLine(const int argc, char * const argv[])
{
    Config c;
    int ret;
    while ((ret = getopt(argc, argv, "hA:B:C:tw:e:")) != -1) {
        switch (ret) {
            case 'h':
                printUsage(argc, argv);
                std::exit(1);
            case 'A':
                c.A = std::string(optarg);
                break;
            case 'B':
                c.B = std::string(optarg);
                break;
            case 'C':
                c.C = std::string(optarg);
                break;
            case 't':
                c.TrackIndividualTime = true;
                break;
            case 'w':
                c.InterationsWarmUp = atol(optarg);
                break;
             case 'e':
                c.InterationsExecution = atol(optarg);
                break;
             default:
                std::cerr << "Error: the input option is invalid" << std::endl;
                printUsage(argc, argv);
                std::exit( 1 );
        }
    }

    if( c.A == "" ) {
        std::cerr << "Error: the input file is not given" << std::endl;
        printUsage(argc, argv);
        std::exit(1);
    }

    if( !FileExist( c.A ) ) {
        std::cerr << "Error: file A Not Found: " << c.A << std::endl;
        std::exit( 1 );
    }

    if( c.B == "" || !FileExist( c.B ) ) {
        c.B = c.A;
        std::cout << "Note: File B is not given or doesn't exist, so use A as B as well." << std::endl;
    }

    return c;
}
#endif
