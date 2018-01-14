#ifndef ALIGNING_H
#define ALIGNING_H

#include "Header.h"
#include "Tools.h"
#include "LocalAligner.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

class Aligning {
private:
    iRead* reads;
    int characterPerRow;
    string genomeFileName;
    int chunkSize;
    long long numberOfReads;

    string getGenome(long long, int);
    string complement(string);


public:
    Aligning(iRead*, int, long long, string, int);
    ~Aligning();

    void breakReadsForThreads(int);
    void combineSamForThreads(int, bool);
    void prepareHeaderFile();
    void executeLocalAlignment(int);
};


#endif // ALIGNING_H
