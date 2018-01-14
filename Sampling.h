#ifndef SAMPLING_H
#define SAMPLING_H

#include "Header.h"
#include "Tools.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>


class Sampling {
private:
    iRead* reads;
    vector<interval>* depth;
    string genomeFileName;
    int characterPerRow;
    int floatingEdge;
    string outputDir;

    string getGenome(long long, long long);

public:
    Sampling(iRead*, vector<interval>*, int, string, int, string);
    ~Sampling();

    string revComplemACGT(string);
    void prepareSampling(string, long long, long long);
    long long prepareSampling(string, bool, int );
    void prepareExonSequence();
    void createIndexBowtie();
};

#endif // SAMPLING_H
