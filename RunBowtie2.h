#ifndef RUNBOWTIE2_H
#define RUNBOWTIE2_H

#include "Header.h"
#include "Tools.h"

class RunBowtie2 {
private:
    iRead* reads;
    vector<interval>* depth;
    vector<intervalVirtualBowtie2> iDepth_virtual;
    string outputDir;

    vector<readNumAndSnap> rearrange(string);
    int isReverse(readNumAndSnap);
    int* flagAnalyzer(int);
    void addIntervalToVirtualDepth(long long, long long);
    void mergeDepthVirtualWithDepth(int);
    void mergeOnIntervalOfDepthVirtualwithDepth(long long, long long, long long);

public:
    RunBowtie2(iRead*, vector<interval>*, string);
    ~RunBowtie2();

    void run(int, string, int, int, int, string);

};

#endif // RUNBOWTIE2_H
