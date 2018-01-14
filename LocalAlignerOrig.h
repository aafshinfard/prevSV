#ifndef LOCALALIGNERORIG_H
#define LOCALALIGNERORIG_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>

class LocalAlignerOrig {

public:
    std::string mAlignmentSeqA;
    std::string mAlignmentSeqB;
    LocalAlignerOrig(std::string, std::string, int, int, int, double, int);
    void process();
    void backtrack();
    void Print();
    ~LocalAlignerOrig();
    int mScore, totalMismatch = 0, totalGap = 0, totalGapr = 0, totalGapfa = 0; //totalMatch = 0
    //std::string cigar;
    //std::string realCigar;
    //void produceCigar();
    int MATCH;
    int MISMATCH;
    int GAP;
    int counter;
    int shift;
    //int NegativeInf;
    //int k;

private:
    std::string mSeqA;
    std::string mSeqB;
    int** mD;
    int weight(size_t, size_t);
    //std::string reverseCigar();
};




#endif // LOCALALIGNERORIG_H
