#ifndef ANCHORING_H
#define ANCHORING_H

#include "Header.h"
#include "Tools.h"

class Anchoring {
private:
    iRead* reads;
    vector<interval>* depth;
    int d;
    int numGap;
    int numberFiles;
    string referenceName;
    string outputName;
    bool isSample;
    int chunkSize;
    int floatingEdge;
    int dashV;
    int numberOfThread;
    string outputDir;

    int* getMyArraySteps(int);
    int* getFragsOrders(int);

    void checkCorrespond(int, int, int, int, int bowtiedindex = -1 );
    void alignThenAnalyze(int);
    void alignThenAnalyzeFirstAndLast();
    void changeFastQ(int);
    vector<readNumAndSnap> rearrange(string);
    vector<int> getWhichCheck(int*, int);
    int* flagAnalyzer(int);
    bool isInTheSameDirection(readNumAndSnap, readNumAndSnap);

    int isReverse(readNumAndSnap);
    inline void changeBinGenome(long long , long long);

    void addIntervalToDepth(long long, long long);
    tuple<long long, long long> calculateSnapUnSample(long long, long long);
    tuple<bool,long long, long long, long long> calculateIntronOrExonAndDistance(long long);

public:
    int correspond;
    int correspondCluster2;
    vector< vector<readNumAndSnap> > notAligned;

    Anchoring(iRead*, vector<interval>*, string, string, bool, int, int, int, int, int, int, string);
    long long getNumberFiles();
    ~Anchoring();
    // OLDER FUNCTIONS
    void prepareFiles();
    void nextStepAlignment();
    void createNewReads(int);
    void writingNotAligned();

    // NEWNER FUNCTIONS
    void prepareFiles(int range1 ,int range2 ,int readsFilesCount ,int slide ,int iteration);
    void nextStepAlignment(int range1, int range2, int readsFilesCount);
    void moreFragClusters_Recursive(vector<readNumAndSnap> ,int ,int ,int ,int );
    inline int correctFragNo(int  , int );
    void writeSingleUniques(vector<readNumAndSnap> ReadVec ,int bowtiedIndex, ofstream* file, long long readIndex);
    bool checkCorrespond1read(vector<readNumAndSnap> ,int , int , int , int  );
    void RecursiveRemainedAligner(int,int);

};

#endif // ANCHORING_H
