#ifndef RNASEQ_H
#define RNASEQ_H

#include "Sampling.h"
#include "Signalling.h"
#include "Anchoring.h"
#include "Assignment.h"
#include <time.h>
#include "RunBowtie2.h"
#include "Tools.h"
#include <iomanip>
#include <ctime>
//#include <pthread.h>
#include <thread>
#include <mutex>

using namespace std;

typedef struct params params_t;
struct ExtInterval{ //-/ extension interval
    FragLoci extFrom;
    FragLoci extTo;
};

struct singleBowtied{
    long long readIndex;
    int fragNum;
    long long snap;
    int flag;
    singleBowtied(long long rI, int fN, long long s, int f):readIndex(rI),fragNum(fN),snap(s),flag(f){}
};

class SVAnchor{
private:
    long long numberOfReads;
    string readName;
    int characterPerRow;
    string genomeName;
    string indexName;
    int chunkSize;
    int numberOfThread;
    int dashV;
    int numGap;
    int anchoringShift;
    //int assignmentShift;
    //bool isSlideNotShift;
    int readLength;
    string outputDir;

    //-/ Ameer Edits
    vector<int> threadNums;
    float misMacthRateInReads = 0.009;
    float indelRateInReads  = 0.001;
    vector<int> AlignmentScores; // temp
    string genomeFileName;
    //long lenVar; // variation in length for local alignment in extension step
    long localAlThreshold;
    long long genomeLength;
    int* extensionThreadIntervals;
    double indelShift;
    double editDistance;
    std::mutex mainMutex;
    std::mutex inClustChangesMutex;
    std::mutex extChangesMutex;
    std::mutex bGRChangeMutex;
    std::mutex bGLChangeMutex;
    std::mutex getGenomeMutex;
    pthread_mutex_t lock;

    bool checkBinGenome(long long , long long, bool, bool);

    void changeBinGenome(long long , long long, bool);

    void changeBinGenomeByAnchoring();

    void keepReadsForExtension();
    void breakReadsForExtension(bool);
    static void* extensionStepWrapper(void *);
    void* extensionStep(void* );
    void extensionStep2(ifstream* , long long * , int* , int *readsFileCount);
    string getGenome(long long, int);
    long long getGenomeLength();
    string revComplemACGT(string );
    int numberOfExtendableFrags(vector<string> , vector<ExtInterval>::iterator , int , int );
    char ifAlignable(long long startLoci, string readString, int flag, int threadNum, bool isToRight, int jumpsCount, bool onlyBGlocal, bool useBinGenome);
    int ifAlignable(long long startLoci, string readString, int flag, bool isToRight, int numGapCount, bool useBinGenome, int *shift);
    void inClusterExtThread(ifstream* , long long * , int* , int *readsFileCount);
    void inClusterExtension(int readsFileCount);
    bool findNewClusters(vector<string> FragsOfRead, int firstFrag , long long firstPos, int lastFrag, long long lastPos, int flag, long long readIndex);
    bool checkInBetween(vector<string> FragsOfRead, int firstFrag , long long firstPos, int lastFrag, long long lastPos, int flag, long long readIndex);
    void inExtendability(vector<string> ,long long ,long long ,int ,int ,int , long long , int, int );
    long long basesRemained(bool);
    void writeBinGenome();

    //void deleteAllJunk(string, int);
    //void printStatus();

    //void readReads(string);
    void writeReads(string);

    vector<vector<singleBowtied>> singleUniqes;
    void readSingleUniques();
    void pruneSingleUniques();
    void writeSingleUniques();
    //void writeReadsNext(string);
    //void readDepth(string);
    //void readDepthOnInterval(string);
    //void writeDepth(string);
    //void writeDepthOnInterval(string);
    //void writeGenes(string);

    //void createFASTQAligned();
    //void createFASTQNotAligend();

    int findCharacterPerRow(string);
    long long findNumberOfRead(string);


public:
    iRead* reads;
    //-/ Ameer edits
    void writeScores();
    vector<vector<ExtInterval>> allExtIntervals;
    vector<int> allExtIndices;
    vector<int> extIndices;
    int numExtendables = 0;


    vector<interval> depth;
    vector<gene> genes;
    vector<iReadNext> readsNext;

    SVAnchor(string, string, string, int, int, int, int, int, int, bool, int, string, long long);
    ~SVAnchor();

    void automation(bool);
};
struct params {
    SVAnchor* ptr;
    pthread_mutex_t mutex;
    pthread_cond_t done;
    int id;
};


#endif // RNASEQ_H
