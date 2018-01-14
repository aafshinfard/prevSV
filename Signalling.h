#ifndef SIGNALLING_H
#define SIGNALLING_H
#include <vector>
#include <Header.h>
#include "Tools.h"
#include <unistd.h>


struct EventInterval{
    long long start = 0;
    long long end = 0;
};
struct PairedEvents{
    long long destination;
    int countMinusMinus = 0;
    //int mmInversed = 0;
    int countMinusPlus = 0;
    //int mpInversed = 0;
    int countPlusMinus = 0;
    //int pmInversed = 0;
    int countPlusPlus = 0;
    //int ppInversed = 0;
    PairedEvents(long long dest):destination(dest){}
};

struct Cluster{

    long long index = 0;
    int length = 0;
    //int d = 0;
    int firstFragment = 0;
    int lastFragment = 0;
    long long firstPosition = 0;
    long long lastPosition = 0;
    int flag = -1;
    short sign = 0; // (left:-1  or  right:+1 or both:0) of a breakpoint
    Cluster(){}
    Cluster(long long ind,int len,int ff, int lf, long long fp, long long lp, int f, short s):
        index(ind),length(len),firstFragment(ff),lastFragment(lf),firstPosition(fp),lastPosition(lp),flag(f),sign(s){}
};
struct ClusterAndConnection{
    Cluster cluster;
    long long connectedEvent=-1;
    short connectedSign = 0;
    int pairDistance = 0;
    //int connectedFlag = 0;
    ClusterAndConnection(){}
    ClusterAndConnection(Cluster clus, long long conEvent, short conSign):cluster(clus),connectedEvent(conEvent),connectedSign(conSign){}
    ClusterAndConnection(Cluster clus, long long conEvent, short conSign, int pairDist):cluster(clus),connectedEvent(conEvent),connectedSign(conSign),pairDistance(pairDist){}
    //ClusterAndConnection(Cluster clus, long long conEvent, short conSign, int conFlag):cluster(clus),connectedEvent(conEvent),connectedSign(conSign),connectedFlag(conFlag){}
};

struct cCluster{
    Cluster cluster;
    int fragNum;
    long long corespEvent;
    cCluster(){}
    cCluster(Cluster c, int f, long long e):cluster(c),fragNum(f),corespEvent(e){}
};

class Signalling
{
public:
    Signalling(iRead* reads , string readName, string genomeName, int chunkSize, int numberOfThread, int numGap, string outputDir,long long numberOfReads1);
    ~Signalling();

    long long getGenomeLength();
    void changeBinGenome(long long , long long );

    void automation();
    void buildEventsBinGenome();
    void buildEvents();
    void buildPairedEvents();


    void printEventIntervals();
    void printPaiedEvents();


    void typeMatching();

    void buildSignal();

    long long correspondingEvent(long long);



    iRead* reads;

    vector<bool> eventBinGenome;
    //vector<bool> bGLeft;
    vector<iRead> unMappedReads;
    vector<iRead> informativeReads;
    //vector<PairedEvents> pairedEvents;

    vector<EventInterval> eventIntervals;
    //vector<vector<Cluster>> eventInformativeLMers;
    vector<vector<ClusterAndConnection>> eventInformativeClusters;
    vector<vector<PairedEvents>> eventPairedEvents;
    vector<int> eventSingleCountsMinus;
    vector<int> eventSingleCountsPlus;
    vector<vector<vector<int>>> readConnectedEvents;
    vector<vector<int>> eventConnectedEvents;

    int *mainSignal;

private:
    string readName;
    string genomeName;
    long long genomeLength;
    int chunkSize;
    int numberOfThread;
    int numGap;
    string outputDir;
    long long numberOfReads1;

};

#endif // SIGNALLING_H


