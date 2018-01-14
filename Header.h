#ifndef HEADER_H
#define HEADER_H
// ================================|
// ================================|
//> Settings:                      |
// ================================|

#define MAINADDRESS "/home/ameer/SVAnchoring/"
//#define MAINADDRESS "/home/ubuntu/SVAnchoring/"
#define BOWTIE "bowtie"
#define BOWTIE2ADDRESS "bowtie2"
#define BOWTIEBUILD "bowtie-build"
//bool multiReadFiles = false;

// ================================|
//> End of Settings:               |
// ================================|

#include "Tools.h"



using namespace std;
struct FragLoci{
    int Fragment = 0;
    int Flag = -1;
    long long Position = 0;

    FragLoci():Fragment(0),Flag(-1),Position(0){}
    FragLoci(int a, int b, long long c) : Fragment(a), Flag(b), Position(c) {}

    bool operator < (const FragLoci& str) const
    {
        return (Fragment < str.Fragment);
    }
    bool operator < (const int& ineteger) const
    {
        return (Fragment < ineteger);
    }
    bool operator > (const int& ineteger) const
    {
        return (Fragment > ineteger);
    }
    int operator - (const FragLoci& str) const
    {
        return (Fragment - str.Fragment);
    }
};

struct iRead {
    long long index = 0;
    int length = 0;
    int d = 0;

    int firstFragment = 0;
    int lastFragment = 0;
    long long firstPosition = 0;
    long long lastPosition = 0;
    //-/ Cluster2 details:
    vector<int> firstFragClus2 = {0};
    vector<int> lastFragClus2 = {0};
    vector<long long> firstPosClus2 = {0};
    vector<long long> lastPosClus2 = {0};
    vector<int> flagClus2 = {-1};

    int clusterCounts = 0;
    //-/ end of Cluster2 details

    int flag = -1; // 0-16 | 1-17 (E-I) | 2-18 (E-E) | 3-19 (Unsample)
                   // 4-20 (unique-unspliced)
                   // 5-21 (unique-spliced)
                   // 6-22 (multi-unspliced)
                   // 7-23 (multi-spliced)
                   // +100 (small-exon)

    string cigar = "~";

#ifdef runInRealMode
#else
    long long* origin = new long long[5];
#endif
};

struct iReadNext {
    long long index = 0;
    int length = 0;
    int d = 0;

    int firstFragment = 0;
    int lastFragment = 0;
    long long firstPosition = 0;
    long long lastPosition = 0;
    //-/ Cluster2:
    vector<int> firstFragClus2 = {0};
    vector<int> lastFragClus2 = {0};
    vector<long long> firstPosClus2 = {0};
    vector<long long> lastPosClus2 = {0};
    vector<int> flagClus2 = {-1};

    int clusterCounts = 0;
    //-/ end of Clus2
    int flag = -1;
    string cigar = "~";
    long long previousIndex;
};

struct onInterval {
    long long index;
    int fragment;
    int distance;
};

struct interval {
    long long start = 0;
    long long end = 0;

    vector<onInterval> rightRead;
    vector<onInterval> leftRead;
    vector<string> rightSeq;
    vector<string> leftSeq;
};

struct intervalVirtualBowtie2 {
    long long start = 0;
    long long end = 0;

    vector<long long> inReadNumber;
};

struct fragInfromation {
    int answer1;
    int answer2;
    int len1;
    int len2;
    int score;
    bool isX;
    string cigar;
};

struct readSplicedInformation {
    long long positionInGenome;
    int positionInRead;
    string cigar;
};

struct splicedFrag {
    int firstFragment = 0;
    int lastFragment = 0;
    long long firstPosition = 0;
    long long lastPosition = 0;
};

struct splicedRead {
    long long index = 0;
    string cigar = "";
    vector<splicedFrag> fragments;
};

struct gene {
    long long start = 0;
    long long end = 0;
    bool isForward;

    vector<interval> inExon;
    vector<splicedRead> inRead;
};

#endif // HEADER_H
