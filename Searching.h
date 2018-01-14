#ifndef SEARCHING_H
#define SEARCHING_H
#include "Aligning.h"
#include "Header.h"
#include "Tools.h"

class Searching {
private:
    long long numberOfReads;

    iRead* reads;
    vector<interval>* depth;
    int characterPerRow;
    string genomeFileName;
    int chunkSize;
    int indel;
    int mismatch;
    int matchScore;
    int indelScore;
    int mismatchScore;

    string getGenome(long long, int);
    string complement(string);
    string changeCigar(string);
    tuple<int, int> findPosition(int, bool, int, int);
    tuple<bool,long long, long long, long long> calculateIntronOrExonAndDistance(long long);
    string globalAlign(string, string);
    int score(char, char);
    string getTheChunkRead(long long, int);

    vector<readSplicedInformation> amirHossein4(string, vector<string>, int, int, int, bool, long long);
    vector<fragInfromation> amirHossein5(string, string, int, vector <int>);

    vector<int> localAlignerForSearch(string, string, vector <int>);
    double similarity_score(char, char, int);
    double find_array_max(double array[], int length, int* ind);


public:
    Searching(iRead*, vector<interval>*, long long, int, string, int, int, int, int, int, int);
    ~Searching();

    void searchingWithNoCutting(int);
    void searchingWithCutting3(int);
    void smallExon(string, long long, int, bool);
    void prepareForSearchingWithCutting2();
    void breakReadsForThreads(int);
    void concatExons();
};

#endif // SEARCHING_H
