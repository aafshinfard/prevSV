#ifndef SAMANALAYZER_H
#define SAMANALAYZER_H


#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

struct readNumAndSnap {
    long long readNum;
    long long snap;
    int flag;
};

class SamAnalyzer{
public:
    SamAnalyzer(string);
    ~SamAnalyzer();
    void parseSAMFile();
    vector<readNumAndSnap> getOutput();

private:
    string samFile;
    vector<readNumAndSnap> output;

    long long resolveHeaderReadNum(string);
};

#endif // SAMANALAYZER_H
