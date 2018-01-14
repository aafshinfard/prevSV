#ifndef PATH_H
#define PATH_H

#include <string>
#include <vector>

using namespace std;

class Path {

public:
    Path();
    Path(int input1,long long input2,int Rows, int flag, long long partLength, int d, string ref, int readLength);
    virtual ~Path();
    void addScore(long long newScore);
    long long getBirth();
    int getBirthRow();
    int getScore();
    void setTable(int RowNumber, long long value);
    vector <long long> getTable();
    bool isTruePath();
    void setTruePath();
    void setIndexTable(int RowNumber, long long value);
    vector <long long> getIndices();
    int getFlag();
    long long getMain();
    vector <string> getCigars();
    void addToCigars(std::string cigar, int number);
    bool clearance();
    bool hasMore(int limit);
    string getRef();
    long long birthIndex, mainIndex;
    long long tinyPosition;
    long long LastNonZeroIndex();//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int getLastRow();//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    long long getLastIndex();//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    void joining();//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    bool getSpliced();//++++++++++++++++++++++++++++++++++++++++++++++++++++
    string getTotalCigar();//++++++++++++++++++++++++++++++++++++++++++++++++++++
    void setTotalCigar(string Cigar);//++++++++++++++++++++++++++++++++++++++++++++++++++++

private:
    bool truePath;
    int Rows;
    int directionFlag;
    int static numberOfPaths;
    bool spliced;//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    int birthRow;
    int score;
    string ref;
    vector <long long> table;
    vector <long long> indexTable;
    vector <string> cigars;
    string totalCigar;//*************************************************
};

#endif // PATH_H
