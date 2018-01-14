/*
 * Path.cpp
 *
 *  Created on: Mar 5, 2015
 *      Author: amin
 */

#include "Path.h"
#include <stdlib.h>
#include <string>

using namespace std;

int Path::numberOfPaths = 0;

Path::Path() {

}

Path::Path(int input1, long long input2, int Rows, int flag, long long partLength, int d, string ref, int readLength) {
    truePath=false;
    spliced = false;
    this->Rows = Rows;
    this->ref = ref;
    score = 0;
    this->directionFlag= flag;
    // Table contains the flags of path nodes that by default are -1 and
    // when turned on are
    for (int k=0;k<Rows;k++){
        table.push_back(-1);
    }
    for (int m=0;m<Rows;m++){
        cigars.push_back("");
    }
    for (int l=0;l<Rows;l++){
        indexTable.push_back(0);
    }
    this->birthRow = input1;
    this->birthIndex = input2;
    this->totalCigar = "";//*************************************************
    if(flag==0 || flag==256){
        this->mainIndex = this->birthIndex - this->birthRow*(partLength-d);
        tinyPosition = this->mainIndex;
    }
    else if(flag==16 || flag==272){
        this->mainIndex = this->birthIndex + this->birthRow*(partLength-d)+partLength;
        tinyPosition = this->mainIndex - readLength;
    }

    numberOfPaths++;

}

Path::~Path() {
    table.clear();
    cigars.clear();
    indexTable.clear();
}

void Path::addScore(long long newScore) {
    this->score = (int)newScore + this->score;
    return;
}

// First index that creates this path
long long Path::getBirth() {
    return birthIndex;
}


//First row that creates this path
int Path::getBirthRow() {
    return birthRow;
}

int Path::getScore() {
    return this->score;
}

void Path::setTable(int RowNumber, long long value) {
    table[RowNumber] = value;
    return;
}

void Path::setIndexTable(int RowNumber, long long value) {
    indexTable[RowNumber] = value;
}

vector<long long> Path::getTable() {
    return table;
}

bool Path::hasMore(int limit) {
    bool flag=false;
    int counter=0;
    for(int i=0;i<indexTable.size();i++){
        if(indexTable[i]!=0)
            counter++;
        if(counter>limit){
            flag = true;
            break;
        }
    }
    return flag;
}

bool Path::isTruePath() {
    return truePath;
}

void Path::setTruePath() {
    truePath=true;
    return;
}

vector<long long> Path::getIndices() {
    return indexTable;
}

int Path::getFlag() {
    return this->directionFlag;
}

long long Path::getMain() {
    return this->mainIndex;
}

void Path::addToCigars(std::string cigar,int number) {
    cigars[number] = cigar;
    return;
}

vector<string> Path::getCigars() {
    return this->cigars;
}

bool Path::clearance() {
    if(this->tinyPosition>0){
        return true;
    }
    return false;
}

string Path::getRef() {
    return this->ref;
}

//***********************************************************
long long Path::LastNonZeroIndex() {
    long long j=0;
    for(int i=0;i<indexTable.size();i++){
        if (indexTable[i]>0){
            j=indexTable[i];
        }
    }
    return j;
}

int Path::getLastRow() {
    int lastrow = 0;
    for(int i=0;i<indexTable.size();i++){
        if(indexTable[i]!=0)
            lastrow = i;
    }
    return lastrow;
}

long long Path::getLastIndex() {
    int lastrow = 0;
    for(int i=0;i<indexTable.size();i++){
        if(indexTable[i]!=0)
            lastrow = i;
    }
    return indexTable[lastrow];
}

void Path::joining() {
    numberOfPaths--;
    spliced = true;
}

bool Path::getSpliced() {
    return spliced;
}


string Path::getTotalCigar(){
    return totalCigar;
}

void Path::setTotalCigar(string Cigar){
    totalCigar = Cigar;
}




