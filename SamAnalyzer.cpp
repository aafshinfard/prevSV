//
//  SamAnalyzer.cpp
//  MetaAligner
//
//  Created by ahmad on 12/1/93.
//  Copyright (c) 1393 ahmad. All rights reserved.
//

#include "SamAnalyzer.h"


SamAnalyzer::SamAnalyzer(string fileName){
    this->samFile = fileName;
}


SamAnalyzer::~SamAnalyzer(){
    vector<readNumAndSnap>().swap(this->output);
}


vector<readNumAndSnap> SamAnalyzer::getOutput() {
    return this->output;
}


long long SamAnalyzer::resolveHeaderReadNum(string input) {
    stringstream test(input);
    string segment;
    getline(test, segment, '_');
    getline(test, segment, '_');
    return stoll(segment.c_str());
}


void SamAnalyzer::parseSAMFile() {
    string line;
    ifstream stream(this->samFile.c_str());
    while (getline(stream, line)) {
        if (line[0] == '@') continue;
        else {
            readNumAndSnap sam;
            istringstream st(line);
            string first;
            st>>first;
            first = string(first.begin()+1, first.end());
            sam.readNum = this->resolveHeaderReadNum(first);
            string second;
            st>>second;
            sam.flag = atoi(second.c_str());
            st>>second>>second;
            sam.snap = stoll(second.c_str());
            output.push_back(sam);
        }
    }
    return;
}
