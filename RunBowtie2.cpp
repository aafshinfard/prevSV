//
//  RunBowtie2.cpp
//  RNA-Seq
//
//  Created by Farid Rashidi on 7/31/16.
//  Copyright © 2016 Farid Rashidi. All rights reserved.
//

#include "RunBowtie2.h"

RunBowtie2::RunBowtie2(iRead* reads, vector<interval>* depth, string outputDir) {
    this->reads = reads;
    this->depth = depth;
    this->outputDir = outputDir;
}


RunBowtie2::~RunBowtie2() {

}


void RunBowtie2::run(int d, string name, int numberOfThread, int mmseed, int lseed, string indexName) {
    string file = outputDir+name+"_reads_d="+convertNumToStr(d)+".fq";
    string file2 = outputDir+name+"_reads_d="+convertNumToStr(d)+".sam";
    string command = string(BOWTIE2ADDRESS) + " -t --threads "+convertNumToStr(numberOfThread)+" --local --reorder -N "+convertNumToStr(mmseed)+" -L "+convertNumToStr(lseed)+" -x "+indexName+" -U "+file+" -S "+file2+" 2>> "+outputDir+"log-bowtie.txt";
    system(command.c_str());

    vector<readNumAndSnap> alignsInThisStep;
    alignsInThisStep = this->rearrange(file2);
    for (long long i = 0; i < alignsInThisStep.size(); i++) {
        long long whereAligned = alignsInThisStep.at(i).snap;
        this->addIntervalToVirtualDepth(whereAligned, whereAligned+150);
    }

    this->mergeDepthVirtualWithDepth(30);

    return;
}


void RunBowtie2::mergeDepthVirtualWithDepth(int threshold) {

    long long number = 0;
    long long i = 0;

    for (vector<intervalVirtualBowtie2>::iterator it=this->iDepth_virtual.begin(); it!=this->iDepth_virtual.end(); ) {

        for (long long j = number; j < this->depth->size(); j++) {

            long long a = this->depth->at(j).start - threshold;
            long long b = this->depth->at(j).end + threshold;

            if (a < it->start && it->start < b) {
                if (j == 0) {
                    number = 0;
                } else {
                    number = j-1;
                }
                this->mergeOnIntervalOfDepthVirtualwithDepth(it->start, it->end, number);
                it = this->iDepth_virtual.erase(it);
                if (i == 0) {
                    it--;
                    i--;
                } else {
                    it -= 2;
                    i -= 2;
                }
                break;
            } else if (a < it->end && it->end < b) {
                if (j == 0) {
                    number = 0;
                } else {
                    number = j-1;
                }
                this->mergeOnIntervalOfDepthVirtualwithDepth(it->start, it->end, number);
                it = this->iDepth_virtual.erase(it);
                if (i == 0) {
                    it--;
                    i--;
                } else {
                    it -= 2;
                    i -= 2;
                }
                break;
            }
        }

        ++it;
        ++i;
    }

    return;
}


void RunBowtie2::mergeOnIntervalOfDepthVirtualwithDepth(long long a, long long b, long long index) {

    long long threshold = 5;

    for (long long i = index; i < this->depth->size(); i++) {
        if (this->depth->at(i).start <= a && this->depth->at(i).end >= b) {
            return;
        }
        if (this->depth->at(i).start > b + threshold || this->depth->at(i).end < a - threshold) {
            continue;
        } else {
            long long tempA = this->depth->at(i).start;
            long long tempB = this->depth->at(i).end;
            this->depth->erase(this->depth->begin()+i);
            this->mergeOnIntervalOfDepthVirtualwithDepth(min(tempA, a), max(tempB, b), i-1);
            return;
        }
    }
    long long insert_index = index;
    while (insert_index < this->depth->size() && this->depth->at(insert_index).start < b) {
        insert_index++;
    }
    this->depth->insert(this->depth->begin()+insert_index, interval());
    this->depth->at(insert_index).start = a;
    this->depth->at(insert_index).end = b;
    return;
}


void RunBowtie2::addIntervalToVirtualDepth(long long a, long long b) {

    if (this->iDepth_virtual.size() == 0) {
        this->iDepth_virtual.push_back(intervalVirtualBowtie2());
        this->iDepth_virtual.at(0).start = a;
        this->iDepth_virtual.at(0).end = b;
        return;
    }
    for (long long i = 0; i < this->iDepth_virtual.size(); i++) {
        if (this->iDepth_virtual.at(i).start <= a && this->iDepth_virtual.at(i).end >= b) {
            return;
        }
        if (this->iDepth_virtual.at(i).start > b || this->iDepth_virtual.at(i).end < a) {
            continue;
        } else {
            long long tempA = this->iDepth_virtual.at(i).start;
            long long tempB = this->iDepth_virtual.at(i).end;
            this->iDepth_virtual.erase(this->iDepth_virtual.begin()+i);
            this->addIntervalToVirtualDepth(min(tempA, a), max(tempB, b));
            return;
        }
    }
    long long insert_index = 0;
    while (insert_index < this->iDepth_virtual.size() && this->iDepth_virtual.at(insert_index).start < b) {
        insert_index++;
    }
    this->iDepth_virtual.insert(this->iDepth_virtual.begin()+insert_index, intervalVirtualBowtie2());
    this->iDepth_virtual.at(insert_index).start = a;
    this->iDepth_virtual.at(insert_index).end = b;
    return;

}


bool compareZippedSnapAndReadNumNew(readNumAndSnap r1, readNumAndSnap r2) {return r1.readNum < r2.readNum;}


vector<readNumAndSnap> RunBowtie2::rearrange(string filename) {
    vector<readNumAndSnap> vec;
    SamAnalyzer* analyze = new SamAnalyzer(filename);
    analyze->parseSAMFile();
    vec = analyze->getOutput();
    sort(vec.begin(), vec.end(), compareZippedSnapAndReadNumNew);
    delete analyze;
    return vec;
}


int RunBowtie2::isReverse(readNumAndSnap r) {
    int* flags = this->flagAnalyzer(r.flag);
    if (flags[4]) {
        delete [] flags;
        return 1;
    } else {
        delete [] flags;
        return 0;
    }
}


int* RunBowtie2::flagAnalyzer(int flag) {
    int * flags =  new int[8];
    for (int i=0;i<8;i++){flags[i] = 0;}
    for (int i = 7; i >= 0; i--) {

        if (flag >= (int)pow(2, i)) {
            flags[i] = 1;
            flag -= (int)pow(2, i) ;
        }
    }
    return flags;
}
