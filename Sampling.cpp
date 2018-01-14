//
//  Sampling.cpp
//  SV-Detection
//
//  Created by Farid Rashidi on 10/21/15.
//  Copyright © 2015 Farid Rashidi. All rights reserved.
//

#include "Sampling.h"

Sampling::Sampling(iRead* reads, vector<interval>* depth, int characterPerRow, string genomeFileName, int floatingEdge, string outputDir) {
    this->reads = reads;
    this->depth = depth;
    this->genomeFileName = genomeFileName;
    this->characterPerRow = characterPerRow;
    this->floatingEdge = floatingEdge;
    this->outputDir = outputDir;
}


Sampling::~Sampling() {

}
inline string Sampling::revComplemACGT(string readSegment){
    for (char& c : readSegment) {
        switch(c) {
        case 'A':
            c = 'T';
            break;
        case 'G':
            c = 'C';
            break;
        case 'C':
            c = 'G';
            break;
        case 'T':
            c = 'A';
            break;
        case 'a':
            c = 't';
            break;
        case 'g':
            c = 'c';
            break;
        case 'c':
            c = 'g';
            break;
        case 't':
            c = 'a';
            break;
        }
    }
    readSegment = string ( readSegment.rbegin(), readSegment.rend() );
    return readSegment;
}

long long Sampling::prepareSampling(string readName, bool pairedEndedReads = false, int readsFilesCount = 1  ) {
    long long counter = 0;
    if(!pairedEndedReads){
        if(readsFilesCount==1){
            string name = readName;
            ifstream ifstr(name.c_str());
            name = outputDir+"an_reads_d=0.fq";
            ofstream ofstr2(name.c_str());

            string line;
            string privous_line;
            string segment;

            while (getline(ifstr, line)) {
                if (counter % 4 == 0){
//                    stringstream originInput(line.substr(2));
#ifdef runInRealMode
#else
                    int p = 0;
//                    while(getline(originInput, segment, '-')) {
//                        this->reads[counter/4].origin[p] = stoll(segment.c_str());
//                        p++;
//                    }
#endif
                    replace(line.begin(), line.end(), ' ', '-');
                    privous_line = line;
                } else if (counter % 4 == 1) {
                    this->reads[(counter-1)/4].length = (int) line.size();
                    this->reads[(counter-1)/4].index = (counter-1)/4+1;
                    ofstr2 << privous_line << "_" << (counter-1)/4+1 << "\n";
                    ofstr2 << line << "\n";
                } else if (counter % 4 == 3) {
                    ofstr2 << "+\n" << line << "\n";
                }
                ++counter;
            }
            ifstr.close();
            ofstr2.close();

        }else{
            long long counter = 0;
            for(int i = 0 ; i < readsFilesCount ; i++){
                string name = readName+"_"+convertNumToStr(i)+".fq";
                ifstream ifstr(name.c_str());
                name = outputDir+"an_reads_d=0_f"+convertNumToStr(i)+".fq";
                ofstream ofstr2(name.c_str());

                string line;
                string privous_line;
                string segment;

                while (getline(ifstr, line)) {
                    if (counter % 4 == 0){
//                        stringstream originInput(line.substr(2));
    #ifdef runInRealMode
    #else
//                        int p = 0;
//                        while(getline(originInput, segment, '-')) {
//                            this->reads[counter/4].origin[p] = stoll(segment.c_str());
//                            p++;
//                        }
    #endif
                        replace(line.begin(), line.end(), ' ', '-');
                        privous_line = line;
                    } else if (counter % 4 == 1) {
                        this->reads[(counter-1)/4].length = (int) line.size();
                        this->reads[(counter-1)/4].index = (counter-1)/4+1;
                        ofstr2 << privous_line << "_" << (counter-1)/4+1 << "\n";
                        ofstr2 << line << "\n";
                    } else if (counter % 4 == 3) {
                        ofstr2 << "+\n" << line << "\n";
                    }
                    ++counter;
                }
                ofstr2.close();
                ifstr.close();
            }


            return (counter-1)/4+1;
        }
    }else if(pairedEndedReads){
        string name1 = readName+"_1.fq";
        string name2 = readName+"_2.fq";
        ifstream ifstr1(name1.c_str());
        ifstream ifstr2(name2.c_str());
        string name = outputDir+"an_reads_d=0.fq";
        ofstream ofstr2(name.c_str());

        string linef1,linef2;
        string privous_line;
        string segment;

        while (getline(ifstr1, linef1) && getline(ifstr2, linef2)) {
            if (counter % 4 == 0){
                stringstream originInput(linef1.substr(2));
#ifdef runInRealMode
#else
                int p = 0;
                while(getline(originInput, segment, '-')) {
                    this->reads[counter/4].origin[p] = stoll(segment.c_str());
                    p++;
                }
#endif
                replace(linef1.begin(), linef1.end(), ' ', '-');
                privous_line = linef1;
            } else if (counter % 4 == 1) {
                this->reads[(counter-1)/4].length = ((int) linef1.size()) * 2;
                this->reads[(counter-1)/4].index = (counter-1)/4+1;
                ofstr2 << privous_line << "_" << (counter-1)/4+1 << "\n";
                ofstr2 << linef1 << revComplemACGT(linef2) << "\n";
            } else if (counter % 4 == 3) {
                ofstr2 << "+\n" << linef1 <<string(linef2.rbegin(),linef2.rend())<< "\n";

            }
            ++counter;
        }
        ifstr1.close();
        ofstr2.close();

    }
    return (counter-1)/4+1;
}


void Sampling::prepareSampling(string readName, long long numberOfRandom, long long numberOfReads) {
    srand(unsigned(time(0)));

    string name = readName;
    ifstream ifstr(name.c_str());
    name = outputDir+"ans_reads_d=0.fq";
    ofstream ofstr(name.c_str());

    string line;
    string privous_line;
    string segment;
    long long counter = 0;
    while (getline(ifstr, line)) {
        if (counter % 4 == 0){
            stringstream originInput(line.substr(2));
#ifdef runInRealMode
#else
            int p = 0;
            while(getline(originInput, segment, '-')) {
                this->reads[counter/4].origin[p] = stoll(segment.c_str());
                p++;
            }
#endif
            replace(line.begin(), line.end(), ' ', '-');
            privous_line = line;
        } else if (counter % 4 == 1) {
            this->reads[(counter-1)/4].length = (int) line.size();
            this->reads[(counter-1)/4].index = (counter-1)/4+1;
            ofstr << privous_line << "_" << (counter-1)/4+1 << "\n";
            ofstr << line << "\n";
        } else if (counter % 4 == 3) {
             ofstr << "+\n" << line << "\n";
        }
        ++counter;
    }
    ifstr.close();
    ofstr.close();
    return;
}


void Sampling::prepareExonSequence() {
    string name = outputDir+"EXON.fa";
    ofstream ofstr(name.c_str());
    ofstr << ">EXON\n";
    string exonSequence;
    int remain = 0;
    for (long long i = 0; i < this->depth->size(); i++) {
        exonSequence = getGenome((this->depth->at(i).start-this->floatingEdge)-1, this->depth->at(i).end+this->floatingEdge - (this->depth->at(i).start-this->floatingEdge) + 1);



        if (remain != 0) {
            if (exonSequence.substr(0,remain).size() < remain) {
                ofstr << exonSequence;
                remain = remain - (int) exonSequence.substr(0,remain).size();
            } else {
                ofstr << exonSequence.substr(0,remain) << "\n";
                exonSequence = exonSequence.substr(remain, exonSequence.size()-remain);
                if (exonSequence.size() != 0) {
                    for (int j = 0; j <= ceil(exonSequence.size()/this->characterPerRow); j++) {
                        if (exonSequence.substr(j*this->characterPerRow, this->characterPerRow).size() == this->characterPerRow) {
                            ofstr << exonSequence.substr(j*this->characterPerRow, this->characterPerRow) << "\n";
                        } else {
                            ofstr << exonSequence.substr(j*this->characterPerRow, this->characterPerRow);
                            remain = this->characterPerRow - (int)exonSequence.substr(j*this->characterPerRow, this->characterPerRow).size();
                        }
                    }
                } else {
                    remain = 0;
                }
            }
        } else {
            for (int j = 0; j <= ceil(exonSequence.size()/this->characterPerRow); j++) {
                if (exonSequence.substr(j*this->characterPerRow, this->characterPerRow).size() == this->characterPerRow) {
                    ofstr << exonSequence.substr(j*this->characterPerRow, this->characterPerRow) << "\n";
                } else {
                    ofstr << exonSequence.substr(j*this->characterPerRow, this->characterPerRow);
                    remain = this->characterPerRow - (int)exonSequence.substr(j*this->characterPerRow, this->characterPerRow).size();
                }
            }
        }


    }
    if (remain != 0) {
        ofstr << "\n";
    }
    ofstr.close();
    return;
}


void Sampling::createIndexBowtie() {
    string command = string(BOWTIEBUILD) + " --quiet "+outputDir+"EXON.fa EXON";
    system(command.c_str());
    return;
}


string Sampling::getGenome(long long first, long long length) {
    struct stat sb;
    off_t len;
    char *p;
    int fd;
    string name = this->genomeFileName;
    ifstream ifstr(name.c_str());
    string genome;
    string firstLine;
    getline(ifstr, firstLine);
    int firstLineCharCount = (int)firstLine.size();
    ifstr.close();
    first += firstLineCharCount+1 + first/this->characterPerRow;
    fd = open(name.c_str(), O_RDONLY);
    if (fd == -1)
        perror ("open");

    if (fstat (fd, &sb) == -1)
        perror ("fstat");

    p = (char *)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

    if (close (fd) == -1)
        perror ("close");

    string temp = "";
    for (len = first; temp.size() < length; ){
        if (p[len] != '\n')
            temp += p[len];

        ++len;
    }
    if (munmap (p, sb.st_size) == -1)
        perror ("munmap");

    close(fd);
    return temp;
}
