//  Anchoring.cpp
//  SV-Detection
//
//  Created by Ameerhosein Afshinfard - BRLab - Sharif Uiversity of echnology
//  Copyright (c) 2016 Ameerhosein Afshinfard. All rights reserved.
//

#include "Anchoring.h"
#include "vardefs.h"


Anchoring::Anchoring(iRead* reads, vector<interval>* depth, string referenceName, string outputName, bool isSample, int numGap,
                     int chunkSize, int d, int floatingEdge, int dashV, int numberOfThread, string outputDir) {
    this->reads = reads;
    this->depth = depth;
    this->correspond = 0;
    this->correspondCluster2 = 0;
    this->referenceName = referenceName;
    this->outputName = outputName;
    this->isSample = isSample;
    this->chunkSize = chunkSize;
    this->d = d;
    this->numGap = numGap;
    this->floatingEdge = floatingEdge;
    this->dashV = dashV;
    this->numberOfThread = numberOfThread;
    this->outputDir = outputDir;

}


Anchoring::~Anchoring() {
    vector< vector<readNumAndSnap> >().swap(this->notAligned);
}

long long Anchoring::getNumberFiles(){
    return this->numberFiles;
}

void Anchoring::prepareFiles(){
    // taghir header read ha: @headerGhabli_tool_shomare
    // peyda kardan bishtarin toole read ha

    //string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+".fq";
    string name = outputDir+this->outputName+"_reads_d=0.fq";
    if(this->d > 0 )
        cout<<"a";
    ifstream ifstr(name.c_str());
    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_my.fq";
    ofstream ofstr(name.c_str());
    string line;
    string privous_line;
    long long counter = 0;
    int readLenMax = 0;
    long long newreadcount=0;
    bool unanched = true;
    if( this->d == 0 ){
        while (getline(ifstr, line)) {

            if (counter % 4 == 0){
                privous_line = line;
                newreadcount++;
            } else if (counter % 4 == 1 ) {
                ofstr << privous_line << "\n";
                ofstr << line.substr(this->d) << "\n";
                if(readLenMax < line.substr(this->d).size()) readLenMax = (int) line.substr(this->d).size();
            } else if (counter % 4 == 3) {
                ofstr << "+\n" << line.substr(this->d) << "\n";
            }
            ++counter;
        }
    }else{
        while (getline(ifstr, line)) {
            string segmorde;
            stringstream test(line);
            getline(test, segmorde, '_');
            getline(test, segmorde, '_');

            if (counter % 4 == 0){
                unanched = (this->reads[stoll(segmorde.c_str())-1].flag == -1);
                if (unanched){
                    privous_line = line;
                    newreadcount++;
                }
            } else if (counter % 4 == 1 ) {
                 if (unanched){
                     ofstr << privous_line << "\n";
                     ofstr << line.substr(this->d) << "\n";
                     if(readLenMax < line.substr(this->d).size()) readLenMax = (int) line.substr(this->d).size();
                 }
            } else if (counter % 4 == 3) {
                if (unanched)
                    ofstr << "+\n" << line.substr(this->d) << "\n";
            }
            ++counter;
        }
    }
    cout<<"\n ++++++ reads in this Loop: "<<newreadcount<<endl;
    ifstr.close();
    ofstr.close();

    // tafkik read ha be chunk ha
    this->numberFiles = (int) (readLenMax+this->chunkSize-20)/this->chunkSize;
    int* myStepArray = this->getMyArraySteps(this->numberFiles);
    ofstream* ofs = new ofstream[this->numberFiles];
    for (int i = 0; i < this->numberFiles; i++) {
        name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(i+1)+".fq";
        ofs[i].open(name.c_str());
    }
    counter = 0;
    int steps = 0;
    int remain = 0;
    line = "";
    privous_line = "";
    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_my.fq";
    ifstream ifstr2(name.c_str());
    while (getline(ifstr2, line)) {
        if (counter % 4 == 0) {
            privous_line = line;
        } else if (counter % 4 == 1) {
            steps = (int) (line.size()+this->chunkSize-20)/this->chunkSize;
            remain = (int) line.size() - steps*this->chunkSize;
            for (int i = 0; i < steps; i++) {
                ofs[myStepArray[i]-1] << privous_line << "\n";
                if (i%2 == 1) {
                    if (i == 1) {
                        ofs[myStepArray[i]-1] << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
                    } else {
                        ofs[myStepArray[i]-1] << line.substr(line.size()-this->chunkSize*(i+1)/2-remain, this->chunkSize) << "\n";
                    }
                } else {
                    ofs[myStepArray[i]-1] << line.substr(this->chunkSize*i/2, this->chunkSize) << "\n";
                }
            }
        } else if (counter % 4 == 3) {
            for (int i = 0; i < steps; i++) {
                if (i%2 == 1) {
                    if (i == 1){
                        ofs[myStepArray[i]-1] << "+\n" << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
                    } else {
                        ofs[myStepArray[i]-1] << "+\n" << line.substr(line.size()-this->chunkSize*(i+1)/2-remain, this->chunkSize) << "\n";
                    }
                } else {
                    ofs[myStepArray[i]-1] << "+\n" << line.substr(this->chunkSize*i/2, this->chunkSize) << "\n";
                }
            }
        }
        ++counter;
    }
    for (int i = 0; i < this->numberFiles; i++) ofs[i].close();
    ifstr2.close();
    delete[] myStepArray;
    delete[] ofs;
    return;
}
void Anchoring::prepareFiles(int range1, int range2, int readsFileCount,int slide=0, int iteration=0){ // range1 to range2 and end-range1 to end-range2
    // taghir header read ha: @headerGhabli_tool_shomare
    // peyda kardan bishtarin toole read ha
     long long newreadcount=0;
    if(range1==1){
        cout<<"\n Finding Longest Read...\n";
        int readLenMax = 0;
        for(int i=0 ; i < readsFileCount ; i++){
            //string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(i)+".fq";
            string name = outputDir+this->outputName+"_reads_d=0_f"+convertNumToStr(i)+".fq";
            ifstream ifstr(name.c_str());
            name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_f"+convertNumToStr(i)+"_my.fq";
            ofstream ofstr(name.c_str());
            string line = "";
            string privous_line = "";
            long long counter = 0;
            bool unanched = true;

            if(this->d == 0 ){
                while (getline(ifstr, line)) {
                    if (counter % 4 == 0){
                        privous_line = line;
                    } else if (counter % 4 == 1) {
                        ofstr << privous_line << "\n";
                        ofstr << line.substr(this->d) << "\n";
                        if(readLenMax < line.substr(this->d).size())
                            readLenMax = (int) line.substr(this->d).size();
                    } else if (counter % 4 == 3) {
                        ofstr << "+\n" << line.substr(this->d) << "\n";
                    }
                    ++counter;
                }
            }else{
                while (getline(ifstr, line)) {
                    string segmorde;
                    stringstream test(line);
                    getline(test, segmorde, '_');
                    getline(test, segmorde, '_');

                    if (counter % 4 == 0){
                        unanched = (this->reads[stoll(segmorde.c_str())-1].flag == -1);
                        if (unanched){
                            privous_line = line;
                            newreadcount++;
                        }
                    } else if (counter % 4 == 1) {
                        if (unanched){
                            ofstr << privous_line << "\n";
                            ofstr << line.substr(this->d) << "\n";
                            if(readLenMax < line.substr(this->d).size())
                                readLenMax = (int) line.substr(this->d).size();
                        }
                    } else if (counter % 4 == 3) {
                        if (unanched)
                            ofstr << "+\n" << line.substr(this->d) << "\n";
                    }
                    ++counter;
                }
            }
            ifstr.close();
        }
        cout<<"\n ++++++ reads in this Loop: "<<newreadcount<<endl;
        if(this->d > 0 )
            cout<<"a";
        this->numberFiles = (int) (readLenMax+this->chunkSize-20)/this->chunkSize;
        if(this->numberFiles % 2 == 1)
            this->numberFiles++;
        // //////////////////////////////////////
        // ///////**********************/////////
        // //////////|  adding +1   |////////////
        // //////////|  to odd      |////////////
        // //////////| numberFiles  |////////////
        // ///////**********************/////////
        // //////////////////////////////////////
        cout<<"\n Longest Read Found Successfully...\n";
    }
    if( range2 > (this->numberFiles/2) ){
        range2 = (this->numberFiles/2);
    }
    int end_ofs = (range2-range1+1)*2 - 1;
    ofstream* ofs = new ofstream[end_ofs + 1];

    // system("rm "+outputDir+"m_reads_d=*_*.fq");
    string name;
    for (int i = range1; i <= range2; i++) {
        name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(i)+".fq";
        ofs[i-range1].open(name.c_str());
    }
    for (int i = range1; i <= range2; i++) {
        name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(i+(range2-range1+1))+".fq";
        ofs[i-range1+(range2-range1+1)].open(name.c_str());
    }
    //int* myStepArray = this->getMyArraySteps(this->numberFiles);
    for(int k=0 ; k < readsFileCount ; k++){
        //string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+to_string(k)+".fq";
        string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_f"+convertNumToStr(k)+"_my.fq";
        ifstream ifstr(name.c_str());
        string line = "";
        string privous_line = "";
        long long counter = 0;
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        int steps;
        int remain;
        while (getline(ifstr, line)) {
            steps = (int) (line.size()+this->chunkSize-20)/this->chunkSize;
            remain = (int) line.size() - steps*this->chunkSize;
            if (counter % 4 == 0) {
                privous_line = line;
            }
            if (counter % 4 == 1){
                if( range1 == 1 ){
                    ofs[ 0 ] << privous_line << "\n";
                    ofs[ 0 ] << line.substr(0, this->chunkSize) << "\n";
                    ofs[ end_ofs ] << privous_line << "\n";
                    ofs[ end_ofs ] << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
                }
                int both = min( range2, steps/2 ) + 1;
                int i;
                for (i = range1 + (range1==1 ? 1 : 0); i < both ; i++){

                    ofs[ i-range1 ] << privous_line << "\n";
                    ofs[ i-range1 ] << line.substr(this->chunkSize*(i-1), this->chunkSize) << "\n";

                    ofs[ end_ofs - (i-range1) ] << privous_line << "\n";
                    ofs[ end_ofs - (i-range1) ] << line.substr(line.size()-this->chunkSize*i -remain, this->chunkSize) << "\n";
                }
                if( i < range2 )
                    if( steps %  2 == 1 ){
                        ofs[ i-range1 ] << privous_line << "\n";
                        ofs[ i-range1 ] << line.substr(this->chunkSize*(i-1), this->chunkSize) << "\n";
                    }
            }
            if (counter % 4 == 3){
                if( range1 == 1 ){
                    ofs[ 0 ] << "+\n" << line.substr(0, this->chunkSize) << "\n";
                    ofs[ end_ofs ] << "+\n" << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
                }
                int both = min( range2, steps/2 ) + 1;
                int i;
                for (i = range1 + (range1==1 ? 1 : 0); i < both ; i++){
                    ofs[ i-range1 ] << "+\n" << line.substr(this->chunkSize*(i-1), this->chunkSize) << "\n";
                    ofs[ end_ofs - (i-range1) ] << "+\n" << line.substr(line.size()-this->chunkSize*i -remain, this->chunkSize) << "\n";
                }
                if( i < range2 )
                    if( steps %  2 == 1 ){
                        ofs[ i-range1 ] << "+\n" << line.substr(this->chunkSize*(i-1), this->chunkSize) << "\n";
                    }
//                +Meta-Aligner:
//                for (int i = range1; i < range2; i++)
//                    if (slide*iteration+chunkSize*(i+1) <= line.size() / 2){
//                        ofs[i-range1]<<"+\n"<<line.substr(slide*iteration+chunkSize*i, chunkSize)<<endl;
//                        ofs[i-range1+(range2-range1)]<<"+\n"<<line.substr(line.size()-slide*iteration+chunkSize*i-chunkSize, chunkSize)<<endl;
//                    }else{
//                        ofs[i-range1]<<"+\nX"<<endl;
//                        ofs[i-range1+(range2-range1)]<<"+\nX"<<endl;
//                    }
            }
            ++counter;
//                for (int i = range1; i < range2 && i < (steps+1)/2; i++){
//                    ofs[i-range1]<<prevline<<"\n";
//                    ofs[i-range1+(range2-range1)]<<prevline<<"\n";
//                    if (slide*iteration+chunkSize*(i+1) <= line.size() / 2 ){
//                        ofs[i-range1]<<line.substr(slide*iteration+chunkSize*i, chunkSize)<<"\n";
//                        ofs[i-range1+(range2-range1)]<<line.substr(line.size()-slide*iteration+chunkSize*i-chunkSize, chunkSize)<<"\n";
//                    }else{
//                        ofs[i-range1]<<"X"<<endl;
//                        ofs[i-range1+(range2-range1)]<<"X"<<endl;
//                    }
//                }
        }
        ifstr.close();
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
        // +++++++++++++++++++++++++++++
    }


    for (int i = range1; i < range2; i++){
        ofs[i-range1].close();
        ofs[i-range1+(range2-range1)].close();
    }
    delete[] ofs;
    return;
}

void Anchoring::createNewReads(int d) {
    string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+".fq";
    ifstream ifstr(name.c_str());
    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(d)+".fq";
    ofstream ofstr(name.c_str());
    string line;
    string test;
    string segment;
    long long lineCounter = 0;
    while (getline(ifstr, line)) {
        if (lineCounter % 4 == 0) {
            stringstream test(line);
            getline(test, segment, '_');
            getline(test, segment, '_');
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<line<<"\n";
            }
        } else if (lineCounter % 4 == 1) {
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<line<<"\n";
            }
        } else if (lineCounter % 4 == 3) {
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<"+\n"<<line<<"\n";
            }
        }
        ++lineCounter;
    }
    ifstr.close();
    ofstr.close();
    return;
}


void Anchoring::writingNotAligned() {
    string name = outputDir+"all_notAlign.txt";
    ofstream ofstr(name.c_str());
    for (long long i = 0; i < this->notAligned.size(); i++) {
        ofstr << this->notAligned.at(i).at(0).readNum << "\t";
        for (int j = 0; j < this->notAligned.at(i).size(); j++) {
            ofstr << this->notAligned.at(i).at(j).snap << "\t";
        }
        ofstr << endl;
    }
    ofstr.close();
    return;
}

void Anchoring::nextStepAlignment() {
    //-/? numberFiles : number of files for different chunks? for example max readlength=150 l=30 => 5
    //    int** myArrayAlignmentCorrespond = new int*[this->numberFiles];
    // .... deleted
    //cout<<"here 1 ";
    int* myArraySteps = this->getMyArraySteps(this->numberFiles); //-/ order of mappings
    //int* myFragOrders = this->getFragsOrders(this->numberFiles); //-/ order of mappings
    //cout<<"here 2 ";

    int half = (int) ceil((double) this->numberFiles/2);
    if (this->numberFiles %2 != 0) {
        half--;
    }
    //pFS("HALF");pFS(convertNumToStr(half));
    //cout<<"here 3 ";
    // First and last file ALIGNMENT and CHECKING
    pFS("\n -- Checking Last and First Fragments ..!");
    //cout<<"here 4 ";
    this->alignThenAnalyzeFirstAndLast();
    pFS("\n -- Last and First Fragments checked..!");

    for (int i = 1; i < half; i++) {
        // ALIGNMENT
        this->changeFastQ(myArraySteps[2*i]);
        this->alignThenAnalyze(myArraySteps[2*i]);
        // CHECKING
        for (int j = 1; j < 2*i; j+=2)
            this->checkCorrespond(myArraySteps[2*i], myArraySteps[j], 2*i, j);
        pFS("\n -- -- checked Fragments for Anchoring:"+convertNumToStr(myArraySteps[2*i])+"and mates");
        // ALIGNMENT
        this->changeFastQ(myArraySteps[2*i+1]);
        this->alignThenAnalyze(myArraySteps[2*i+1]);
        // CHECKING
        for (int j = 0; j < 2*i+1; j+=2)
            this->checkCorrespond(myArraySteps[2*i+1], myArraySteps[j], 2*i+1, j);
        pFS("\n -- -- checked Fragments for Anchoring:"+convertNumToStr(myArraySteps[2*i+1])+"and mates");
    }
    if (this->numberFiles %2 != 0) {
        // ALIGNMENT
        this->changeFastQ(myArraySteps[this->numberFiles-1]);
        this->alignThenAnalyze(myArraySteps[this->numberFiles-1]);
        // CHECKING
        for (int j = 0; j < this->numberFiles-1; j++) {
            this->checkCorrespond(myArraySteps[this->numberFiles-1], myArraySteps[j], this->numberFiles-1, j);
        }
    }
    pFS("\n -- Recursive Remained Checking...");
    // proccessing reads in which lefts and rights are insconsistent
    RecursiveRemainedAligner(1,half);
    RecursiveRemainedAligner(half+1,this->numberFiles);

    delete [] myArraySteps;
    //delete [] myFragOrders;
    return;
}
void Anchoring::nextStepAlignment(int range1, int range2, int readsFilesCount) {

    int* myArraySteps = this->getMyArraySteps(this->numberFiles); //-/ order of mappings
    int half = (int) ceil((double) this->numberFiles/2);
    if (this->numberFiles %2 != 0) {
        half--;
    }
    bool corresponded = false;
    if(range1==1){
        // First and last file ALIGNMENT and CHECKING
        pFS("\n -- Checking Last and First Fragments ..!");
        this->alignThenAnalyzeFirstAndLast();
        pFS("\n -- Last and First Fragments checked..!");
    }
    int both = min(half , (range2-1)+1);
    for (int i = (range1-1+(range1==1?1:0)) ; i < both ; i++) {
        // ALIGNMENT
        this->changeFastQ(myArraySteps[2*i]);
        this->alignThenAnalyze(myArraySteps[2*i]);
        // CHECKING
        for (int j = 1; j < 2*i; j+=2)
            this->checkCorrespond(myArraySteps[2*i], myArraySteps[j], 2*i, j);
        pFS("\n -- -- checked Fragments for Anchoring:"+convertNumToStr(myArraySteps[2*i])+"and mates");
        // ALIGNMENT
        this->changeFastQ(myArraySteps[2*i+1]);
        this->alignThenAnalyze(myArraySteps[2*i+1]);
        // CHECKING
        for (int j = 0; j < 2*i+1; j+=2)
            this->checkCorrespond(myArraySteps[2*i+1], myArraySteps[j], 2*i+1, j);
        pFS("\n -- -- checked Fragments for Anchoring:"+convertNumToStr(myArraySteps[2*i+1])+"and mates");
    }
    if (this->numberFiles %2 != 0 && range2 >= half) {
        // ALIGNMENT
        pFS("\n -- Mid Fragment "+convertNumToStr(half+1)+" with all others...");
        this->changeFastQ(myArraySteps[this->numberFiles-1]);
        this->alignThenAnalyze(myArraySteps[this->numberFiles-1]);
        // CHECKING
        for (int j = 0; j < this->numberFiles-1; j++) {
            this->checkCorrespond(myArraySteps[this->numberFiles-1], myArraySteps[j], this->numberFiles-1, j);
        }
    }

    if( range2 >= half ){
        pFS("\n -- Recursive Remained Checking...");
        // proccessing reads in which lefts and rights are insconsistent
        RecursiveRemainedAligner(1,half);
        RecursiveRemainedAligner(half+1,this->numberFiles);
    }
    delete [] myArraySteps;
    //delete [] myFragOrders;
    return;
}

void Anchoring::RecursiveRemainedAligner(int first, int last){
    int sign = (last <= this->numberFiles/2 ? 1:-1 );
    if( first >= last ){

        return;
    }

    int Len = last-first+1;
    int FirstIndx = ( sign > 0 ? 2*(first-1) : (this->numberFiles-first)*2+1  );
    int LastIndx = ( sign > 0 ? 2*(last-1) : (this->numberFiles-last)*2+1  );
    int* myArraySteps = this->getMyArraySteps( this->numberFiles );

    //pFS("\n checking recremained aligner:\n");
    for (int i = 0; i < (Len/2); i++) {
        // CHECKING this Left with rights
        for (int j = 0; j < i; j++){
            //pFS("checking frag "+convertNumToStr(myArraySteps[FirstIndx + 2*i*sign])+" with frag "+convertNumToStr(myArraySteps[ LastIndx - 2*j*sign])+"for"+convertNumToStr(notAligned.size())+"\n");
            this->checkCorrespond(myArraySteps[FirstIndx + 2*i*sign], myArraySteps[ LastIndx - 2*j*sign], FirstIndx + 2*i*sign , LastIndx - 2*j*sign );
        }
        // CHECKING this Right with lefts
        for (int j = 0; j < i+1; j++){
            //pFS("checking frag "+convertNumToStr(myArraySteps[LastIndx-2*i*sign])+" with frag "+convertNumToStr(myArraySteps[2*j*sign+FirstIndx])+"for"+convertNumToStr(notAligned.size())+"\n");
            this->checkCorrespond(myArraySteps[LastIndx-2*i*sign], myArraySteps[2*j*sign+FirstIndx], LastIndx-2*i*sign , 2*j*sign + FirstIndx );

        }
    }
    //-/ remained
    delete [] myArraySteps;
    recstarts = true;
    RecursiveRemainedAligner(first,first+Len/2-1);
    RecursiveRemainedAligner(first+Len/2,last);

    return;
}


tuple<long long, long long> Anchoring::calculateSnapUnSample(long long position1, long long position2) {
    long long result = 0;
    for (long long i = 0; i < this->depth->size(); i++) {
        result += this->depth->at(i).end+this->floatingEdge - (this->depth->at(i).start-this->floatingEdge) + 1;
        if (position1 <= result) {
            if (position2 <= result) {
                return make_tuple((position1-result+this->depth->at(i).end+this->floatingEdge), (position2-result+this->depth->at(i).end+5));
            }
            long long temp = position1-result+this->depth->at(i).end+this->floatingEdge;
            for (long long j = i+1; j < this->depth->size(); j++) {
                result += this->depth->at(j).end+this->floatingEdge - (this->depth->at(j).start-this->floatingEdge) + 1;
                if (position2 <= result) {
                    return make_tuple((temp), (position2-result+this->depth->at(j).end+this->floatingEdge));
                }
            }
        }
    }
    return make_tuple(0, 0);
}


tuple<bool, long long, long long, long long> Anchoring::calculateIntronOrExonAndDistance(long long position) {
    bool isInExon = false;
    long long index = 0;
    long long fromRight = 0;
    long long fromLeft = 0;

    long long i;
    for (i = 1; i < this->depth->size(); i++) {
        if (this->depth->at(i).start <= position && position <= this->depth->at(i).end) {
            isInExon = true;
            index = i;
            fromLeft = abs(position - this->depth->at(i).start);
            fromRight = abs(position - this->depth->at(i).end);
            return make_tuple(isInExon, index, fromLeft, fromRight);
        } else if (position < this->depth->at(i).start && this->depth->at(i-1).end < position) {
            isInExon = false;
            index = i;
            fromLeft = abs(position - this->depth->at(i-1).end);
            fromRight = abs(position - this->depth->at(i).start);
            return make_tuple(isInExon, index, fromLeft, fromRight);
        }
    }

    if (position < this->depth->at(0).start) {
        isInExon = false;
        index = 0;
        fromLeft = 0;
        fromRight = abs(position - this->depth->at(0).start);
        return make_tuple(isInExon, index, fromLeft, fromRight);
    } else if (this->depth->at(0).start <= position && position <= this->depth->at(0).end) {
        isInExon = true;
        index = 0;
        fromLeft = abs(position - this->depth->at(0).start);;
        fromRight = abs(position - this->depth->at(0).end);
        return make_tuple(isInExon, index, fromLeft, fromRight);
    } else if (this->depth->at(this->depth->size()-1).end < position) {
        isInExon = false;
        index = this->depth->size();
        fromLeft = abs(position - this->depth->at(this->depth->size()-1).end);
        fromRight = 0;
        return make_tuple(isInExon, index, fromLeft, fromRight);
    }
    return make_tuple(isInExon, index, fromLeft, fromRight);
}


void Anchoring::addIntervalToDepth(long long a, long long b) {
    //http://codereview.stackexchange.com/questions/46463/adding-intervals-to-an-interval-store
    if (this->depth->size() == 0) {
        this->depth->push_back(interval());
        this->depth->at(0).start = a;
        this->depth->at(0).end = b;
        return;
    }
    for (long long i = 0; i < this->depth->size(); i++) {
        if (this->depth->at(i).start <= a && this->depth->at(i).end >= b) {
            return;
        }
        if (this->depth->at(i).start > b || this->depth->at(i).end < a) {
            continue;
        } else {
            long long tempA = this->depth->at(i).start;
            long long tempB = this->depth->at(i).end;
            this->depth->erase(this->depth->begin()+i);
            this->addIntervalToDepth(min(tempA, a), max(tempB, b));
            return;
        }
    }
    long long insert_index = 0;
    while (insert_index < this->depth->size() && this->depth->at(insert_index).start < b) {
        insert_index++;
    }
    this->depth->insert(this->depth->begin()+insert_index, interval());
    this->depth->at(insert_index).start = a;
    this->depth->at(insert_index).end = b;
    return;

}
void Anchoring::moreFragClusters_Recursive(vector<readNumAndSnap> ReadVec,int first,int last,int index1,int index2){

    int sign = (last <= this->numberFiles/2 ? 1:-1 );
    if( first == last ){
        return;
    }
    int Len = last-first+1;
    int* myArraySteps = this->getMyArraySteps( this->numberFiles );
//    int FirstIndx = ( sign > 0 ? 2*(first-1) : (this->numberFiles-first)*2+1  );
    int FirstIndx = std::distance(myArraySteps, std::find(myArraySteps, myArraySteps + this->numberFiles, first));
//    int LastIndx = ( sign > 0 ? 2*(last-1) : (this->numberFiles-last)*2+1  );
    int LastIndx = std::distance(myArraySteps, std::find(myArraySteps, myArraySteps + this->numberFiles, last));

    bool anyCoresspond = false;
    for (int i = 0; i < (Len/2); i++) {
        // CHECKING this Left with rights
        for (int j = 0; j < 2*i && !anyCoresspond ; j++){
            if(this->checkCorrespond1read(ReadVec,first+i                                                                               , last-j                                                                               ,
                                                  distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+i)) , distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, last-j)) )){
                // more frags in this scope
                anyCoresspond = true;
                if( i > 1 )// left more
                    moreFragClusters_Recursive( ReadVec,first,first+i-1,index1,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+i-1)));
                if( j > 1 )// right more
                    moreFragClusters_Recursive( ReadVec,last-j+1,last,FirstIndx + 2*(i+1)*sign,index2 );
            }
        }
        // CHECKING this Right with lefts
        for (int j = 0; j < 2*i+1 && !anyCoresspond ; j++)
            if(this->checkCorrespond1read(ReadVec,last-i,                                                                                   first+j,
                                                  distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, last-i)), distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+j)) )){
                anyCoresspond = true;
                if( j > 1 )// left more
                    moreFragClusters_Recursive( ReadVec,first,first+j-1,index1,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+j-1)));
                if( i > 1 )// right more
                    moreFragClusters_Recursive( ReadVec,last-i+1,last,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, last-i+1)),index2 );

            }
    }
    // for odd len vectors:
    if( Len%2 != 0 && !anyCoresspond ){
        // check len/2 + 1 th with all other
        for(int j=0 ; j < Len/2 && !anyCoresspond ; j++ ){
            // with left
            if(this->checkCorrespond1read(ReadVec,first+j,                                                                               (first+last)/2,
                                                  distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+j)), distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, (first+last)/2)) )
                    && !anyCoresspond){
                anyCoresspond = true;
                if( j > 1 )// left more
                    moreFragClusters_Recursive( ReadVec,first,first+j-1,index1,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, first+j-1)));
                // +and half right
                moreFragClusters_Recursive( ReadVec,(first+last)/2 + 1,last,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, (first+last)/2 +1)),index2);

            }
            // with right
            if(this->checkCorrespond1read(ReadVec,last-j,                                                                               (first+last)/2,
                                                  distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, last-j)), distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, (first+last)/2)) )
                    && !anyCoresspond){
                anyCoresspond = true;
                if( j > 1 )// left more
                    moreFragClusters_Recursive( ReadVec,last-j+1,last,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, last-j+1)),index2);
                // +and half left
                moreFragClusters_Recursive( ReadVec,first, (first+last)/2 - 1,index1,distance(myArraySteps, find(myArraySteps, myArraySteps + this->numberFiles, (first+last)/2 - 1)));

            }
        }
    }

    if(! anyCoresspond){
        int next_indx1 = std::distance(myArraySteps, std::find(myArraySteps, myArraySteps + this->numberFiles, first+Len/2-1));
        int next_indx2 = std::distance(myArraySteps, std::find(myArraySteps, myArraySteps + this->numberFiles, first+Len/2));
        moreFragClusters_Recursive(ReadVec,       first, first+Len/2-1,       index1, next_indx1 );
        moreFragClusters_Recursive(ReadVec, first+Len/2, last         ,   next_indx2, index2      );
    }
    delete [] myArraySteps;
    //
    return;
}
bool Anchoring::checkCorrespond1read(vector<readNumAndSnap> ReadVec,int step1, int step2, int index1, int index2 ){
    //-/ change needed:
    // /////////////////////////////
    // /////////////////////////////
    //
    // we need to change corresponding details for reads and add new clusters in a vector here
    // ///////////////////////////////////////////////////////

    if (step1 > step2) {
        int temp = step1;
        step1 = step2;
        step2 = temp;
        temp = index1;
        index1 = index2;
        index2 = temp;
    }
    readNumAndSnap read1;
    readNumAndSnap read2;
    if (index1 < ReadVec.size() && index2 < ReadVec.size()) {
        read1 = ReadVec.at(index1);
        read2 = ReadVec.at(index2);
        int realStep1 = correctFragNo(this->reads[read1.readNum-1].length,step1);
        int realStep2 = correctFragNo(this->reads[read1.readNum-1].length,step2);
        int readLength = this->reads[read1.readNum-1].length - this->d;
        int chunkCount = (readLength + this->d + chunkSize - 20) / chunkSize ;
        int chunkCountPrime = (readLength + chunkSize - 20) / chunkSize ;
        int steps = (int) (readLength+this->chunkSize-20)/this->chunkSize;
        int remain = (int) readLength - steps*this->chunkSize;

        bool st2IsLast = step2 == this->numberFiles;
        int anyRemain = ( st2IsLast
                          ?
                              isReverse( read1 )*remain
                            :
                              0 );

        if ( abs(read1.snap-read2.snap - pow(-1,isReverse(read1))*(this->chunkSize*(realStep1-realStep2))+anyRemain)
             < this->numGap*abs(realStep2-realStep1)
             && isInTheSameDirection(read1, read2)){
            //changeBinGenome(read1.snap,read1.snap +this->chunkSize);
            //changeBinGenome(read2.snap,read2.snap +this->chunkSize);

            // lets change here :

            ++this->correspond;
            if(this->isSample) {

                this->reads[read1.readNum-1].flagClus2.push_back(read1.flag);
                this->reads[read1.readNum-1].firstPosClus2.push_back(read1.snap - (this->d <= chunkSize/2 ?
                                                                                       this->d*(read1.flag == 1 ? 1 : -1 )
                                                                                     :
                                                                                       (chunkSize-this->d)*(read1.flag == 1 ? -1 : 1 )
                                                                                       ));
                if( st2IsLast )
                    if(chunkCount == chunkCountPrime)
                        this->reads[read1.readNum-1].lastPosClus2.push_back(read2.snap - (read2.flag == 1 ? this->d : 0 ));
                    else
                        this->reads[read1.readNum-1].lastPosClus2.push_back(read2.snap - ( (chunkSize-this->d)*(read1.flag == 1 ? -1 : 0 )  ));

                else{
                    this->reads[read1.readNum-1].lastPosClus2.push_back(read2.snap - (this->d <= chunkSize/2 ?
                                                                                          this->d*(read2.flag == 1 ? 1 : -1 )
                                                                                        :
                                                                                          (chunkSize-this->d)*(read1.flag == 1 ? -1 : 1 )
                                                                                          ));
                }


            } /*else {
                tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
                this->reads[read1.readNum-1].flag = read1.flag+3;
             }*/


            if( st2IsLast )
                this->reads[read1.readNum-1].lastFragClus2.push_back(realStep2 + (chunkCount == chunkCountPrime ? 0 : 1));
            else
                this->reads[read1.readNum-1].lastFragClus2.push_back(realStep2 + (this->d <= chunkSize/2 ? 0 : 1 ) );
                this->reads[read2.readNum-1].lastFragment = realStep2 + (this->d <= chunkSize/2 ? 0 : 1 );
            //-/c this->reads[read1.readNum-1].firstFragClus2.push_back(step1);
            this->reads[read1.readNum-1].firstFragClus2.push_back(correctFragNo(this->reads[read1.readNum-1].length,step1));
            //-/c this->reads[read1.readNum-1].lastFragClus2.push_back(step2);
            this->reads[read1.readNum-1].lastFragClus2.push_back(correctFragNo(this->reads[read1.readNum-1].length,step2));
            //this->reads[read1.readNum-1].dClus2 = this->d;
            return true;
        }
    }else pFS("E3");
    return false;
}

inline int Anchoring::correctFragNo(int readLen , int step){

    int numberOfFrags = (readLen+chunkSize-20) / chunkSize ;
    int mid = ((this->numberFiles + (this->numberFiles%2) ) / 2 ) ;
    return (step < mid ? step : numberOfFrags - (this->numberFiles-step));
}

inline void Anchoring::writeSingleUniques(vector<readNumAndSnap> ReadVec ,int bowtiedIndex, ofstream* file, long long readIndex){
    int* myStepArray = this->getMyArraySteps(this->numberFiles);
    //*file << "\n " << readIndex;
    for(int i=0 ; i < bowtiedIndex; i++)
        if( ReadVec.at(i).snap > 0 )
            // if not included in any cluster in reads[readIndex] --> can't check now, we have to check after extension loop
            *file<<readIndex<<"\t"<<correctFragNo(reads[readIndex].length,myStepArray[i])<<"\t"<<ReadVec.at(i).snap<<"\t"<<ReadVec.at(i).flag<<endl;// <<bowtiedIndex
    delete myStepArray;
}

void Anchoring::checkCorrespond(int step1, int step2, int index1, int index2, int bowtiedindex) {
    if (step1 > step2) {
        int temp = step1;
        step1 = step2;
        step2 = temp;
        temp = index1;
        index1 = index2;
        index2 = temp;
    }
    readNumAndSnap read1;
    readNumAndSnap read2;
    vector< vector<readNumAndSnap> > notAlignedInThisStep;
    string name = outputDir+"/SingleUniques.txt";
    ofstream singleUniques(name.c_str());
    int readLength = 0;
    int steps = 0;
    int remain = 0;

    for (long long i = 0; i < this->notAligned.size(); i++) {

        if (index1 < this->notAligned.at(i).size() && index2 < this->notAligned.at(i).size()) {//-//Question?Else?

            read1 = this->notAligned.at(i).at(index1);
            read2 = this->notAligned.at(i).at(index2);

            int realStep1 = correctFragNo(this->reads[read1.readNum-1].length,step1);
            int realStep2 = correctFragNo(this->reads[read1.readNum-1].length,step2);

            readLength = this->reads[read1.readNum-1].length - this->d;
            int chunkCount = (readLength + this->d + chunkSize - 20) / chunkSize ;
            int chunkCountPrime = (readLength + chunkSize - 20) / chunkSize ;
            steps = (int) (readLength+this->chunkSize-20)/this->chunkSize;
            remain = (int) readLength - steps*this->chunkSize;

            bool st2IsLast = step2 == this->numberFiles;
            int anyRemain = ( st2IsLast
                              ?
                                  isReverse( read1 )*remain
                                :
                                  0 );
            if ( abs(read1.snap-read2.snap - pow(-1,isReverse(read1))*(this->chunkSize*(realStep1-realStep2)) + anyRemain)
                < this->numGap*abs(realStep2-realStep1)
                && isInTheSameDirection(read1, read2)) {

                //changeBinGenome(read1.snap,read1.snap +this->chunkSize);
                //changeBinGenome(read2.snap,read2.snap +this->chunkSize);
                ++this->correspond;
                if(this->isSample) {
                    this->reads[read1.readNum-1].firstPosition = read1.snap - (this->d <= chunkSize/2 ?
                                                                                   this->d*(read1.flag == 1 ? 1 : -1 )
                                                                                 :
                                                                                   (chunkSize-this->d)*(read1.flag == 1 ? -1 : 1 )
                                                                                                       );
                    if( st2IsLast )
                        if(chunkCount == chunkCountPrime)
                            this->reads[read2.readNum-1].lastPosition = read2.snap - (read2.flag == 1 ? this->d : 0 );
                        else
                            this->reads[read2.readNum-1].lastPosition = read2.snap - ( (chunkSize-this->d)*(read1.flag == 1 ? -1 : 0 )  );

                    else{
                        this->reads[read2.readNum-1].lastPosition = read2.snap - (this->d <= chunkSize/2 ?
                                                                                      this->d*(read2.flag == 1 ? 1 : -1 )
                                                                                    :
                                                                                      (chunkSize-this->d)*(read1.flag == 1 ? -1 : 1 )
                                                                                      );

                    }
                    this->reads[read1.readNum-1].flag = read1.flag;
                } else {
                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
                    this->reads[read1.readNum-1].flag = read1.flag+3;
                 }

                //-/c this->reads[read1.readNum-1].firstFragment = step1;
                this->reads[read1.readNum-1].firstFragment = realStep1 + (this->d <= chunkSize/2 ? 0 : 1);
                //-/c this->reads[read1.readNum-1].lastFragment = step2;
                if( st2IsLast )
                    this->reads[read2.readNum-1].lastFragment = realStep2 + (chunkCount == chunkCountPrime ? 0 : 1);
                else
                    this->reads[read2.readNum-1].lastFragment = realStep2 + (this->d <= chunkSize/2 ? 0 : 1 );

                if(bowtiedindex != -1)
                    writeSingleUniques(this->notAligned.at(i),bowtiedindex,&singleUniques, read1.readNum-1);
                // ------------------------------------------------------------------------------------------***
                // ////////////////////////////////                    /////////////       //////////////////***
                // ////////////////////////////////     NEW CODES      /////////////  V3   //////////////////***
                // //////////////////////////////////////////////////////////////////////////////////////////***
                // /////////////////////////////// All Clusters of all reads
                // ------------------------------------------------------------------------------------------***
                //moreFragClusters_Recursive(this->notAligned.at(i),step1,step2,index1,index2);
                bool leftRemained = ( step1 > 1+1 ? true : false ); // another cluster needs at least two more frags
                bool rightRemained = ( step2 < this->numberFiles-(1) ? true : false );
                int* myStepArray = this->getMyArraySteps(this->numberFiles);

                if( leftRemained ){
                    int next_indx = std::distance(myStepArray, std::find(myStepArray, myStepArray + this->numberFiles, step1-1));
                    //if(recstarts == true) pFS("\n Clust2check for "+convertNumToStr(next_indx)+"\n");
                    moreFragClusters_Recursive(this->notAligned.at(i),      1,          step1-1,      0 , next_indx);
                }
                if( rightRemained ){
                    int next_indx = std::distance(myStepArray, std::find(myStepArray, myStepArray + this->numberFiles, step2+1));
                    //if(recstarts == true) pFS("\n Clust2check for "+convertNumToStr(next_indx)+convertNumToStr(step2+1)+"\n");
                    moreFragClusters_Recursive(this->notAligned.at(i),step2+1,this->numberFiles,next_indx, 1     );
                }

                // ------------------------------------------------------------------------------------------***
                // ////////////////////////////////                    /////////////       //////////////////***
                // ////////////////////////////////     NEW CODES      /////////////  V2   //////////////////***
                // //////////////////////////////////////////////////////////////////////////////////////////***
                // //////////////////////////////// two clusters of all reads with recursive calling
                // ------------------------------------------------------------------------------------------***
                // version 3 is better :D

                // ------------------------------------------------------------------------------------------***
                // ////////////////////////////////                    /////////////       //////////////////***
                // ////////////////////////////////     NEW CODES      /////////////  V1   //////////////////***
                // //////////////////////////////////////////////////////////////////////////////////////////***
                // //////////////////////////////// two clusters of all reads w. checking of left and right
                // ------------------------------------------------------------------------------------------***


                //-/code: to add second "fragment cluster" proccessing here:
                //-/
//                bool leftRemained = ( step1 > 1+1 ? true : false ); // ( if()? then : else )     //-/ at least two fragments for another Cluster
//                bool rightRemained = ( step2 < this->numberFiles-(1+1) ? true : false );         //-/ at least two fragments for another Cluster
//                //-/ ************* ---------------------------------------------------- ***************
//                //-/ ************* but 1 fragment is enough for entering extension step ***************
//                //-/ *************                                                      ***************
//                //-/ int leftRemained = ( step1 > 1 ? ( step1 > 2 ? 2 : 1 ): 0);

//                bool _2ndClusterFound = false;
//                // /////////////////////////                    ////////////////////////
//                // /////////////////////////   LEFT REMAINED    ////////////////////////
//                // /////////////////////////                    ////////////////////////
//                if( leftRemained ){
//                    //-/ Extension to cover repeat regions if(multipleLeft)
//                    //-/ if(still lefremained) then: check for other clusters
//                    vector<int> indices;
//                    //-/? clear ? right ?
//                    indices.clear();
//                    for( int ii = index1 % 2 ; ii < index1 ; ii+=2 )
//                        indices.push_back(ii);
//                    for( int ii = 0 , jj = indices.size()-1  ; ii < indices.size()/2 && !_2ndClusterFound; ii+=1 , jj-=1 ){
//                        for( int k = ii ; k > 0 && !_2ndClusterFound ; k-- ){
//                            //-/ check this->notAligned.at(i).at(indices(ii)) and this->notAligned.at(i).at(indices(jj+k))
//                            //-/

//                            read3 = this->notAligned.at(i).at(indices[ii]);
//                            read4 = this->notAligned.at(i).at(indices[jj+k]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*(ii - (jj+k) )))
//                                 < this->numGap*abs( (jj+k) - ii )
//                                 && isInTheSameDirection(read3, read4)){
//                                //-/ change this :
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = ii+1;
//                                this->reads[read3.readNum-1].lastFragClus2 = jj+k+1;
//                                this->reads[read3.readNum-1].d = this->d;

//                            }
//                            //-/ check this->notAligned.at(i).at(indices(ii-k)) and this->notAligned.at(i).at(indices(jj))
//                            //-/
//                            read3 = this->notAligned.at(i).at(indices[ii-k]);
//                            read4 = this->notAligned.at(i).at(indices[jj]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*( (ii-k) - jj )) && !_2ndClusterFound )
//                                 < this->numGap*abs( jj - (ii-k) ) /* && !_2ndClusterFound */
//                                 && isInTheSameDirection(read3, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = ii-k+1;
//                                this->reads[read3.readNum-1].lastFragClus2 = jj+1;
//                                this->reads[read3.readNum-1].d = this->d;
//                            }

//                        }

//                        //-/ check for ii and jj:
//                        if(!_2ndClusterFound){
//                            read3 = this->notAligned.at(i).at(indices[ii]);
//                            read4 = this->notAligned.at(i).at(indices[jj]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*(ii - jj )))
//                                 < this->numGap*abs( jj - ii )
//                                 && isInTheSameDirection(read3, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = ii+1;
//                                this->reads[read3.readNum-1].lastFragClus2 = jj+1;
//                                this->reads[read3.readNum-1].d = this->d;
//                            }
//                        }
//                    }

//                    if( ! indices.size()%2 && !_2ndClusterFound ){ // check this (mid fragment with all remaineds for LEFT)
//                        // /////////////////////////////////////////////////////////////
//                        // ///////////////// check mid fragment with all other fragments
//                        read5 = this->notAligned.at(i).at( indices[indices.size()/2] );
//                        //read5 = this->notAligned.at(i).at( indices[indices.size()/2 ]);
//                        for( int t = 0 ; t < indices.size()/2 && !_2ndClusterFound; t++ ){
//                            // t and indices.size()/2
//                            read3 = this->notAligned.at(i).at( indices[t]);
//                            read4 = this->notAligned.at(i).at( indices[indices.size()-t-1]);

//                            if ( abs(read3.snap-read5.snap - pow(-1,isReverse(read3))*(this->chunkSize*( t - indices.size()/2 )))
//                                 < this->numGap*abs( indices.size()/2 - t )
//                                 && isInTheSameDirection(read3, read5)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read5.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = t+1;
//                                this->reads[read3.readNum-1].lastFragClus2 = indices.size()/2 +1;
//                                this->reads[read3.readNum-1].d = this->d;
//                            }
//                            //  indices.size()-i-1 and indices.size()/2
//                            if ( abs(read5.snap-read4.snap - pow(-1,isReverse(read5))*(this->chunkSize*( indices.size()/2 - (indices.size()-t-1) )))
//                                 < this->numGap*abs( (indices.size()-t-1) - indices.size()/2 )
//                                 && isInTheSameDirection(read5, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read5.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read5.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = indices.size()/2 +1;
//                                this->reads[read3.readNum-1].lastFragClus2 = (indices.size()-t-1)+1 ;
//                                this->reads[read3.readNum-1].d = this->d;

//                            }


//                        }
//                        // /////////////////end->check mid fragment with all other fragments

//                    }
//                }


//                // /////////////////////////                    ////////////////////////
//                // /////////////////////////   RIGHT REMAINED   ////////////////////////
//                // /////////////////////////                    ////////////////////////
//                if( rightRemained && !_2ndClusterFound){
//                    vector<int> indices;
//                    //-/? clear ? right ?
//                    indices.clear();
//                    int readLen = this->notAligned.at(i).size();
//                    //-/ risky:
//                    // if(index2 % 2 == 0 ) index2+=3;
//                    //-/ risky:
//                    for( int ii = index2 % 2 ; ii < index2 ; ii+=2 )
//                        indices.push_back(ii);
//                    for( int ii = 0 , jj = indices.size()-1  ; ii < indices.size()/2 && !_2ndClusterFound; ii+=1 , jj-=1 ){
//                        for( int k = ii ; k > 0 && !_2ndClusterFound ; k-- ){
//                            //-/ check this->notAligned.at(i).at(indices(ii)) and this->notAligned.at(i).at(indices(jj+k))
//                            //-/
//                            read3 = this->notAligned.at(i).at(indices[ii]);
//                            read4 = this->notAligned.at(i).at(indices[jj+k]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*(ii - (jj+k) )))
//                                 < this->numGap*abs( (jj+k) - ii )
//                                 && isInTheSameDirection(read3, read4)){
//                                //-/ change this :
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = readLen -(jj+k);
//                                this->reads[read3.readNum-1].lastFragClus2 = readLen - ii;
//                                this->reads[read3.readNum-1].d = this->d;

//                            }
//                            //-/ check this->notAligned.at(i).at(indices(ii-k)) and this->notAligned.at(i).at(indices(jj))
//                            //-/
//                            read3 = this->notAligned.at(i).at(indices[ii-k]);
//                            read4 = this->notAligned.at(i).at(indices[jj]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*( (ii-k) - jj )))
//                                 < this->numGap*abs( jj - (ii-k) ) /* && !_2ndClusterFound */
//                                 && isInTheSameDirection(read3, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = readLen - jj;
//                                this->reads[read3.readNum-1].lastFragClus2 = readLen - (ii-k);
//                                this->reads[read3.readNum-1].d = this->d;
//                            }

//                        }

//                        //-/ check for ii and jj:
//                        if(!_2ndClusterFound){
//                            read3 = this->notAligned.at(i).at(indices[ii]);
//                            read4 = this->notAligned.at(i).at(indices[jj]);
//                            if ( abs(read3.snap-read4.snap - pow(-1,isReverse(read3))*(this->chunkSize*(ii - jj )))
//                                 < this->numGap*abs( jj - ii )
//                                 && isInTheSameDirection(read3, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = readLen - jj;
//                                this->reads[read3.readNum-1].lastFragClus2 = readLen - ii;
//                                this->reads[read3.readNum-1].d = this->d;
//                            }
//                        }
//                    }
//                    if( ! indices.size()%2 && !_2ndClusterFound ){ // check this (mid fragment with all remaineds for RIGHT)
//                        // /////////////////////////////////////////////////////////////
//                        // ///////////////// check mid fragment with all other fragments
//                        int readLen = this->notAligned.at(i).size();
//                        read5 = this->notAligned.at(i).at( indices[indices.size()/2 ]);
//                        for( int t = 0 ; t < indices.size()/2 && !_2ndClusterFound; t++ ){
//                            // t and indices.size()/2
//                            read3 = this->notAligned.at(i).at(indices[t]);
//                            read4 = this->notAligned.at(i).at(indices[indices.size()-t-1]);

//                            if ( abs(read3.snap-read5.snap - pow(-1,isReverse(read3))*(this->chunkSize*( t - indices.size()/2 )))
//                                 < this->numGap*abs( indices.size()/2 - t )
//                                 && isInTheSameDirection(read3, read5)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read5.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read3.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read3.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = readLen - (indices.size()/2);
//                                this->reads[read3.readNum-1].lastFragClus2 = readLen - t ;
//                                this->reads[read3.readNum-1].d = this->d;
//                            }
//                            //  indices.size()-i-1 and indices.size()/2
//                            if ( abs(read5.snap-read4.snap - pow(-1,isReverse(read5))*(this->chunkSize*( indices.size()/2 - (indices.size()-t-1) )))
//                                 < this->numGap*abs( (indices.size()-t-1) - indices.size()/2 )
//                                 && isInTheSameDirection(read5, read4)){
//                                _2ndClusterFound = true;
//                                ++this->correspondCluster2;
//                                if(this->isSample) {
//                                    this->reads[read3.readNum-1].firstPosClus2 = read4.snap;
//                                    this->reads[read3.readNum-1].lastPosClus2 = read5.snap;
//                                    this->reads[read3.readNum-1].flagClus2 = read5.flag;
//                                } /*else {
//                                                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                                                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                                                 }*/
//                                this->reads[read3.readNum-1].firstFragClus2 = readLen - (indices.size()-t-1) ;
//                                this->reads[read3.readNum-1].lastFragClus2 = readLen - indices.size()/2 ;
//                                this->reads[read3.readNum-1].d = this->d;

//                            }


//                        }
//                        // /////////////////end->check mid fragment with all other fragments

//                    }
//                }

                //-/ end of second "fragment cluster"
                // ------------------------------------------------------------------------------------------***
                // ////////////////////////////////                    //////////////////////////////////////***
                // ////////////////////////////////  END of NEW CODES  //////////////////////////////////////***
                // //////////////////////////////////////////////////////////////////////////////////////////***
                // ------------------------------------------------------------------------------------------***
                // Remained :
                //      - Add to depth for clust2 and other functions
                //      - Extension
                //      - More than two Clusters
                // Idea :
                //      - Clustering Fragments Anchors => x 100 100 300 100 x 100 400  x 400 x x
                // ------------------------------------------------------------------------------------------***
//                if (this->isSample) {
//                    if (isReverse(read1)) {
//                        this->addIntervalToDepth(read2.snap, read1.snap+this->chunkSize);
//                    } else {
//                        this->addIntervalToDepth(read1.snap, read2.snap+this->chunkSize);
//                    }
//                }

            } else {
                notAlignedInThisStep.push_back(this->notAligned.at(i));
            }
        }
    }
    singleUniques.close();
    vector< vector<readNumAndSnap> >().swap(this->notAligned);
    this->notAligned = vector< vector<readNumAndSnap> >(notAlignedInThisStep.begin(), notAlignedInThisStep.end());
    vector< vector<readNumAndSnap> >().swap(notAlignedInThisStep);
    return;
}

void Anchoring::alignThenAnalyzeFirstAndLast() {

    string command;
    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".sam 2>> "+outputDir+"log-bowtie.txt";
    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".fq -S "+outputDir+this->outputName+"_reads_d="+convertNu.mToStr(this->d)+"_"+convertNumToStr(1)+".sam 2>> "+outputDir+"log-bowtie.txt";
    system(command.c_str());
    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam 2>> "+outputDir+"log-bowtie.txt";
    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".fq -S "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam 2>> "+outputDir+"log-bowtie.txt";
    system(command.c_str());

    vector<readNumAndSnap> list1;
    vector<readNumAndSnap> list2;

    list1 = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".sam");
    list2 = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam");

    int readLength = 0;

    int steps = 0;
    int remain = 0;

    for (long long i = 0; i < list1.size(); i++) {
        readLength = this->reads[list1.at(i).readNum-1].length - this->d;
        //cout<<"RL:"<<readLength<<endl;
        steps = (int) (readLength+this->chunkSize-20)/this->chunkSize;
        remain = (int) readLength - steps*this->chunkSize;
        if( abs(list2.at(i).snap-list1.at(i).snap - pow(-1,isReverse(list1.at(i)))*(this->chunkSize*(steps-1))+isReverse(list1.at(i))*remain)
            < this->numGap*(steps-1)
            && isInTheSameDirection(list1.at(i), list2.at(i))) {
            //changeBinGenome(list1.at(i).snap,list1.at(i).snap +this->chunkSize);
            //changeBinGenome(list2.at(i).snap,list2.at(i).snap +this->chunkSize);
            ++this->correspond;
            int chunkCount = (this->reads[list1.at(i).readNum-1].length + chunkSize - 20) / chunkSize ;
            int chunkCountPrime = (readLength + chunkSize - 20) / chunkSize ;
            if(this->isSample){
                this->reads[list1.at(i).readNum-1].firstPosition = list1.at(i).snap - (this->d <= chunkSize/2 ?
                                                                                           this->d*(list1.at(i).flag == 1 ? 1 : -1 )
                                                                                         :
                                                                                           (chunkSize-this->d)*(list1.at(i).flag == 1 ? -1 : 1 ) );
                if( chunkCount == chunkCountPrime )
                    this->reads[list1.at(i).readNum-1].lastPosition = list2.at(i).snap - (list2.at(i).flag == 1 ? this->d : 0 );
                else
                    this->reads[list1.at(i).readNum-1].lastPosition = list2.at(i).snap + (list2.at(i).flag == 1 ? chunkSize - this->d : 0 ) ;


                this->reads[list1.at(i).readNum-1].flag = list1.at(i).flag;
            } else {
                tie(this->reads[list1.at(i).readNum-1].firstPosition, this->reads[list1.at(i).readNum-1].lastPosition) = this->calculateSnapUnSample(list1.at(i).snap, list2.at(i).snap);
                this->reads[list1.at(i).readNum-1].flag = list1.at(i).flag+3;
            }

            this->reads[list1.at(i).readNum-1].firstFragment = 1 + (this->d <= chunkSize/2 ? 0 : 1);
            this->reads[list1.at(i).readNum-1].lastFragment = chunkCount;

            this->reads[list1.at(i).readNum-1].d = this->d;

            if (this->isSample) {
                if (isReverse(list1.at(i))) {
                    this->addIntervalToDepth(list2.at(i).snap, list1.at(i).snap+this->chunkSize);
                } else {
                    this->addIntervalToDepth(list1.at(i).snap, list2.at(i).snap+this->chunkSize);
                }
            }

        } else {
            vector<readNumAndSnap> tempVec;
            tempVec.push_back(list1.at(i));
            tempVec.push_back(list2.at(i));
            this->notAligned.push_back(tempVec);
            vector<readNumAndSnap>().swap(tempVec);
        }
    }
    vector<readNumAndSnap>().swap(list1);
    vector<readNumAndSnap>().swap(list2);
    return;
}


void Anchoring::alignThenAnalyze(int step) {
    string command;
    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+"_ch.fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam 2>> "+outputDir+"log-bowtie.txt";
    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+"_ch.fq -S "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam 2>> "+outputDir+"log-bowtie.txt";
    system(command.c_str());

    vector<readNumAndSnap> alignsInThisStep;
    alignsInThisStep = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam");
    long long j = 0;
    for (long long i = 0; i < this->notAligned.size(); i++) {
        if( j >= alignsInThisStep.size() )
            break;
        if (this->notAligned.at(i).at(0).readNum == alignsInThisStep.at(j).readNum) {
            this->notAligned[i].push_back(alignsInThisStep.at(j));
            j++;
            if( j >= alignsInThisStep.size() )
                break;
        }
    }
    return;
}


void Anchoring::changeFastQ(int n) {
    string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(n)+".fq";
    ifstream ifstr(name.c_str());
    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(n)+"_ch.fq";
    ofstream ofstr(name.c_str());
    string line;
    string segment;
    long long lineCounter = 0;
    while (getline(ifstr, line)) {
        if (lineCounter % 4 == 0) {
            stringstream test(line);
            getline(test, segment, '_');
            getline(test, segment, '_');
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<line<<"\n";
            }
        } else if (lineCounter % 4 == 1) {
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<line<<"\n";
            }
        } else if (lineCounter % 4 == 3) {
            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
                ofstr<<"+\n"<<line<<"\n";
            }
        }
        ++lineCounter;
    }
    ifstr.close();
    ofstr.close();
    return;
}


vector<int> Anchoring::getWhichCheck(int* myArraySteps, int k) {
    vector<int> tempArray(k);
    vector<int> tempArray2(k);
    for (int j = 0; j < k; j++) {
        tempArray[j] = abs(myArraySteps[j]-myArraySteps[k]);
        tempArray2[j] = j;
    }
    sort(begin(tempArray2),end(tempArray2), [&](int i1, int i2) {return tempArray[i1] > tempArray[i2];});
    vector<int>().swap(tempArray);
    return tempArray2;
}


bool compareZippedSnapAndReadNum(readNumAndSnap r1, readNumAndSnap r2) {return r1.readNum < r2.readNum;}


vector<readNumAndSnap> Anchoring::rearrange(string filename) {
    vector<readNumAndSnap> vec;
    SamAnalyzer* analyze = new SamAnalyzer(filename);
    analyze->parseSAMFile();
    vec = analyze->getOutput();
    sort(vec.begin(), vec.end(), compareZippedSnapAndReadNum);
    delete analyze;
    return vec;
}


int Anchoring::isReverse(readNumAndSnap r) {
    int* flags = this->flagAnalyzer(r.flag);
    if (flags[4]) {
        delete [] flags;
        return 1;
    } else {
        delete [] flags;
        return 0;
    }
}


bool Anchoring::isInTheSameDirection(readNumAndSnap r1, readNumAndSnap r2) {
    int* flags1 = this->flagAnalyzer(r1.flag);
    int* flags2 = this->flagAnalyzer(r2.flag);
    if (flags1[4] == flags2[4]) {
        delete [] flags1;
        delete [] flags2;
        return true;
    } else {
        delete [] flags1;
        delete [] flags2;
        return false;
    }
}


int* Anchoring::flagAnalyzer(int flag) {
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


int* Anchoring::getMyArraySteps(int fragMax) {
    int *myArraySteps = new int[fragMax];
    for (int i = 0; i < fragMax; i+=2) {
        myArraySteps[i] = i/2+1;
    }
    for (int i = 1; i < fragMax; i+=2) {
        myArraySteps[i] = fragMax-(i/2);
    }
    return myArraySteps;
}

int* Anchoring::getFragsOrders(int fragMax) {
    int *myFragOrders = new int[fragMax];
    int j=1;
    for (int i = 0; i < fragMax/2; i++) {
        myFragOrders[i] = j++;
        myFragOrders[fragMax-1-i] = j++;
    }
    if( fragMax%2 == 1)
        myFragOrders[fragMax/2] = j;
    return myFragOrders;
}
