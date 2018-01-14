////
////  Anchoring.cpp
////  RNA-Seq
////
////  Created by Farid Rashidi on 2015-09-27.
////  Copyright © 2015 Farid Rashidi. All rights reserved.
////

//#include "Anchoring.h"


//Anchoring::Anchoring(iRead* reads, vector<interval>* depth, string referenceName, string outputName, bool isSample, int numGap, int chunkSize, int d, int floatingEdge, int dashV, int numberOfThread, string outputDir) {
//    this->reads = reads;
//    this->depth = depth;
//    this->correspond = 0;
//    this->referenceName = referenceName;
//    this->outputName = outputName;
//    this->isSample = isSample;
//    this->chunkSize = chunkSize;
//    this->d = d;
//    this->numGap = numGap;
//    this->floatingEdge = floatingEdge;
//    this->dashV = dashV;
//    this->numberOfThread = numberOfThread;
//    this->outputDir = outputDir;
//}


//Anchoring::~Anchoring() {
//    vector< vector<readNumAndSnap> >().swap(this->notAligned);
//}


//void Anchoring::prepareFiles() {
//    // taghir header read ha: @headerGhabli_tool_shomare
//    // peyda kardan bishtarin toole read ha
//    string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+".fq";
//    ifstream ifstr(name.c_str());
//    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_my.fq";
//    ofstream ofstr(name.c_str());
//    string line;
//    string privous_line;
//    int counter = 0;
//    int readLenMax = 0;
//    while (getline(ifstr, line)) {
//        if (counter % 4 == 0){
//            privous_line = line;
//        } else if (counter % 4 == 1) {
//            ofstr << privous_line << "\n";
//            ofstr << line.substr(this->d) << "\n";
//            if(readLenMax < line.substr(this->d).size()) readLenMax = (int) line.substr(this->d).size();
//        } else if (counter % 4 == 3) {
//            ofstr << "+\n" << line.substr(this->d) << "\n";
//        }
//        ++counter;
//    }
//    ifstr.close();
//    ofstr.close();

//    // tafkik read ha be chunk ha
//    this->numberFiles = (int) (readLenMax+this->chunkSize-20)/this->chunkSize;
//    int* myStepArray = this->getMyArraySteps(this->numberFiles);
//    ofstream* ofs = new ofstream[this->numberFiles];
//    for (int i = 0; i < this->numberFiles; i++) {
//        name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(i+1)+".fq";
//        ofs[i].open(name.c_str());
//    }
//    counter = 0;
//    int steps = 0;
//    int remain = 0;
//    line = "";
//    privous_line = "";
//    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_my.fq";
//    ifstream ifstr2(name.c_str());
//    while (getline(ifstr2, line)) {
//        if (counter % 4 == 0) {
//            privous_line = line;
//        } else if (counter % 4 == 1) {
//            steps = (int) (line.size()+this->chunkSize-20)/this->chunkSize;
//            remain = (int) line.size() - steps*this->chunkSize;
//            for (int i = 0; i < steps; i++) {
//                ofs[myStepArray[i]-1] << privous_line << "\n";
//                if (i%2 == 1) {
//                    if (i == 1) {
//                        ofs[myStepArray[i]-1] << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
//                    } else {
//                        ofs[myStepArray[i]-1] << line.substr(line.size()-this->chunkSize*(i+1)/2-remain, this->chunkSize) << "\n";
//                    }
//                } else {
//                    ofs[myStepArray[i]-1] << line.substr(this->chunkSize*i/2, this->chunkSize) << "\n";
//                }
//            }
//        } else if (counter % 4 == 3) {
//            for (int i = 0; i < steps; i++) {
//                if (i%2 == 1) {
//                    if (i == 1) {
//                        ofs[myStepArray[i]-1] << "+\n" << line.substr(line.size()-this->chunkSize-remain, this->chunkSize+remain) << "\n";
//                    } else {
//                        ofs[myStepArray[i]-1] << "+\n" << line.substr(line.size()-this->chunkSize*(i+1)/2-remain, this->chunkSize) << "\n";
//                    }
//                } else {
//                    ofs[myStepArray[i]-1] << "+\n" << line.substr(this->chunkSize*i/2, this->chunkSize) << "\n";
//                }
//            }
//        }
//        ++counter;
//    }
//    for (int i = 0; i < this->numberFiles; i++) ofs[i].close();
//    ifstr2.close();
//    delete[] myStepArray;
//    delete[] ofs;
//    return;
//}


//void Anchoring::createNewReads(int d) {
//    string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+".fq";
//    ifstream ifstr(name.c_str());
//    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(d)+".fq";
//    ofstream ofstr(name.c_str());
//    string line;
//    string test;
//    string segment;
//    long long lineCounter = 0;
//    while (getline(ifstr, line)) {
//        if (lineCounter % 4 == 0) {
//            stringstream test(line);
//            getline(test, segment, '_');
//            getline(test, segment, '_');
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<line<<"\n";
//            }
//        } else if (lineCounter % 4 == 1) {
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<line<<"\n";
//            }
//        } else if (lineCounter % 4 == 3) {
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<"+\n"<<line<<"\n";
//            }
//        }
//        ++lineCounter;
//    }
//    ifstr.close();
//    ofstr.close();
//    return;
//}


//void Anchoring::writingNotAligned() {
//    string name = outputDir+"all_notAlign.txt";
//    ofstream ofstr(name.c_str());
//    for (long long i = 0; i < this->notAligned.size(); i++) {
//        ofstr << this->notAligned.at(i).at(0).readNum << "\t";
//        for (int j = 0; j < this->notAligned.at(i).size(); j++) {
//            ofstr << this->notAligned.at(i).at(j).snap << "\t";
//        }
//        ofstr << endl;
//    }
//    ofstr.close();
//    return;
//}


//void Anchoring::nextStepAlignment() {
//    int** myArrayAlignmentCorrespond = new int*[this->numberFiles];
//    for (int i = 0; i < this->numberFiles; ++i) {
//        myArrayAlignmentCorrespond[i] = new int[this->numberFiles+1];
//    }
//    for(int i = 0; i<this->numberFiles; i++) {
//        for(int j=0; j<this->numberFiles; j++) {
//            if (i < j) {
//                myArrayAlignmentCorrespond[i][j]=0;
//            } else {
//                myArrayAlignmentCorrespond[i][j]=1;
//            }
//        }
//    }

//    int* myArraySteps = this->getMyArraySteps(this->numberFiles);

//    int half = (int) ceil((double) this->numberFiles/2);
//    if (this->numberFiles %2 != 0) {
//        half--;
//    }


//    // First and last file ALIGNMENT and CHECKING
//    this->alignThenAnalyzeFirstAndLast();
//    myArrayAlignmentCorrespond[0][1] = 1;


//    for (int i = 1; i < half; i++) {
//        // ALIGNMENT
//        this->changeFastQ(myArraySteps[2*i]);
//        this->alignThenAnalyze(myArraySteps[2*i]);
//        // CHECKING
//        vector<int> whichNext(2*i);
//        whichNext = this->getWhichCheck(myArraySteps, 2*i);
//        for (int j = 0; j < i; j++) {
//            this->checkCorrespond(myArraySteps[2*i], myArraySteps[whichNext[j]], 2*i, whichNext[j]);
//            myArrayAlignmentCorrespond[whichNext[j]][2*i] = 1;
//        }
//        vector<int>().swap(whichNext);


//        // ALIGNMENT
//        this->changeFastQ(myArraySteps[2*i+1]);
//        this->alignThenAnalyze(myArraySteps[2*i+1]);
//        // CHECKING
//        vector<int> whichNext2(2*i+1);
//        whichNext2 = this->getWhichCheck(myArraySteps, 2*i+1);
//        for (int j = 0; j < i+1; j++) {
//            this->checkCorrespond(myArraySteps[2*i+1], myArraySteps[whichNext2[j]], 2*i+1, whichNext2[j]);
//            myArrayAlignmentCorrespond[whichNext2[j]][2*i+1] = 1;
//        }
//        vector<int>().swap(whichNext2);

//        // REMAINED CHECKING
//        for (int k = 0; k < 2*i+1; k++) {
//            if (myArrayAlignmentCorrespond[k][2*i]!=1) {
//                this->checkCorrespond(myArraySteps[2*i], myArraySteps[k], 2*i, k);
//                myArrayAlignmentCorrespond[k][2*i] = 1;
//            }
//            if (myArrayAlignmentCorrespond[k][2*i+1]!=1) {
//                this->checkCorrespond(myArraySteps[2*i+1], myArraySteps[k], 2*i+1, k);
//                myArrayAlignmentCorrespond[k][2*i+1] = 1;
//            }
//        }
//    }


//    if (this->numberFiles %2 != 0) {
//        // ALIGNMENT
//        this->changeFastQ(myArraySteps[this->numberFiles-1]);
//        this->alignThenAnalyze(myArraySteps[this->numberFiles-1]);
//        // CHECKING
//        vector<int> whichNext3(this->numberFiles-1);
//        whichNext3 = this->getWhichCheck(myArraySteps, this->numberFiles-1);
//        for (int j = 0; j < this->numberFiles-1; j++) {
//            this->checkCorrespond(myArraySteps[this->numberFiles-1], myArraySteps[whichNext3[j]], this->numberFiles-1, whichNext3[j]);
//            myArrayAlignmentCorrespond[whichNext3[j]][this->numberFiles-1] = 1;
//        }
//        vector<int>().swap(whichNext3);
//    }

//    for (int i = 0; i < this->numberFiles; ++i) {
//        delete [] myArrayAlignmentCorrespond[i];
//    }
//    delete [] myArrayAlignmentCorrespond;
//    delete [] myArraySteps;
//    return;
//}


//void Anchoring::findExonIntron() {
//    long long counter = 0;
//    bool gotoMainLoop;

//    bool isInExon1;
//    long long index1;
//    long long fromRight1;
//    long long fromLeft1;

//    bool isInExon2;
//    long long index2;
//    long long fromRight2;
//    long long fromLeft2;

//    int* myStepArray = this->getMyArraySteps(this->numberFiles);
//    vector< vector<readNumAndSnap> > notAlignedInThisStep;

//    for (long long i = 0; i < this->notAligned.size() ; i++) {
//        gotoMainLoop = false;
//        for (int j = 0; j < this->notAligned.at(i).size() && !gotoMainLoop; j++) {
//            if(this->notAligned.at(i).at(j).snap != 0) {
//                tie(isInExon1, index1, fromLeft1, fromRight1) = this->calculateIntronOrExonAndDistance(this->notAligned.at(i).at(j).snap);
//                if (isInExon1 == true && !isReverse(this->notAligned.at(i).at(j))) {
//                    for (int k = 0; k < this->notAligned.at(i).size() && !gotoMainLoop; k++) {
//                        if (myStepArray[j] != myStepArray[k]) {
//                            if(this->notAligned.at(i).at(k).snap != 0) {
//                                tie(isInExon2, index2, fromLeft2, fromRight2) = this->calculateIntronOrExonAndDistance(this->notAligned.at(i).at(k).snap);
//                                if (myStepArray[j] < myStepArray[k]) {
//                                    if(isInExon2 == false
//                                       && index2 == index1+1
//                                       && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))
//                                       && fromRight1 < 150) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 1;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        this->addIntervalToDepth(this->notAligned.at(i).at(k).snap, this->notAligned.at(i).at(k).snap+this->chunkSize);

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                } else {
//                                    if(isInExon2 == false
//                                       && index2 == index1
//                                       && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))
//                                       && fromLeft1 < 150) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 1;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        this->addIntervalToDepth(this->notAligned.at(i).at(k).snap, this->notAligned.at(i).at(k).snap+this->chunkSize);

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                }

//                            }
//                        }
//                    }

//                } else if (isInExon1 == true && isReverse(this->notAligned.at(i).at(j))) {
//                    for (int k = 0; k < this->notAligned.at(i).size() && !gotoMainLoop; k++) {
//                        if (myStepArray[j] != myStepArray[k]) {
//                            if(this->notAligned.at(i).at(k).snap != 0) {
//                                tie(isInExon2, index2, fromLeft2, fromRight2) = this->calculateIntronOrExonAndDistance(this->notAligned.at(i).at(k).snap);
//                                if (myStepArray[j] < myStepArray[k]) {
//                                    if(isInExon2 == false
//                                       && index2 == index1
//                                       && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))
//                                       && fromLeft1 < 120) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 17;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        this->addIntervalToDepth(this->notAligned.at(i).at(k).snap, this->notAligned.at(i).at(k).snap+this->chunkSize);

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                } else {
//                                    if(isInExon2 == false
//                                       && index2 == index1+1
//                                       && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))
//                                       && fromRight1 < 120) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 17;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        this->addIntervalToDepth(this->notAligned.at(i).at(k).snap, this->notAligned.at(i).at(k).snap+this->chunkSize);

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                }


//                            }
//                        }
//                    }
//                }
//            }
//        }
//        if (!gotoMainLoop) {
//            notAlignedInThisStep.push_back(this->notAligned.at(i));
//        }
//    }

////    cout << "******************** E_I:      d = " << this->d << " : " << counter << "\n";

//    vector< vector<readNumAndSnap> >().swap(this->notAligned);
//    this->notAligned = vector< vector<readNumAndSnap> >(notAlignedInThisStep.begin(), notAlignedInThisStep.end());
//    vector< vector<readNumAndSnap> >().swap(notAlignedInThisStep);
//    delete[] myStepArray;
//    return;
//}


//void Anchoring::findTwoAdjacentExon() {
//    long long counter = 0;
//    bool gotoMainLoop;

//    bool isInExon1;
//    long long index1;
//    long long fromRight1;
//    long long fromLeft1;

//    bool isInExon2;
//    long long index2;
//    long long fromRight2;
//    long long fromLeft2;

//    int* myStepArray = this->getMyArraySteps(this->numberFiles);
//    vector< vector<readNumAndSnap> > notAlignedInThisStep;

//    for (long long i = 0; i < this->notAligned.size() ; i++) {
//        gotoMainLoop = false;
//        for (int j = 0; j < this->notAligned.at(i).size() && !gotoMainLoop; j++) {
//            if(this->notAligned.at(i).at(j).snap != 0) {
//                tie(isInExon1, index1, fromLeft1, fromRight1) = this->calculateIntronOrExonAndDistance(this->notAligned.at(i).at(j).snap);
//                for (int k = j+1; k < this->notAligned.at(i).size() && !gotoMainLoop; k++) {
//                    if(this->notAligned.at(i).at(k).snap != 0) {
//                        tie(isInExon2, index2, fromLeft2, fromRight2) = this->calculateIntronOrExonAndDistance(this->notAligned.at(i).at(k).snap);

//                        if (isInExon1 == true && isInExon2 == true) {
//                            if (!isReverse(this->notAligned.at(i).at(j))) {
//                                if (myStepArray[j] < myStepArray[k]) {
//                                    if (abs(abs(myStepArray[j]-myStepArray[k])*this->chunkSize - (fromRight1 + fromLeft2)) < 10
//                                        && index1+1 == index2
//                                        && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 2;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                } else {
//                                    if (abs(abs(myStepArray[j]-myStepArray[k])*this->chunkSize - (fromLeft1 + fromRight2)) < 10
//                                        && index1 == index2+1
//                                        && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 2;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                }

//                            } else {
//                                if (myStepArray[j] < myStepArray[k]) {
//                                    if (abs(abs(myStepArray[j]-myStepArray[k])*this->chunkSize - (fromLeft1 + fromRight2)) < 10
//                                        && index1 == index2+1
//                                        && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 18;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                } else {
//                                    if (abs(abs(myStepArray[j]-myStepArray[k])*this->chunkSize - (fromRight1 + fromLeft2)) < 10
//                                        && index1+1 == index2
//                                        && isInTheSameDirection(this->notAligned.at(i).at(j), this->notAligned.at(i).at(k))) {

//                                        ++this->correspond;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstPosition = this->notAligned.at(i).at(k).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastPosition = this->notAligned.at(i).at(j).snap;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].firstFragment = myStepArray[k];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].lastFragment = myStepArray[j];
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].flag = 18;
//                                        this->reads[this->notAligned.at(i).at(j).readNum-1].d = this->d;

//                                        gotoMainLoop = true;
//                                        counter++;
//                                    }
//                                }

//                            }
//                        }
//                    }
//                }
//            }
//        }
//        if (!gotoMainLoop) {
//            notAlignedInThisStep.push_back(this->notAligned.at(i));
//        }
//    }

////    cout << "-------------------- E_E:      d = " << this->d << " : " << counter << "\n";

//    vector< vector<readNumAndSnap> >().swap(this->notAligned);
//    this->notAligned = vector< vector<readNumAndSnap> >(notAlignedInThisStep.begin(), notAlignedInThisStep.end());
//    vector< vector<readNumAndSnap> >().swap(notAlignedInThisStep);
//    delete[] myStepArray;
//    return;
//}




//tuple<long long, long long> Anchoring::calculateSnapUnSample(long long position1, long long position2) {
//    long long result = 0;
//    for (long long i = 0; i < this->depth->size(); i++) {
//        result += this->depth->at(i).end+this->floatingEdge - (this->depth->at(i).start-this->floatingEdge) + 1;
//        if (position1 <= result) {
//            if (position2 <= result) {
//                return make_tuple((position1-result+this->depth->at(i).end+this->floatingEdge), (position2-result+this->depth->at(i).end+5));
//            }
//            long long temp = position1-result+this->depth->at(i).end+this->floatingEdge;
//            for (long long j = i+1; j < this->depth->size(); j++) {
//                result += this->depth->at(j).end+this->floatingEdge - (this->depth->at(j).start-this->floatingEdge) + 1;
//                if (position2 <= result) {
//                    return make_tuple((temp), (position2-result+this->depth->at(j).end+this->floatingEdge));
//                }
//            }
//        }
//    }
//    return make_tuple(0, 0);
//}


//tuple<bool, long long, long long, long long> Anchoring::calculateIntronOrExonAndDistance(long long position) {
//    bool isInExon = false;
//    long long index = 0;
//    long long fromRight = 0;
//    long long fromLeft = 0;

//    long long i;
//    for (i = 1; i < this->depth->size(); i++) {
//        if (this->depth->at(i).start <= position && position <= this->depth->at(i).end) {
//            isInExon = true;
//            index = i;
//            fromLeft = abs(position - this->depth->at(i).start);
//            fromRight = abs(position - this->depth->at(i).end);
//            return make_tuple(isInExon, index, fromLeft, fromRight);
//        } else if (position < this->depth->at(i).start && this->depth->at(i-1).end < position) {
//            isInExon = false;
//            index = i;
//            fromLeft = abs(position - this->depth->at(i-1).end);
//            fromRight = abs(position - this->depth->at(i).start);
//            return make_tuple(isInExon, index, fromLeft, fromRight);
//        }
//    }

//    if (position < this->depth->at(0).start) {
//        isInExon = false;
//        index = 0;
//        fromLeft = 0;
//        fromRight = abs(position - this->depth->at(0).start);
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    } else if (this->depth->at(0).start <= position && position <= this->depth->at(0).end) {
//        isInExon = true;
//        index = 0;
//        fromLeft = abs(position - this->depth->at(0).start);;
//        fromRight = abs(position - this->depth->at(0).end);
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    } else if (this->depth->at(this->depth->size()-1).end < position) {
//        isInExon = false;
//        index = this->depth->size();
//        fromLeft = abs(position - this->depth->at(this->depth->size()-1).end);
//        fromRight = 0;
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    }
//    return make_tuple(isInExon, index, fromLeft, fromRight);
//}


//void Anchoring::addIntervalToDepth(long long a, long long b) {
//    //http://codereview.stackexchange.com/questions/46463/adding-intervals-to-an-interval-store
//    if (this->depth->size() == 0) {
//        this->depth->push_back(interval());
//        this->depth->at(0).start = a;
//        this->depth->at(0).end = b;
//        return;
//    }
//    for (long long i = 0; i < this->depth->size(); i++) {
//        if (this->depth->at(i).start <= a && this->depth->at(i).end >= b) {
//            return;
//        }
//        if (this->depth->at(i).start > b || this->depth->at(i).end < a) {
//            continue;
//        } else {
//            long long tempA = this->depth->at(i).start;
//            long long tempB = this->depth->at(i).end;
//            this->depth->erase(this->depth->begin()+i);
//            this->addIntervalToDepth(min(tempA, a), max(tempB, b));
//            return;
//        }
//    }
//    long long insert_index = 0;
//    while (insert_index < this->depth->size() && this->depth->at(insert_index).start < b) {
//        insert_index++;
//    }
//    this->depth->insert(this->depth->begin()+insert_index, interval());
//    this->depth->at(insert_index).start = a;
//    this->depth->at(insert_index).end = b;
//    return;

//}


//void Anchoring::checkCorrespond(int step1, int step2, int index1, int index2) {
//    if (step1 > step2) {
//        int temp = step1;
//        step1 = step2;
//        step2 = temp;
//        temp = index1;
//        index1 = index2;
//        index2 = temp;
//    }
//    readNumAndSnap read1;
//    readNumAndSnap read2;
//    vector< vector<readNumAndSnap> > notAlignedInThisStep;


//    for (long long i = 0; i < this->notAligned.size(); i++) {

//        if (index1 < this->notAligned.at(i).size() && index2 < this->notAligned.at(i).size()) {
//            read1 = this->notAligned.at(i).at(index1);
//            read2 = this->notAligned.at(i).at(index2);

//            if ( abs(read1.snap-read2.snap - pow(-1,isReverse(read1))*(this->chunkSize*(step1-step2)))
//                < this->numGap*abs(step2-step1)
//                && isInTheSameDirection(read1, read2)) {

//                ++this->correspond;
//                if(this->isSample) {
//                    this->reads[read1.readNum-1].firstPosition = read1.snap;
//                    this->reads[read1.readNum-1].lastPosition = read2.snap;
//                    this->reads[read1.readNum-1].flag = read1.flag;
//                } else {
//                    tie(this->reads[read1.readNum-1].firstPosition, this->reads[read1.readNum-1].lastPosition) = this->calculateSnapUnSample(read1.snap, read2.snap);
//                    this->reads[read1.readNum-1].flag = read1.flag+3;
//                }
//                this->reads[read1.readNum-1].firstFragment = step1;
//                this->reads[read1.readNum-1].lastFragment = step2;
//                this->reads[read1.readNum-1].d = this->d;

//                if (this->isSample) {
//                    if (isReverse(read1)) {
//                        this->addIntervalToDepth(read2.snap, read1.snap+this->chunkSize);
//                    } else {
//                        this->addIntervalToDepth(read1.snap, read2.snap+this->chunkSize);
//                    }
//                }

//            } else {
//                notAlignedInThisStep.push_back(this->notAligned.at(i));
//            }
//        }
//    }

//    vector< vector<readNumAndSnap> >().swap(this->notAligned);
//    this->notAligned = vector< vector<readNumAndSnap> >(notAlignedInThisStep.begin(), notAlignedInThisStep.end());
//    vector< vector<readNumAndSnap> >().swap(notAlignedInThisStep);
//    return;
//}


//void Anchoring::alignThenAnalyzeFirstAndLast() {
//    string command;
//    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".fq -S "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    system(command.c_str());
//    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".fq -S "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    system(command.c_str());

//    vector<readNumAndSnap> list1;
//    vector<readNumAndSnap> list2;

//    list1 = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(1)+".sam");
//    list2 = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(this->numberFiles)+".sam");

//    int readLength = 0;

//    int steps = 0;
//    int remain = 0;

//    for (long long i = 0; i < list1.size(); i++) {
//        readLength = this->reads[list1.at(i).readNum-1].length - this->d;
//        steps = (int) (readLength+this->chunkSize-20)/this->chunkSize;
//        remain = (int) readLength - steps*this->chunkSize;

//        if (abs( list2.at(i).snap - list1.at(i).snap -
//                pow(-1,isReverse(list1.at(i)))*(readLength-this->chunkSize)
//                -remain*(1-isReverse(list1.at(i))) )
//            < this->numGap*(steps-1)
//            && isInTheSameDirection(list1.at(i), list2.at(i))) {

//            ++this->correspond;
//            if(this->isSample) {
//                this->reads[list1.at(i).readNum-1].firstPosition = list1.at(i).snap;
//                this->reads[list1.at(i).readNum-1].lastPosition = list2.at(i).snap;
//                this->reads[list1.at(i).readNum-1].flag = list1.at(i).flag;
//            } else {
//                tie(this->reads[list1.at(i).readNum-1].firstPosition, this->reads[list1.at(i).readNum-1].lastPosition) = this->calculateSnapUnSample(list1.at(i).snap, list2.at(i).snap);
//                this->reads[list1.at(i).readNum-1].flag = list1.at(i).flag+3;
//            }


//            this->reads[list1.at(i).readNum-1].firstFragment = 1;
//            this->reads[list1.at(i).readNum-1].lastFragment = this->numberFiles;

//            this->reads[list1.at(i).readNum-1].d = this->d;

//            if (this->isSample) {
//                if (isReverse(list1.at(i))) {
//                    this->addIntervalToDepth(list2.at(i).snap, list1.at(i).snap+this->chunkSize);
//                } else {
//                    this->addIntervalToDepth(list1.at(i).snap, list2.at(i).snap+this->chunkSize);
//                }
//            }

//        } else {
//            vector<readNumAndSnap> tempVec;
//            tempVec.push_back(list1.at(i));
//            tempVec.push_back(list2.at(i));
//            this->notAligned.push_back(tempVec);
//            vector<readNumAndSnap>().swap(tempVec);
//        }
//    }
//    vector<readNumAndSnap>().swap(list1);
//    vector<readNumAndSnap>().swap(list2);
//    return;
//}


//void Anchoring::alignThenAnalyze(int step) {
//    string command;
//    command = string(BOWTIE) + " -S "+this->referenceName+" "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+"_ch.fq -v "+convertNumToStr(this->dashV)+" -m 1 -t --mm -p "+convertNumToStr(this->numberOfThread)+" > "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    //command = string(BOWTIE2ADDRESS) + " -t --threads 4 --local --reorder -N 1 -L 20 -x CHR18 -U "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+"_ch.fq -S "+outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam 2>> "+outputDir+"log-bowtie.txt";
//    system(command.c_str());

//    vector<readNumAndSnap> alignsInThisStep;
//    alignsInThisStep = this->rearrange(outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(step)+".sam");
//    long long j = 0;
//    for (long long i = 0; i < this->notAligned.size(); i++) {
//        if (this->notAligned.at(i).at(0).readNum == alignsInThisStep.at(j).readNum) {
//            this->notAligned[i].push_back(alignsInThisStep.at(j));
//            j++;
//        }
//    }
//    return;
//}


//void Anchoring::changeFastQ(int n) {
//    string name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(n)+".fq";
//    ifstream ifstr(name.c_str());
//    name = outputDir+this->outputName+"_reads_d="+convertNumToStr(this->d)+"_"+convertNumToStr(n)+"_ch.fq";
//    ofstream ofstr(name.c_str());
//    string line;
//    string test;
//    string segment;
//    long long lineCounter = 0;
//    while (getline(ifstr, line)) {
//        if (lineCounter % 4 == 0) {
//            stringstream test(line);
//            getline(test, segment, '_');
//            getline(test, segment, '_');
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<line<<"\n";
//            }
//        } else if (lineCounter % 4 == 1) {
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<line<<"\n";
//            }
//        } else if (lineCounter % 4 == 3) {
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                ofstr<<"+\n"<<line<<"\n";
//            }
//        }
//        ++lineCounter;
//    }
//    ifstr.close();
//    ofstr.close();
//    return;
//}


//vector<int> Anchoring::getWhichCheck(int* myArraySteps, int k) {
//    vector<int> tempArray(k);
//    vector<int> tempArray2(k);
//    for (int j = 0; j < k; j++) {
//        tempArray[j] = abs(myArraySteps[j]-myArraySteps[k]);
//        tempArray2[j] = j;
//    }
//    sort(begin(tempArray2),end(tempArray2), [&](int i1, int i2) {return tempArray[i1] > tempArray[i2];});
//    vector<int>().swap(tempArray);
//    return tempArray2;
//}


//bool compareZippedSnapAndReadNum(readNumAndSnap r1, readNumAndSnap r2) {return r1.readNum < r2.readNum;}


//vector<readNumAndSnap> Anchoring::rearrange(string filename) {
//    vector<readNumAndSnap> vec;
//    SamAnalyzer* analyze = new SamAnalyzer(filename);
//    analyze->parseSAMFile();
//    vec = analyze->getOutput();
//    sort(vec.begin(), vec.end(), compareZippedSnapAndReadNum);
//    delete analyze;
//    return vec;
//}


//int Anchoring::isReverse(readNumAndSnap r) {
//    int* flags = this->flagAnalyzer(r.flag);
//    if (flags[4]) {
//        delete [] flags;
//        return 1;
//    } else {
//        delete [] flags;
//        return 0;
//    }
//}


//bool Anchoring::isInTheSameDirection(readNumAndSnap r1, readNumAndSnap r2) {
//    int* flags1 = this->flagAnalyzer(r1.flag);
//    int* flags2 = this->flagAnalyzer(r2.flag);
//    if (flags1[4] == flags2[4]) {
//        delete [] flags1;
//        delete [] flags2;
//        return true;
//    } else {
//        delete [] flags1;
//        delete [] flags2;
//        return false;
//    }
//}


//int* Anchoring::flagAnalyzer(int flag) {
//    int * flags =  new int[8];
//    for (int i=0;i<8;i++){flags[i] = 0;}
//    for (int i = 7; i >= 0; i--) {

//        if (flag >= (int)pow(2, i)) {
//            flags[i] = 1;
//            flag -= (int)pow(2, i) ;
//        }
//    }
//    return flags;
//}


//int* Anchoring::getMyArraySteps(int fragMax) {
//    int *myArraySteps = new int[fragMax];
//    for (int i = 0; i < fragMax; i+=2) {
//        myArraySteps[i] = i/2+1;
//    }
//    for (int i = 1; i < fragMax; i+=2) {
//        myArraySteps[i] = fragMax-(i/2);
//    }
//    return myArraySteps;
//}

