////
////  Searching.cpp
////  RNA-Seq
////
////  Created by Farid Rashidi on 10/19/15.
////  Copyright © 2015 Farid Rashidi. All rights reserved.
////

//#include "Searching.h"


//Searching::Searching(iRead* reads, vector<interval>* depth, long long numberOfReads, int characterPerRow, string genomeFileName, int chunkSize, int indel, int mismatch, int indelScore, int mismatchScore, int matchScore) {
//    this->reads = reads;
//    this->depth = depth;
//    this->numberOfReads = numberOfReads;
//    this->characterPerRow = characterPerRow;
//    this->genomeFileName = genomeFileName;
//    this->chunkSize = chunkSize;
//    this->indel = indel;
//    this->mismatch = mismatch;
//    this->matchScore = matchScore;
//    this->mismatchScore = mismatchScore;
//    this->indelScore = indelScore;
//}


//Searching::~Searching() {

//}


//void Searching::searchingWithCutting3(int floatingEdge) {

//    string referenceSeqs = "";

//    for (long long i = 0; i < this->depth->size(); i++) {

//        if (i % 500 == 0) {
//            cout << i << endl;
//        }


//        //cout << "((((((((((((( LEFT SEQ )))))))))))))   " << i << endl;
//        referenceSeqs = "";
//        if (i == 0) {
//            referenceSeqs += this->getGenome(this->depth->at(0).start-1-3*floatingEdge, 4*floatingEdge);
//        } else if (i == 1) {
//            referenceSeqs += this->getGenome(this->depth->at(0).end-1-3*floatingEdge, 6*floatingEdge);
//            referenceSeqs += this->getGenome(this->depth->at(1).start-1-3*floatingEdge, 4*floatingEdge);
//        } else {
//            for (long long j = i-2; j <= i-1; j++) {
//                referenceSeqs += this->getGenome(this->depth->at(j).end-1-3*floatingEdge, 6*floatingEdge);
//            }
//            referenceSeqs += this->getGenome(this->depth->at(i).start-1-3*floatingEdge, 4*floatingEdge);
//        }


//        vector<readSplicedInformation> myAnswer1;
//        myAnswer1 = this->amirHossein4(referenceSeqs, this->depth->at(i).leftSeq, 3, floatingEdge, 3, true, i);

//        for (int j = 0; j < myAnswer1.size(); j++) {
//            string cigar = convertNumToStr(myAnswer1.at(j).positionInGenome);
//            cigar += "_" + convertNumToStr(myAnswer1.at(j).positionInRead);
//            cigar += "_" + myAnswer1.at(j).cigar;
//            this->reads[this->depth->at(i).leftRead.at(j).index-1].cigar = cigar;
//            smallExon(cigar, i, j, true);//###########################################

//        }
//        vector<readSplicedInformation>().swap(myAnswer1);


//        //cout << "((((((((((((( RIGHT SEQ )))))))))))))   " << i << endl;
//        referenceSeqs = "";
//        if (i == this->depth->size()-1) {
//            referenceSeqs += this->getGenome(this->depth->at(i).end-1-floatingEdge, 4*floatingEdge);
//        } else if (i == this->depth->size()-2) {
//            referenceSeqs += this->getGenome(this->depth->at(i).end-1-floatingEdge, 6*floatingEdge);
//            referenceSeqs += this->getGenome(this->depth->at(i+1).start-1-3*floatingEdge, 4*floatingEdge);
//        } else {
//            referenceSeqs += this->getGenome(this->depth->at(i).end-1-floatingEdge, 4*floatingEdge);
//            for (long long j = i+1; j <= i+2; j++) {
//                referenceSeqs += this->getGenome(this->depth->at(j).start-1-3*floatingEdge, 6*floatingEdge);
//            }
//        }

//        vector<readSplicedInformation> myAnswer2;
//        myAnswer2 = this->amirHossein4(referenceSeqs, this->depth->at(i).rightSeq, 3, floatingEdge, 3, false, i);
//        for (int j = 0; j < myAnswer2.size(); j++) {
//            int start = 0;
//            if (this->reads[this->depth->at(i).rightRead.at(j).index-1].flag < 6) {
//                start = (this->depth->at(i).rightRead.at(j).fragment-1)*this->chunkSize + this->reads[this->depth->at(i).rightRead.at(j).index-1].d;
//            } else {
//                start = this->reads[this->depth->at(i).rightRead.at(j).index-1].length-(this->depth->at(i).rightRead.at(j).fragment*this->chunkSize);
//            }
//            myAnswer2.at(j).positionInRead += start;
//        }
//        for (int j = 0; j < myAnswer2.size(); j++) {
//            string cigar = convertNumToStr(myAnswer2.at(j).positionInGenome);
//            cigar += "_" + convertNumToStr(myAnswer2.at(j).positionInRead);
//            cigar += "_" + myAnswer2.at(j).cigar;
//            this->reads[this->depth->at(i).rightRead.at(j).index-1].cigar = cigar;

//        }
//        vector<readSplicedInformation>().swap(myAnswer2);

//    }

//    return;
//}


//void Searching::smallExon(string cigar, long long i, int j, bool isLeft) {
//    string referenceSeqs = "";
//    string readSeq = "";
//    stringstream test(cigar);
//    string positionInGenomeString;
//    string positionInReadString;
//    string input2;
//    string thirdPart;

//    getline(test, positionInGenomeString, '_');
//    long long positionInGenome = stoll(positionInGenomeString);
//    getline(test, positionInReadString, '_');
//    int positionInRead = stoi(positionInReadString);
//    getline(test, input2, '_');

//    int readLength = 0;
//    int genomeLength = 0;
//    int intronLength = 0;
//    int smallPart = 0;
//    vector<int> lengthVec;
//    if (isLeft) {
//        int p = 0;
//        for (int k = 0; k < input2.length(); k++) {
//            switch (input2[k]) {
//                case 'M':
//                    readLength += stoi(input2.substr(p, k-p));
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'm':
//                    readLength += stoi(input2.substr(p, k-p));
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'd':
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'i':
//                    readLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'N':
//                    //genomeLength += stoi(input2.substr(p, k-p));
//                    intronLength = stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case '-': {
//                    smallPart = stoi(input2.substr(p, k-p));
//                    long long a = positionInGenome+genomeLength;
//                    long long b = a + intronLength;
//                    for (long long m = i-2; m < i+1; m++) {
//                        if (a < this->depth->at(m).start && this->depth->at(m).end < b) {
//                            referenceSeqs += this->getGenome(a-1, (int)(this->depth->at(m).start-a));
//                            lengthVec.push_back((int)(this->depth->at(m).start-a));
//                            a = this->depth->at(m).end;
//                        }
//                    }
//                    if (b > a) {
//                        referenceSeqs += this->getGenome(a-1, (int)(b-a));
//                        lengthVec.push_back((int)(b-a));
//                    }
//                    readSeq = this->depth->at(i).leftSeq[j].substr(positionInRead+readLength, smallPart);
//                    p = k+1;
//                    //LocalAligner...
//                    //vector<readSplicedInformation> myAnswer1;
//                    //myAnswer1 = this->amirHossein4(referenceSeqs, this->depth->at(i).leftSeq, 3, floatingEdge, 3, true, i);
//                    break;
//                }
//                default:
//                    break;
//            }
//        }
//    }
//    else{
//        int p = 0;
//        for (int k = 0; k < input2.length(); k++) {
//            switch (input2[k]) {
//                case 'M':
//                    readLength += stoi(input2.substr(p, k-p));
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'm':
//                    readLength += stoi(input2.substr(p, k-p));
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'd':
//                    genomeLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'i':
//                    readLength += stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case 'N':
//                    //genomeLength += stoi(input2.substr(p, k-p));
//                    intronLength = stoi(input2.substr(p, k-p));
//                    p = k+1;
//                    break;
//                case '-': {
//                    smallPart = stoi(input2.substr(p, k-p));
//                    long long a = positionInGenome+genomeLength;
//                    long long b = a + intronLength;
//                    for (long long m = i; m < i+3; m++) {
//                        if (a < this->depth->at(m).start && this->depth->at(m).end < b) {
//                            referenceSeqs += this->getGenome(a-1, (int)(this->depth->at(m).start-a));
//                            lengthVec.push_back((int)(this->depth->at(m).start-a));
//                            a = this->depth->at(m).end;
//                        }
//                    }
//                    if (b > a) {
//                        referenceSeqs += this->getGenome(a-1, (int)(b-a));
//                        lengthVec.push_back((int)(b-a));
//                    }
//                    readSeq = this->depth->at(i).leftSeq[j].substr(positionInRead+readLength, smallPart);
//                    p = k+1;
//                    //localAligner
//                    vector<int > localAligned;
//                    if (readSeq.length() < 7) {
//                        vector<int > localAligned;
//                        //smithWaterman
//                        //start va en va score va in ke to kodom part oftade ro mide yani introni ke az entehaye exon localAligned[5]+i shoro beshe.
//                        //bayad be depth ezafe she
//                        //agar reverse bashe bas read barax she
//                        localAligned = localAlignerForSearch(referenceSeqs, readSeq, lengthVec);
//                    }
//                    else{
//                        vector<fragInfromation> localAligned;
//                        //dar in ja ye vector az fraginformation darim ke chizaye mokhtalef darim,localAligned[i].isX+i shoro mishe
//                        //bayad be depth ezafe she
//                        //age reverse bashe bas read barax she
//                        localAligned = amirHossein5(referenceSeqs,readSeq,3,lengthVec);



//                    }
//                    //vector<readSplicedInformation> myAnswer1;
//                    //myAnswer1 = this->amirHossein4(referenceSeqs, this->depth->at(i).leftSeq, 3, floatingEdge, 3, true, i);
//                    break;
//                }
//                default:
//                    break;
//            }
//        }
//    }

//}


//bool compareFragInformationScore(fragInfromation r1, fragInfromation r2) {return r1.score > r2.score;}


//bool compareFragInformationAnswer1(fragInfromation r1, fragInfromation r2) {return r1.answer1 > r2.answer1;}


//vector<int> Searching::localAlignerForSearch(string seq_a, string seq_b , vector <int> lengthVec){

//    int mu    = 5;//mismatch
//    int delta = 3;//indel
//    int ind;

//    int N_a = (int)seq_a.length();                     // get the actual lengths of the sequences
//    int N_b = (int)seq_b.length();

//    // initialize H
//    double H[N_a+1][N_b+1];
//    for(int i=0;i<=N_a;i++){
//        for(int j=0;j<=N_b;j++){
//            H[i][j]=0.;
//        }
//    }

//    double temp[4];
//    int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking

//    // here comes the actual algorithm

//    for(int i=1;i<=N_a;i++){
//        for(int j=1;j<=N_b;j++){
//            temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1],mu);
//            temp[1] = H[i-1][j]-delta;
//            temp[2] = H[i][j-1]-delta;
//            temp[3] = 0.;
//            H[i][j] = find_array_max(temp,4,&ind);
//            switch(ind){
//                case 0:                                  // score in (i,j) stems from a match/mismatch
//                    I_i[i][j] = i-1;
//                    I_j[i][j] = j-1;
//                    break;
//                case 1:                                  // score in (i,j) stems from a deletion in sequence A
//                    I_i[i][j] = i-1;
//                    I_j[i][j] = j;
//                    break;
//                case 2:                                  // score in (i,j) stems from a deletion in sequence B
//                    I_i[i][j] = i;
//                    I_j[i][j] = j-1;
//                    break;
//                case 3:                                  // (i,j) is the beginning of a subsequence
//                    I_i[i][j] = i;
//                    I_j[i][j] = j;
//                    break;
//            }
//        }
//    }

//    // search H for the maximal score
//    double H_max = 0.;
//    int i_max=0,j_max=0;
//    for(int i=1;i<=N_a;i++){
//        for(int j=1;j<=N_b;j++){
//            if(H[i][j]>H_max){
//                H_max = H[i][j];
//                i_max = i;
//                j_max = j;
//            }
//        }
//    }

//    //cout<<H_max<<endl;

//    // Backtracking from H_max
//    int endReference = i_max,startReference=0;
//    int endRead = j_max,startRead = 0;
//    int current_i=i_max,current_j=j_max;
//    int next_i=I_i[current_i][current_j];
//    int next_j=I_j[current_i][current_j];
//    int tick=0;
//    char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

//    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j!=0) && (next_i!=0)){

//        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
//        else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A

//        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
//        else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B

//        current_i = next_i;
//        current_j = next_j;
//        next_i = I_i[current_i][current_j];
//        next_j = I_j[current_i][current_j];
//        tick++;
//    }


//    startRead = current_j;
//    startReference = current_i;
//    int sum =0;
//    int i;
//    for (i=0; i<lengthVec.size(); i++) {

//        sum+=lengthVec[i];

//        if (startReference < sum && endReference >sum) {
//            if (endReference - sum < sum - startReference) {
//                endReference = sum;
//                break;
//            }
//            else{
//                startReference = sum;
//                i+=1;
//                break;
//            }
//        }
//        if (endReference < sum) {
//            break;
//        }
//    }
//    vector<int> output;
//    output.push_back(startRead);
//    output.push_back(endRead);
//    output.push_back(startReference);
//    output.push_back(endReference);
//    output.push_back(H_max);
//    output.push_back(i);
//    return output;

//}


//double Searching::similarity_score(char a, char b, int mu) {

//    double result;
//    if(a==b){
//        result=1.;
//    }
//    else{
//        result=-mu;
//    }
//    return result;
//}


//double Searching::find_array_max(double array[], int length, int* ind) {

//    double max = array[0];            // start with max = first element
//    ind = 0;

//    for(int i = 1; i<length; i++){
//        if(array[i] > max){
//            max = array[i];
//            ind = &i;
//        }
//    }
//    return max;                    // return highest value in array
//}


//vector<readSplicedInformation> Searching::amirHossein4(string s1, vector<string> s, int numberOfChar, int floatingEdge, int numberOfExon, bool isLeft, long long indexExon) {
//    // s1 reference
//    // s2 read

//    int minLength = 7;
//    int match = 9;
//    int mismatch = 6;
//    int gap_open = 5;


//    string str1;
//    for(int i = 0; i < s1.length(); i++) {
//        if (s1[i] == 'A' || s1[i] == 'a')
//            str1 += '0';
//        if (s1[i] == 'T' || s1[i] == 't')
//            str1 += '1';
//        if (s1[i] == 'C' || s1[i] == 'c')
//            str1 += '2';
//        if (s1[i] == 'G' || s1[i] == 'g')
//            str1 += '3';
//    }
//    //###################################
//    int referenceDepth[480] = {};
//    vector<vector<int>> index1;
//    for (int i = 0; i < pow(4,numberOfChar); i++) {
//        vector<int> tempVec1;
//        index1.push_back(tempVec1);
//        vector<int>().swap(tempVec1);
//    }
//    string tempString;
//    int tempInt = 0;
//    for (int i = 0; i < str1.length()-numberOfChar+1; i++) {
//        for(int j = 0; j < numberOfChar; j++) {
//            tempString = str1[i+j];
//            tempInt += atoi(tempString.c_str())*pow(4, numberOfChar-j-1);
//        }
//        index1[tempInt].push_back(i);
//        tempInt = 0;
//    }

//    int leng1 = 0;
//    int leng2 = 0;
//    int score = 0;
//    string cigar = "";

//    vector<readSplicedInformation> results;

//    for (int z = 0; z < s.size(); z++) {
//        string s2 = s.at(z);
//        string str2;
//        for(int i = 0; i < s2.length(); i++) {
//            if (s2[i] == 'A' || s2[i] == 'a')
//                str2 += '0';
//            if (s2[i] == 'T' || s2[i] == 't')
//                str2 += '1';
//            if (s2[i] == 'C' || s2[i] == 'c')
//                str2 += '2';
//            if (s2[i] == 'G' || s2[i] == 'g')
//                str2 += '3';
//        }
//        vector<fragInfromation> answer;
//        vector<vector<int>> index2;
//        tempInt = 0;
//        for(int i = 0; i < str2.length()-numberOfChar+1; i++) {
//            for(int j = 0; j < numberOfChar; j++) {
//                tempString = str2[i+j];
//                tempInt += atoi(tempString.c_str())*pow(4, numberOfChar-j-1);
//            }
//            index2.push_back(index1[tempInt]);
//            tempInt = 0;
//        }

//        for (int i = 0; i < index2.size(); i++) {
//            if (index2.at(i).size() == 0) {
//                continue;
//            }
//            for (int j = 0; j < index2.at(i).size(); j++) {
//                int answer1 = i;
//                int answer2 = index2.at(i).at(j);
//                int m = i;
//                int n = j;
//                score  = numberOfChar*match;
//                for (int n = 0; n < numberOfChar; n++) {
//                    cigar += "M";
//                }
//                leng1 = numberOfChar;
//                leng2 = numberOfChar;
//                if (i < index2.size() - 1) {
//                    while (1) {
//                        if ((m+1) <= (index2.size()-1)) {
//                            long pos = find(index2.at(m+1).begin(), index2.at(m+1).end(), index2.at(m).at(n)+1) - index2.at(m+1).begin();
//                            if (pos < index2.at(m+1).size()) {
//                                m += 1;
//                                n = (int) pos;
//                                score += match;
//                                cigar += "M";
//                                leng1 += 1;
//                                leng2 += 1;
//                                continue;
//                            }
//                        }
//                        if ((m+numberOfChar+1) <= (index2.size()-1)) { //mismatch
//                            long pos = find(index2.at(m+numberOfChar+1).begin(), index2.at(m+numberOfChar+1).end(), index2.at(m).at(n)+numberOfChar+1) - index2.at(m+numberOfChar+1).begin();
//                            if (pos < index2.at(m+numberOfChar+1).size()) {
//                                m += numberOfChar+1;
//                                n = (int) pos;
//                                score += numberOfChar*match-1*mismatch;
//                                cigar += "m";
//                                for (int n = 0; n < numberOfChar; n++) {
//                                    cigar += "M";
//                                }
//                                leng1 += numberOfChar+1;
//                                leng2 += numberOfChar+1;
//                                continue;
//                            }
//                        }
//                        bool goToWhile = false;
//                        for (int k = 1; k < numberOfChar+1 && !goToWhile; k++) {
//                            if ((m+k) <= (index2.size()-1)) { //delition : reference yedone bishtar dare
//                                long pos = find(index2.at(m+k).begin(), index2.at(m+k).end(), index2.at(m).at(n)+k+1) - index2.at(m+k).begin();
//                                if (pos < index2.at(m+k).size()) {
//                                    m += k;
//                                    n = (int) pos;
//                                    score += k*match - gap_open;
//                                    cigar = cigar.substr(0, cigar.size()-(numberOfChar - k));
//                                    cigar += "d";
//                                    for (int n = 0; n < numberOfChar; n++) {
//                                        cigar += "M";
//                                    }
//                                    leng1 += k+1;
//                                    leng2 += k;
//                                    goToWhile = true;
//                                }
//                            }
//                        }
//                        if (goToWhile) {
//                            continue;
//                        }
//                        goToWhile = false;
//                        for (int k = 1; k < numberOfChar+1 && !goToWhile; k++) {
//                            if ((m+k+1) <= (index2.size()-1)) { //insertion
//                                long pos = find(index2.at(m+k+1).begin(), index2.at(m+k+1).end(), index2.at(m).at(n)+k) - index2.at(m+k+1).begin();
//                                if (pos < index2.at(m+k+1).size()) {
//                                    m += k+1;
//                                    n = (int) pos;
//                                    score += k*match - gap_open;
//                                    cigar = cigar.substr(0, cigar.size()-(numberOfChar - k));
//                                    cigar += "i";
//                                    for (int n = 0; n < numberOfChar; n++) {
//                                        cigar += "M";
//                                    }
//                                    leng1 += k;
//                                    leng2 += k+1;
//                                    goToWhile = true;
//                                }
//                            }
//                        }
//                        if (goToWhile) {
//                            continue;
//                        }
//                        break;
//                    }
//                    if (leng1 > minLength && leng2 > minLength) {
//                        bool isOK = true;
//                        int chek_counter=0;
//                        while (chek_counter < answer.size()) {
//                            if (answer2 >= answer[chek_counter].answer2 && answer2+leng2 <= answer[chek_counter].answer2+answer[chek_counter].len2 && answer1 >= answer[chek_counter].answer1 && answer1+leng1 <= answer[chek_counter].answer1+answer[chek_counter].len1) {
//                                if (score < answer[chek_counter].score){
//                                    isOK = false;
//                                    break;
//                                }
//                            }
//                            if (answer2 <= answer[chek_counter].answer2 && answer2+leng2 >= answer[chek_counter].answer2+answer[chek_counter].len2 && answer1 <= answer[chek_counter].answer1 && answer1+leng1 >= answer[chek_counter].answer1+answer[chek_counter].len1) {
//                                if (score > answer[chek_counter].score){
//                                    answer.erase(answer.begin()+chek_counter);
//                                    chek_counter-=1;
//                                }
//                            }
//                            chek_counter+=1;
//                        }

//                        bool isX = false;
//                        if (isLeft) {
//                            if (answer2 > s1.size()-4*floatingEdge) {
//                                isX = true;
//                                long long position = 0;
//                                //****************************************************************************
//                                if (this->reads[this->depth->at(indexExon).leftRead.at(z).index-1].flag % 100 < 10) {
//                                    position = this->reads[this->depth->at(indexExon).leftRead.at(z).index-1].firstPosition;
//                                } else {
//                                    position = this->reads[this->depth->at(indexExon).leftRead.at(z).index-1].lastPosition;
//                                }
//                                if (abs(abs((int)(this->depth->at(indexExon).start - position) - ((int)s1.size()-floatingEdge-answer2))-answer1) > 5) {
//                                    isOK = false;
//                                }

//                            }
//                        } else {
//                            if (answer2 < 4*floatingEdge) {
//                                isX = true;
//                                long long position = 0;
//                                if (this->reads[this->depth->at(indexExon).rightRead.at(z).index-1].flag % 100 < 10) {
//                                    position = this->reads[this->depth->at(indexExon).rightRead.at(z).index-1].lastPosition+this->chunkSize;
//                                } else {
//                                    position = this->reads[this->depth->at(indexExon).rightRead.at(z).index-1].firstPosition+this->chunkSize;
//                                }
//                                if (abs(abs((int)(this->depth->at(indexExon).start - position) - (floatingEdge-answer2))-answer1) > 5) {
//                                    isOK = false;
//                                }
//                            }
//                        }

//                        if (isOK) {
//                            answer.push_back(fragInfromation());
//                            answer.back().answer1 = answer1;
//                            answer.back().answer2 = answer2;
//                            answer.back().len1 = leng1;
//                            answer.back().len2 = leng2;
//                            answer.back().score = score;
//                            answer.back().isX = isX;
//                            answer.back().cigar = this->changeCigar(cigar);
//                        }
//                    }
//                }
//            }
//        }
//        vector<vector<int>>().swap(index2);

//        int number1 = 0;
//        int number2 = 0;
//        long long positionInGenome = 0;
//        string cigar = "";


//        vector<fragInfromation> myAnswer;
//        if (answer.size() == 0) {
//            results.push_back(readSplicedInformation());
//            results.back().positionInGenome = positionInGenome;
//            results.back().positionInRead = 0;
//            results.back().cigar = cigar;
//            continue;
//        }

//        sort(answer.begin(), answer.end(), compareFragInformationScore);

//        myAnswer.push_back(answer.at(0));
//        int start1 = answer.at(0).answer1;
//        int start2 = answer.at(0).answer2;
//        int end1 = answer.at(0).answer1 + answer.at(0).len1;
//        int end2 = answer.at(0).answer2 + answer.at(0).len2;
//        int threshold = 3;
//        bool isX = answer.at(0).isX;

//        for (int i = 1; i < answer.size(); i++) {
//            if (answer.at(i).answer1 + answer.at(i).len1 < start1+threshold && answer.at(i).answer2 + answer.at(i).len2 < start2+threshold) {
//                if (isX){
//                    start1 = answer.at(i).answer1;
//                    start2 = answer.at(i).answer2;
//                    isX = answer.at(i).isX;
//                    myAnswer.push_back(answer.at(i));
//                }
//                else if(answer.at(i).isX){
//                    start1 = answer.at(i).answer1;
//                    start2 = answer.at(i).answer2;
//                    myAnswer.push_back(answer.at(i));
//                }
//            } else if (answer.at(i).answer1 > end1-threshold && answer.at(i).answer2 > end2-threshold) {
//                if (isX){
//                    end1 = answer.at(i).answer1 + answer.at(i).len1;
//                    end2 = answer.at(i).answer2 + answer.at(i).len2;
//                    isX = answer.at(i).isX;
//                    myAnswer.push_back(answer.at(i));
//                }
//                else if(answer.at(i).isX){
//                    end1 = answer.at(i).answer1 + answer.at(i).len1;
//                    end2 = answer.at(i).answer2 + answer.at(i).len2;
//                    myAnswer.push_back(answer.at(i));
//                }
//            }
//        }
//        //#####################################

//        sort(myAnswer.begin(), myAnswer.end(), compareFragInformationAnswer1);

//        for (int i=0; i < myAnswer.size(); i++) {
//            for (int j=myAnswer[i].answer2; j < myAnswer[i].answer2+myAnswer[i].len2; j++) {
//                referenceDepth[j]+=myAnswer[i].score;
//            }
//        }

//        if (myAnswer.size() > 0) {
//            if (isLeft) {
//                tie(number1, number2) = this->findPosition(myAnswer.at(0).answer2, isLeft, numberOfExon, floatingEdge);
//                if (number1 == numberOfExon-1){
//                    //-3*floatingedge ro kardam +floatingedge
//                    positionInGenome = this->depth->at(indexExon).start + floatingEdge + number2;
//                } else {
//                    // - 3*floatingedge ro +kardam
//                    positionInGenome = this->depth->at(indexExon - numberOfExon + 1 + number1).end + 3*floatingEdge + number2;
//                }
//            } else {
//                tie(number1, number2) = this->findPosition(myAnswer.at(0).answer2, isLeft, numberOfExon, floatingEdge);
//                if (number1 == 0){
//                    //-1 ro +3* kardam
//                    positionInGenome = this->depth->at(indexExon).end + 3*floatingEdge + number2;
//                } else {
//                    //-3* ro +3* kardam
//                    positionInGenome = this->depth->at(indexExon - numberOfExon + 1 + number1).start + 3*floatingEdge + number2;
//                }
//            }
//            for (int u = 0; u < myAnswer.size()-1; u++) {
//                int distanceInRead = myAnswer.at(u+1).answer1-myAnswer.at(u).answer1-myAnswer.at(u).len1;
//                int distanceInReference = myAnswer.at(u+1).answer2-myAnswer.at(u).answer2-myAnswer.at(u).len2;
//                cigar += myAnswer.at(u).cigar;
//                if (abs(distanceInRead-distanceInReference) < threshold) {
//                    string referenceTemp = s1.substr(myAnswer.at(u).answer2+myAnswer.at(u).len2,distanceInReference);
//                    string readTemp = s.at(z).substr(myAnswer.at(u).answer1+myAnswer.at(u).len1,distanceInRead);
//                    cigar += this->globalAlign(referenceTemp, readTemp);
//                } else {
//                    //?????????????????????????????????????????????????????????????
//                    if (distanceInRead >= 7) {
//                        this->reads[this->depth->at(indexExon).rightRead.at(z).index-1].flag += 100;
//                    }
//                    cigar += convertNumToStr(distanceInReference) + "N" + convertNumToStr(distanceInRead) + "-";
//                }
//            }
//            cigar += myAnswer.at(myAnswer.size()-1).cigar;
//        }

//        results.push_back(readSplicedInformation());
//        results.back().positionInGenome = positionInGenome;
//        results.back().positionInRead = myAnswer.at(0).answer1;
//        results.back().cigar = cigar;

//        vector<fragInfromation>().swap(answer);
//        vector<fragInfromation>().swap(myAnswer);

//    }
//    //####################################################
//    /*vector<int> shift;
//     if (isLeft) {
//     for (int i=0; i < numberOfExon; i++) {
//     int j=0;
//     int depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1];
//     if (i==numberOfExon-1) {
//     while (depthScore >0) {
//     j-=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j+1);
//     }
//     else{
//     while (depthScore >0) {
//     j+=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j-1);
//     }

//     }
//     }
//     else{
//     for (int i=0; i < numberOfExon; i++) {
//     int j=0;
//     int depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1];
//     if (i==numberOfExon-1) {
//     while (depthScore >0) {
//     j-=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j+1);
//     }
//     else{
//     while (depthScore >0) {
//     j+=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j-1);
//     }

//     }
//     }
//     */

//    vector<vector<int>>().swap(index1);
//    return results;
//}


//vector<fragInfromation> Searching::amirHossein5(string s1, string s2, int numberOfChar , vector <int> lengthVec) {
//    // s1 reference
//    // s2 read

//    int minLength = 7;
//    int match = 9;
//    int mismatch = 6;
//    int gap_open = 5;


//    string str1;
//    for(int i = 0; i < s1.length(); i++) {
//        if (s1[i] == 'A' || s1[i] == 'a')
//            str1 += '0';
//        if (s1[i] == 'T' || s1[i] == 't')
//            str1 += '1';
//        if (s1[i] == 'C' || s1[i] == 'c')
//            str1 += '2';
//        if (s1[i] == 'G' || s1[i] == 'g')
//            str1 += '3';
//    }
//    //###################################
//    vector<vector<int>> index1;
//    for (int i = 0; i < pow(4,numberOfChar); i++) {
//        vector<int> tempVec1;
//        index1.push_back(tempVec1);
//        vector<int>().swap(tempVec1);
//    }
//    string tempString;
//    int tempInt = 0;
//    for (int i = 0; i < str1.length()-numberOfChar+1; i++) {
//        for(int j = 0; j < numberOfChar; j++) {
//            tempString = str1[i+j];
//            tempInt += atoi(tempString.c_str())*pow(4, numberOfChar-j-1);
//        }
//        index1[tempInt].push_back(i);
//        tempInt = 0;
//    }

//    int leng1 = 0;
//    int leng2 = 0;
//    int score = 0;
//    string cigar = "";


//    string str2;
//    for(int i = 0; i < s2.length(); i++) {
//        if (s2[i] == 'A' || s2[i] == 'a')
//            str2 += '0';
//        if (s2[i] == 'T' || s2[i] == 't')
//            str2 += '1';
//        if (s2[i] == 'C' || s2[i] == 'c')
//            str2 += '2';
//        if (s2[i] == 'G' || s2[i] == 'g')
//            str2 += '3';
//    }
//    vector<fragInfromation> answer;
//    vector<vector<int>> index2;
//    tempInt = 0;
//    for(int i = 0; i < str2.length()-numberOfChar+1; i++) {
//        for(int j = 0; j < numberOfChar; j++) {
//            tempString = str2[i+j];
//            tempInt += atoi(tempString.c_str())*pow(4, numberOfChar-j-1);
//        }
//        index2.push_back(index1[tempInt]);
//        tempInt = 0;
//    }

//    for (int i = 0; i < index2.size(); i++) {
//        if (index2.at(i).size() == 0) {
//            continue;
//        }
//        for (int j = 0; j < index2.at(i).size(); j++) {
//            int answer1 = i;
//            int answer2 = index2.at(i).at(j);
//            int m = i;
//            int n = j;
//            score  = numberOfChar*match;
//            for (int n = 0; n < numberOfChar; n++) {
//                cigar += "M";
//            }
//            leng1 = numberOfChar;
//            leng2 = numberOfChar;
//            if (i < index2.size() - 1) {
//                while (1) {
//                    if ((m+1) <= (index2.size()-1)) {
//                        long pos = find(index2.at(m+1).begin(), index2.at(m+1).end(), index2.at(m).at(n)+1) - index2.at(m+1).begin();
//                        if (pos < index2.at(m+1).size()) {
//                            m += 1;
//                            n = (int) pos;
//                            score += match;
//                            cigar += "M";
//                            leng1 += 1;
//                            leng2 += 1;
//                            continue;
//                        }
//                    }
//                    if ((m+numberOfChar+1) <= (index2.size()-1)) { //mismatch
//                        long pos = find(index2.at(m+numberOfChar+1).begin(), index2.at(m+numberOfChar+1).end(), index2.at(m).at(n)+numberOfChar+1) - index2.at(m+numberOfChar+1).begin();
//                        if (pos < index2.at(m+numberOfChar+1).size()) {
//                            m += numberOfChar+1;
//                            n = (int) pos;
//                            score += numberOfChar*match-1*mismatch;
//                            cigar += "m";
//                            for (int n = 0; n < numberOfChar; n++) {
//                                cigar += "M";
//                            }
//                            leng1 += numberOfChar+1;
//                            leng2 += numberOfChar+1;
//                            continue;
//                        }
//                    }
//                    bool goToWhile = false;
//                    for (int k = 1; k < numberOfChar+1 && !goToWhile; k++) {
//                        if ((m+k) <= (index2.size()-1)) { //delition : reference yedone bishtar dare
//                            long pos = find(index2.at(m+k).begin(), index2.at(m+k).end(), index2.at(m).at(n)+k+1) - index2.at(m+k).begin();
//                            if (pos < index2.at(m+k).size()) {
//                                m += k;
//                                n = (int) pos;
//                                score += k*match - gap_open;
//                                cigar = cigar.substr(0, cigar.size()-(numberOfChar - k));
//                                cigar += "d";
//                                for (int n = 0; n < numberOfChar; n++) {
//                                    cigar += "M";
//                                }
//                                leng1 += k+1;
//                                leng2 += k;
//                                goToWhile = true;
//                            }
//                        }
//                    }
//                    if (goToWhile) {
//                        continue;
//                    }
//                    goToWhile = false;
//                    for (int k = 1; k < numberOfChar+1 && !goToWhile; k++) {
//                        if ((m+k+1) <= (index2.size()-1)) { //insertion
//                            long pos = find(index2.at(m+k+1).begin(), index2.at(m+k+1).end(), index2.at(m).at(n)+k) - index2.at(m+k+1).begin();
//                            if (pos < index2.at(m+k+1).size()) {
//                                m += k+1;
//                                n = (int) pos;
//                                score += k*match - gap_open;
//                                cigar = cigar.substr(0, cigar.size()-(numberOfChar - k));
//                                cigar += "i";
//                                for (int n = 0; n < numberOfChar; n++) {
//                                    cigar += "M";
//                                }
//                                leng1 += k;
//                                leng2 += k+1;
//                                goToWhile = true;
//                            }
//                        }
//                    }
//                    if (goToWhile) {
//                        continue;
//                    }
//                    break;
//                }
//                if (leng1 > minLength && leng2 > minLength) {
//                    bool isOK = true;
//                    int chek_counter=0;
//                    while (chek_counter < answer.size()) {
//                        if (answer2 >= answer[chek_counter].answer2 && answer2+leng2 <= answer[chek_counter].answer2+answer[chek_counter].len2 && answer1 >= answer[chek_counter].answer1 && answer1+leng1 <= answer[chek_counter].answer1+answer[chek_counter].len1) {
//                            if (score < answer[chek_counter].score){
//                                isOK = false;
//                                break;
//                            }
//                        }
//                        if (answer2 <= answer[chek_counter].answer2 && answer2+leng2 >= answer[chek_counter].answer2+answer[chek_counter].len2 && answer1 <= answer[chek_counter].answer1 && answer1+leng1 >= answer[chek_counter].answer1+answer[chek_counter].len1) {
//                            if (score > answer[chek_counter].score){
//                                answer.erase(answer.begin()+chek_counter);
//                                chek_counter-=1;
//                            }
//                        }
//                        chek_counter+=1;
//                    }
//                    int isX = 0;
//                    int sum =0;
//                    int startReference = answer2;
//                    int endReference = answer2+leng2;
//                    int r;
//                    for (r=0; r<lengthVec.size(); r++) {

//                        sum+=lengthVec[r];

//                        if (startReference < sum && endReference >sum) {
//                            if (endReference - sum < sum - startReference) {
//                                endReference = sum;
//                                break;
//                            }
//                            else{
//                                startReference = sum;
//                                r+=1;
//                                break;
//                            }
//                        }
//                        if (endReference < sum) {
//                            break;
//                        }
//                    }
//                    isX=r;
//                    answer.push_back(fragInfromation());
//                    answer.back().answer1 = answer1;
//                    answer.back().answer2 = answer2;
//                    answer.back().len1 = leng1;
//                    answer.back().len2 = leng2;
//                    answer.back().score = score;
//                    answer.back().isX = isX;
//                    //answer.back().cigar = this->changeCigar(cigar);
//                }
//            }
//        }
//    }
//    vector<vector<int>>().swap(index2);



//    sort(answer.begin(), answer.end(), compareFragInformationScore);


//    vector<fragInfromation> myAnswer;
//    myAnswer.push_back(answer.at(0));
//    int start1 = answer.at(0).answer1;
//    int start2 = answer.at(0).answer2;
//    int end1 = answer.at(0).answer1 + answer.at(0).len1;
//    int end2 = answer.at(0).answer2 + answer.at(0).len2;
//    int threshold = 3;
//    bool isX = answer.at(0).isX;

//    for (int i = 1; i < answer.size(); i++) {
//        if (answer.at(i).answer1 + answer.at(i).len1 < start1+threshold && answer.at(i).answer2 + answer.at(i).len2 < start2+threshold) {
//            if (isX >= answer.at(i).isX){
//                start1 = answer.at(i).answer1;
//                start2 = answer.at(i).answer2;
//                isX = answer.at(i).isX;
//                myAnswer.push_back(answer.at(i));
//            }
//        } else if (answer.at(i).answer1 > end1-threshold && answer.at(i).answer2 > end2-threshold) {
//            if (isX <= answer.at(i).isX){
//                end1 = answer.at(i).answer1 + answer.at(i).len1;
//                end2 = answer.at(i).answer2 + answer.at(i).len2;
//                isX = answer.at(i).isX;
//                myAnswer.push_back(answer.at(i));
//            }
//        }
//    }
//    //#####################################

//    sort(myAnswer.begin(), myAnswer.end(), compareFragInformationAnswer1);


//    /*if (myAnswer.size() > 0) {
//     for (int u = 0; u < myAnswer.size()-1; u++) {
//     int distanceInRead = myAnswer.at(u+1).answer1-myAnswer.at(u).answer1-myAnswer.at(u).len1;
//     int distanceInReference = myAnswer.at(u+1).answer2-myAnswer.at(u).answer2-myAnswer.at(u).len2;
//     cigar += myAnswer.at(u).cigar;
//     if (abs(distanceInRead-distanceInReference) < threshold) {
//     string referenceTemp = s1.substr(myAnswer.at(u).answer2+myAnswer.at(u).len2,distanceInReference);
//     string readTemp = s.at(z).substr(myAnswer.at(u).answer1+myAnswer.at(u).len1,distanceInRead);
//     cigar += this->globalAlign(referenceTemp, readTemp);
//     } else {
//     //?????????????????????????????????????????????????????????????
//     if (distanceInRead >= 7) {
//     this->reads[this->depth->at(indexExon).rightRead.at(z).index-1].flag += 100;
//     }
//     cigar += this->convertNumToStr(distanceInReference) + "N" + this->convertNumToStr(distanceInRead) + "-";
//     }
//     }
//     cigar += myAnswer.at(myAnswer.size()-1).cigar;
//     }*/


//    vector<fragInfromation>().swap(answer);
//    //vector<fragInfromation>().swap(myAnswer);



//    //####################################################
//    /*vector<int> shift;
//     if (isLeft) {
//     for (int i=0; i < numberOfExon; i++) {
//     int j=0;
//     int depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1];
//     if (i==numberOfExon-1) {
//     while (depthScore >0) {
//     j-=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j+1);
//     }
//     else{
//     while (depthScore >0) {
//     j+=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j-1);
//     }

//     }
//     }
//     else{
//     for (int i=0; i < numberOfExon; i++) {
//     int j=0;
//     int depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1];
//     if (i==numberOfExon-1) {
//     while (depthScore >0) {
//     j-=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j+1);
//     }
//     else{
//     while (depthScore >0) {
//     j+=1;
//     depthScore = referenceDepth[480-(6*i+1)*floatingEdge-1+j];
//     }
//     shift.push_back(j-1);
//     }

//     }
//     }
//     */

//    vector<vector<int>>().swap(index1);
//    return myAnswer;

//}


//string Searching::changeCigar(string cigarTemp) {
//    string result = "";
//    char temp = cigarTemp[0];
//    int counter = 0;
//    for (int i = 0; i < cigarTemp.size(); i++) {
//        if (temp == cigarTemp[i]) {
//            counter++;
//        } else {
//            result += convertNumToStr(counter) + temp;
//            counter = 1;
//        }
//        temp = cigarTemp[i];
//        if (i == cigarTemp.size()-1) {
//            result += convertNumToStr(counter) + temp;
//        }
//    }
//    return result;
//}


//tuple<int, int> Searching::findPosition(int position, bool isLeft, int size, int floatingEdge) {
//    int sum = 0;
//    for (int i = 0; i < size; i++) {
//        if (position > sum) {
//            if(i==0 && !isLeft){
//                sum += 4*floatingEdge;
//            }
//            else if(i==size-1 && isLeft){
//                sum += 4*floatingEdge;
//            }
//            else{
//                sum += 6*floatingEdge;
//            }
//        } else {
//            //be jaye i , i-1 gozashtam
//            return make_tuple(i-1, position-sum);
//        }
//    }
//    return make_tuple(size-1, position-sum);
//}


//string Searching::getTheChunkRead(long long readNum, int chunkNum) {
//    string chunkACGT = "";
//    string name = "reads.fq";
//    ifstream ifstr(name.c_str());
//    string line;
//    long long lineCounter = 0;
//    int remain = 0;
//    int steps = 0;

//    while (getline(ifstr, line)) {
//        if (lineCounter % 4 == 1) {
//            if ((lineCounter-1)/4 == readNum-1) {
//                line = line.substr(this->reads[readNum-1].d);
//                steps = (int) (line.size()+this->chunkSize-20)/this->chunkSize;
//                remain = (int) line.size() - steps*this->chunkSize;
//                if (chunkNum == steps) {
//                    chunkACGT = line.substr((chunkNum-1)*this->chunkSize,this->chunkSize+remain);
//                } else {
//                    chunkACGT = line.substr((chunkNum-1)*this->chunkSize,this->chunkSize);
//                }
//                return chunkACGT;
//            }
//        }
//        lineCounter++;
//    }

//    return chunkACGT;
//}


//void Searching::prepareForSearchingWithCutting2() {
//    string chunkACGT = "";
//    string name = "reads.fq";
//    ifstream ifstr(name.c_str());
//    string line;
//    string tempLine;
//    long long lineCounter = 0;
//    long long i = 0;
//    int remain = 0;
//    int steps = 0;

//    bool isInExon;
//    long long indexNumber;
//    long long fromRight;
//    long long fromLeft;

//    int temp;

//    while (getline(ifstr, line)) {
//        if (lineCounter % 4 == 1) {
//            i = (lineCounter-1)/4;


//            if (this->reads[i].flag == 0) {
//                int lastSteps = (int) (this->reads[i].length-this->reads[i].d+this->chunkSize-20)/this->chunkSize;
//                if (this->reads[i].firstFragment != 1 && lastSteps != this->reads[i].lastFragment) {
//                    tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[i].firstPosition);

//                    // Forward Left
//                    if (fromLeft < this->chunkSize*(this->reads[i].firstFragment-1) /*&& fromLeft % this->chunkSize != 0*/) {
//                        temp = this->reads[i].firstFragment - (int)(fromLeft)/this->chunkSize - 1;
//                        chunkACGT = "";
//                        tempLine = line.substr(this->reads[i].d);
//                        chunkACGT += line.substr(0, this->reads[i].d);
//                        this->depth->at(indexNumber).leftRead.push_back(onInterval());
//                        this->depth->at(indexNumber).leftRead.back().index = this->reads[i].index;
//                        this->depth->at(indexNumber).leftRead.back().fragment = temp;
//                        this->depth->at(indexNumber).leftRead.back().distance = (int) fromLeft - (reads[i].firstFragment-temp-1)*chunkSize;
//                        for (int iter = 1; iter <= temp; iter++) {
//                            chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize);
//                        }
//                        this->depth->at(indexNumber).leftSeq.push_back(string(chunkACGT));
//                    }

//                    // Forward Right
//                    if (fromRight < this->reads[i].length-this->reads[i].d-this->chunkSize*(this->reads[i].firstFragment-1) /*&& fromRight % this->chunkSize != 0*/) {
//                        temp = this->reads[i].firstFragment + (int)(fromRight)/this->chunkSize;
//                        tempLine = line.substr(this->reads[i].d);
//                        steps = (int) (tempLine.size()+this->chunkSize-20)/this->chunkSize;
//                        remain = (int) tempLine.size() - steps*this->chunkSize;
//                        if (temp > steps) {
//                            temp = steps;
//                        }
//                        this->depth->at(indexNumber).rightRead.push_back(onInterval());
//                        this->depth->at(indexNumber).rightRead.back().index = this->reads[i].index;
//                        this->depth->at(indexNumber).rightRead.back().fragment = temp;
//                        this->depth->at(indexNumber).rightRead.back().distance = (int) fromRight - (temp-reads[i].firstFragment)*chunkSize;
//                        chunkACGT = "";
//                        for (int iter = temp; iter <= steps; iter++) {
//                            if (iter == steps) {
//                                chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize+remain);
//                            } else {
//                                chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize);
//                            }
//                        }
//                        this->depth->at(indexNumber).rightSeq.push_back(string(chunkACGT));
//                    }
//                }
//            }




//            else if (this->reads[i].flag == 16) {
//                int lastSteps = (int) (this->reads[i].length-this->reads[i].d+this->chunkSize-20)/this->chunkSize;
//                if (this->reads[i].firstFragment != 1 && lastSteps != this->reads[i].lastFragment) {
//                    tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[i].firstPosition);

//                    // Reverse Left
//                    if (fromLeft < this->reads[i].length-this->reads[i].d-this->chunkSize*(this->reads[i].firstFragment) /*&& fromLeft % this->chunkSize != 0*/) {
//                        temp = this->reads[i].firstFragment + (int)(fromLeft)/this->chunkSize + 1;
//                        tempLine = line.substr(this->reads[i].d);
//                        steps = (int) (tempLine.size()+this->chunkSize-20)/this->chunkSize;
//                        remain = (int) tempLine.size() - steps*this->chunkSize;
//                        if (temp > steps) {
//                            temp = steps;
//                        }
//                        this->depth->at(indexNumber).leftRead.push_back(onInterval());
//                        this->depth->at(indexNumber).leftRead.back().index = this->reads[i].index;
//                        this->depth->at(indexNumber).leftRead.back().fragment = temp;
//                        this->depth->at(indexNumber).leftRead.back().distance = (int) fromLeft - (temp-reads[i].firstFragment-1)*chunkSize;
//                        chunkACGT = "";
//                        for (int iter = temp; iter <= steps; iter++) {
//                            if (iter == steps) {
//                                chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize+remain);
//                            } else {
//                                chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize);
//                            }
//                        }
//                        chunkACGT = this->complement(chunkACGT);
//                        reverse(chunkACGT.begin(),chunkACGT.end());
//                        this->depth->at(indexNumber).leftSeq.push_back(string(chunkACGT));
//                    }

//                    // Reverse Right
//                    if (fromRight < this->chunkSize*(this->reads[i].firstFragment) /*&& fromRight % this->chunkSize != 0*/) {
//                        temp = this->reads[i].firstFragment - (int)(fromRight)/this->chunkSize;
//                        this->depth->at(indexNumber).rightRead.push_back(onInterval());
//                        this->depth->at(indexNumber).rightRead.back().index = this->reads[i].index;
//                        this->depth->at(indexNumber).rightRead.back().fragment = temp;
//                        this->depth->at(indexNumber).rightRead.back().distance = (int) fromRight - (reads[i].firstFragment-temp)*chunkSize;
//                        tempLine = line.substr(this->reads[i].d);
//                        chunkACGT = "";
//                        chunkACGT += line.substr(0, this->reads[i].d);
//                        for (int iter = 1; iter <= temp; iter++) {
//                            chunkACGT += tempLine.substr((iter-1)*this->chunkSize,this->chunkSize);
//                        }
//                        chunkACGT = this->complement(chunkACGT);
//                        reverse(chunkACGT.begin(),chunkACGT.end());
//                        this->depth->at(indexNumber).rightSeq.push_back(string(chunkACGT));
//                    }
//                }
//            }
//        }



//        lineCounter++;
//    }

//    return;
//}


//void Searching::searchingWithNoCutting(int t_id) {
//    int gapPen = -8, misPen = -6, matchPen = 2;
//    double editDistance = 0.3;
//    string name = "u_reads_d=0_p=" + convertNumToStr(t_id) + ".fq";
//    ifstream ifstr(name.c_str());
//    if (!ifstr.is_open()) {
//        string name = "CAN NOT OPEN THE FILE " + convertNumToStr(t_id);
//        perror(name.c_str());
//        return;
//    }
//    string line;
//    string segment;
//    long long lineCounter = 0;
//    long long index = 0;
//    int remain = 0;
//    int steps = 0;
//    int threshold = 0;

//    string genomeString;
//    string readString;

//    bool isInExon;
//    long long indexNumber;
//    long long fromRight;
//    long long fromLeft;

//    while (getline(ifstr, line)) {

//        if (lineCounter % 4 == 0) {

//            stringstream test(line);
//            getline(test, segment, '_');
//            getline(test, segment, '_');
//            index = stoll(segment.c_str())-1;

//        } else if (lineCounter % 4 == 1) {

//            line = line.substr(this->reads[index].d);
//            steps = (int) (line.size()+this->chunkSize-20)/this->chunkSize;
//            remain = (int) line.size() - steps*this->chunkSize;

//            if (this->reads[index].flag == 0 || this->reads[index].flag == 2) {

//                if (this->reads[index].lastFragment < steps) {
//                    for (int i = steps-this->reads[index].lastFragment; i >= 1 ; i--) {

//                        tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[index].lastPosition+i*this->chunkSize);
//                        if (isInExon == true) {
//                            if(i != steps-this->reads[index].lastFragment && fromRight>this->chunkSize+3) {
//                                readString = line.substr((this->reads[index].lastFragment+i-1)*this->chunkSize,this->chunkSize);
//                                genomeString = this->getGenome(this->reads[index].lastPosition+i*this->chunkSize-1,this->chunkSize);

//                                threshold = this->matchScore*this->chunkSize + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].lastFragment += i;
//                                    this->reads[index].lastPosition += i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;

//                            } else if (fromRight>this->chunkSize+remain+3) {
//                                readString = line.substr((this->reads[index].lastFragment+i-1)*this->chunkSize,this->chunkSize+remain);
//                                genomeString = this->getGenome(this->reads[index].lastPosition+i*this->chunkSize-1,this->chunkSize+remain);

//                                threshold = this->matchScore*(this->chunkSize+remain) + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].lastFragment += i;
//                                    this->reads[index].lastPosition += i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;
//                            }
//                        }
//                    }
//                }
//                if (1 < this->reads[index].firstFragment) {
//                    for (int i = this->reads[index].firstFragment-1; i >= 1; i--) {

//                        tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[index].firstPosition-i*this->chunkSize);
//                        if (isInExon == true) {
//                            if(fromLeft>this->chunkSize+3) {
//                                readString = line.substr((this->reads[index].firstFragment-i-1)*this->chunkSize,this->chunkSize);
//                                genomeString = this->getGenome(this->reads[index].firstPosition-i*this->chunkSize-1,this->chunkSize);

//                                threshold = this->matchScore*this->chunkSize + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].firstFragment -= i;
//                                    this->reads[index].firstPosition -= i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;
//                            }
//                        }
//                    }
//                }
//            }



//            if (this->reads[index].flag == 16 || this->reads[index].flag == 18) {

//                if (this->reads[index].lastFragment < steps) {
//                    for (int i = steps-this->reads[index].lastFragment; i >= 1 ; i--) {

//                        tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[index].lastPosition-i*this->chunkSize);
//                        if (isInExon == true) {
//                            if(i != steps-this->reads[index].lastFragment && fromLeft>this->chunkSize+3) {
//                                readString = line.substr((this->reads[index].lastFragment+i-1)*this->chunkSize,this->chunkSize);
//                                genomeString = this->getGenome(this->reads[index].lastPosition-i*this->chunkSize-1,this->chunkSize);
//                                genomeString = this->complement(genomeString);
//                                reverse(genomeString.begin(),genomeString.end());

//                                threshold = this->matchScore*this->chunkSize + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].lastFragment += i;
//                                    this->reads[index].lastPosition -= i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;

//                            } else if (fromLeft>this->chunkSize+remain+3) {
//                                readString = line.substr((this->reads[index].lastFragment+i-1)*this->chunkSize,this->chunkSize+remain);
//                                genomeString = this->getGenome(this->reads[index].lastPosition-i*this->chunkSize-1,this->chunkSize+remain);
//                                genomeString = this->complement(genomeString);
//                                reverse(genomeString.begin(),genomeString.end());

//                                threshold = this->matchScore*(this->chunkSize+remain) + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].lastFragment += i;
//                                    this->reads[index].lastPosition -= i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;
//                            }
//                        }
//                    }
//                }
//                if (1 < this->reads[index].firstFragment) {
//                    for (int i = this->reads[index].firstFragment-1; i >= 1; i--) {

//                        tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(this->reads[index].firstPosition+i*this->chunkSize);
//                        if (isInExon == true) {
//                            if(fromRight>this->chunkSize+3) {
//                                readString = line.substr((this->reads[index].firstFragment-i-1)*this->chunkSize,this->chunkSize);
//                                genomeString = this->getGenome(this->reads[index].firstPosition+i*this->chunkSize-1,this->chunkSize);
//                                genomeString = this->complement(genomeString);
//                                reverse(genomeString.begin(),genomeString.end());

//                                threshold = this->matchScore*this->chunkSize + this->indel*(this->indelScore-this->matchScore) + this->mismatch*(this->mismatchScore-this->matchScore);
//                                int indelShift = line.size()*0.1;
//                                LocalAligner *localAligner = new LocalAligner(genomeString, readString, gapPen, misPen, matchPen, indelShift, editDistance*line.size());
//                                localAligner->process();
//                                localAligner->backtrack();
//                                if(localAligner->mScore > threshold) {
//                                    this->reads[index].firstFragment -= i;
//                                    this->reads[index].firstPosition += i*this->chunkSize;
//                                    //delete localAligner;
//                                    break;
//                                }
//                                //delete localAligner;

//                            }
//                        }
//                    }
//                }

//            }

//        }
//        ++lineCounter;
//    }

//    ifstr.close();
//    return;
//}


//void Searching::concatExons() {

//    for (vector<interval>::iterator current = this->depth->begin(); current != this->depth->end()-1;) {

//        vector<interval>::iterator next = current + 1;
//        if (current->end >= next->start) {
//            long long end = max(current->end, next->end);
//            current->end = end;
//            this->depth->erase(next);
//        } else {
//            ++current;
//        }

//    }
//    return;
//}


//void Searching::breakReadsForThreads(int thread_num) {
//    ifstream ifstr("u_reads_d=0.fq");
//    string line;
//    int lineCounter = 0;
//    long long numberOfReadsInEachFile = this->numberOfReads / (thread_num);
//    ofstream * ofs = new ofstream[thread_num];
//    for (int i = 0; i < thread_num; i++){
//        string name = "u_reads_d=0_p=" + convertNumToStr(i) + ".fq";
//        ofs[i].open(name.c_str());
//    }

//    while (getline(ifstr, line)) {

//        if (numberOfReadsInEachFile == 0)
//            ofs[0]<<line<<endl;
//        else if	(lineCounter/(4*numberOfReadsInEachFile) < (thread_num-1))
//            ofs[lineCounter/(4*numberOfReadsInEachFile)]<<line<<endl;
//        else
//            ofs[thread_num-1]<<line<<endl;

//        ++lineCounter;
//    }

//    for (int i = 0 ; i < thread_num; i++)
//        ofs[i].close();
//    ifstr.close();
//    delete[] ofs;
//    return;
//}


//tuple<bool, long long, long long, long long> Searching::calculateIntronOrExonAndDistance(long long position) {
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


//string Searching::globalAlign(string s1, string s2) {
//    int up = 0, left = 1, up_left = 2;

//    int **score_mat;
//    score_mat = new int*[s1.length()+1];
//    for (int i = 0; i < s1.length()+1; ++i) {
//        score_mat[i] = new int[s2.length()+1];
//    }
//    for (int i = 0; i < s1.length()+1; ++i) {
//        for (int j = 0; j < s2.length()+1; ++j) {
//            score_mat[i][j] = 0;
//        }
//    }
//    int **path_mat;
//    path_mat = new int*[s1.length()+1];
//    for (int i = 0; i < s1.length()+1; ++i) {
//        path_mat[i] = new int[s2.length()+1];
//    }
//    for (int i = 0; i < s1.length()+1; ++i) {
//        for (int j = 0; j < s2.length()+1; ++j) {
//            path_mat[i][j] = 0;
//        }
//    }

//    string aligned1, aligned2;

//    for (int i = 0; i < s1.length()+1; i++) {
//        for (int j = 0; j < s2.length()+1; j++) {
//            if (i == 0) {
//                path_mat[i][j] = left;
//            }
//            if (j == 0) {
//                path_mat[i][j] = up;
//            }
//            if (i != 0 && j != 0) {
//                int v_up_left = score_mat[i-1][j-1] + score(s1[i-1], s2[j-1]);
//                int v_up = score_mat[i-1][j] + score(s1[i-1], '-');
//                int v_left = score_mat[i][j-1] + score('-', s2[j-1]);
//                int max_v = max(max(v_up_left, v_up), v_left);
//                if (v_up_left > v_left && v_up_left > v_up) {
//                    path_mat[i][j] = up_left;
//                } else if (v_left > v_up) {
//                    path_mat[i][j] = left;
//                } else {
//                    path_mat[i][j] = up;
//                }
//                score_mat[i][j] = max_v;
//            }
//        }
//    }

//    int i = (int) s1.length();
//    int j = (int) s2.length();
//    while (i > 0 || j > 0) {
//        if (path_mat[i][j] == up_left) {
//            aligned1 = s1[i - 1] + aligned1;
//            aligned2 = s2[j - 1] + aligned2;
//            i -= 1;
//            j -= 1;
//        } else if (path_mat[i][j] == up) {
//            aligned1 = s1[i - 1] + aligned1;
//            aligned2 = '-' + aligned2;
//            i -= 1;
//        } else {
//            aligned1 = '-' + aligned1;
//            aligned2 = s2[j - 1] + aligned2;
//            j -= 1;
//        }
//    }

//    string cigar = "";
//    for (int i=0; i < aligned1.length(); i++) {
//        if (aligned1[i] == aligned2[i]) {
//            cigar += 'M';
//        }
//        else if (aligned1[i] == '-'){
//            cigar += 'd';
//        }
//        else if (aligned2[i] == '-'){
//            cigar += 'i';
//        }
//        else{
//            cigar += 'm';
//        }
//    }
//    cigar = this->changeCigar(cigar);

//    for (int i = 0; i < s1.length()+1; ++i) {
//        delete [] score_mat[i];
//    }
//    delete [] score_mat;
//    for (int i = 0; i < s1.length()+1; ++i) {
//        delete [] path_mat[i];
//    }
//    delete [] path_mat;

//    return cigar;
//}


//int Searching::score(char a, char b) {
//    if (a == b) {
//        return 1;
//    } else {
//        return 0;
//    }
//}


//string Searching::getGenome(long long first, int length) {
//    struct stat sb;
//    off_t len;
//    char *p;
//    int fd;
//    string name = this->genomeFileName;
//    ifstream ifstr(name.c_str());
//    string genome;
//    string firstLine;
//    getline(ifstr, firstLine);
//    int firstLineCharCount = (int)firstLine.size();
//    ifstr.close();
//    first += firstLineCharCount+1 + first/this->characterPerRow;
//    fd = open(name.c_str(), O_RDONLY);
//    if (fd == -1)
//        perror ("open");

//    if (fstat (fd, &sb) == -1)
//        perror ("fstat");

//    p = (char *)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

//    if (close (fd) == -1)
//        perror ("close");

//    string temp = "";
//    for (len = first; temp.size() < length; ) {
//        if (p[len] != '\n')
//            temp += p[len];
//        ++len;
//    }

//    if (munmap (p, sb.st_size) == -1)
//        perror ("munmap");

//    //close(fd);
//    return temp;
//}


//string Searching::complement(string read){
//    string temp = "";
//    for (int i=0;i<read.size();i++){
//        switch	(read[i]){
//            case 'A':
//                temp+='T';
//                break;
//            case 'T':
//                temp+='A';
//                break;
//            case 'C':
//                temp+='G';
//                break;
//            case 'G':
//                temp+='C';
//                break;
//            case 'a':
//                temp+='t';
//                break;
//            case 't':
//                temp+='a';
//                break;
//            case 'c':
//                temp+='g';
//                break;
//            case 'g':
//                temp+='c';
//                break;
//            default:
//                temp += 'N';
//        }
//    }
//    return temp;
//}
