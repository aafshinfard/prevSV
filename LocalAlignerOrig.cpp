#include "LocalAlignerOrig.h"
using namespace std;
LocalAlignerOrig::LocalAlignerOrig(std::string seqA, std::string seqB, int gapPen, int misPen, int matchPen, double shiftt, int scoreInit = 0)
{
    this->MATCH = matchPen;
    this->GAP = gapPen;
    this->counter = 0;
    this->MISMATCH = misPen;
    this->shift = (int)(shiftt*seqA.size()+1);
    for (int i=0;i < seqA.size(); i++)
        seqA[i] = toupper(seqA[i]);
    for (int i=0;i < seqB.size(); i++)
        seqB[i] = toupper(seqB[i]);
    this->mSeqA = seqA;
    this->mSeqB = seqB;
    this->mD = new int*[this->mSeqA.size() + 1];
    for (int i = 0; i <= this->mSeqA.size(); i++)
        this->mD[i] = new int[this->mSeqB.size() + 1];
    this->mD[0][0] = scoreInit;
}

LocalAlignerOrig::~LocalAlignerOrig()
{
    for (int i = 0; i <= this->mSeqA.size(); i++)
        delete[] this->mD[i];
    delete this->mD;
}

void LocalAlignerOrig::process()
{
    int scoreDiag = 0;
    int scoreLeft = 0;
    int scoreUp = 0;
    int maximH = 0;
    int maximV = 0;
    int maxim = 0;
    for (int i = 1; i <= this->mSeqA.size(); i++)
    {
        for (int j = 1; j <= this->mSeqB.size(); j++)
        {
            scoreDiag = mD[i-1][j-1] + weight(i, j);
            for (int k = 1; k <= this->mSeqA.size(); k++)
            {
                mD[k][0] = std::max(0,mD[0][0]+k*GAP);
            }
            for (int k = 1; k <= this->mSeqB.size(); k++)
            {
                mD[0][k] = std::max(0,mD[0][0]+k*GAP);
            }
            maximV=0;
            for(int k = 1; k <= std::min(shift,i-1); k++)
            {
                scoreUp = mD[i-k][j] + k*GAP;
                maximV = std::max(maximV,scoreUp);
            }
            maximH=0;
            for(int k = 1; k <= std::min(shift,j-1); k++)
            {
                scoreLeft = mD[i][j-k] + k*GAP;
                maximH = std::max(maximH,scoreLeft);
            }
            maxim = std::max(maximH,maximV);
            mD[i][j] = std::max(std::max(0,scoreDiag),maxim);
        }
    }
}
void LocalAlignerOrig::Print(){

    for(int i=0;i<=this->mSeqA.size();i++){
        for(int j=0;j<=this->mSeqB.size();j++){
            cout << j << ": " << mD[i][j]<< "\t";
        }
        cout << endl;
    }
}

 void LocalAlignerOrig::backtrack(){
     size_t i = this->mSeqA.size();
     int max = -1000 * MATCH , index = 0, flag = 0;
     for (int k = 0; k <= this->mSeqB.size(); k++)
     {
         if(mD[i][k] > max)
         {
             max = mD[i][k];
             index = k;
         }
     }
     size_t j = index;
     this->mScore = mD[i][j];
     std::ostringstream os;

     while (i > 0 && j > 0) {

         if (mD[i][j] == mD[i-1][j] + GAP) {
             mAlignmentSeqA += this->mSeqA[i-1];
             mAlignmentSeqB += "-";
             totalGapfa ++;
             i--;
             continue;
         }
         else if (mD[i][j] == mD[i-1][j-1] + weight(i, j)) {  //MATCH-MISMATCH
             mAlignmentSeqA += this->mSeqA[i-1];
             mAlignmentSeqB += this->mSeqB[j-1];
             if(weight(i,j)==MISMATCH)
                 totalMismatch ++;
             i--;
             j--;
             continue;
         }
         else
         {
             mAlignmentSeqA += "-";
             totalGapr ++;
             mAlignmentSeqB += this->mSeqB[j-1];
             j--;
             continue;
         }
     }
     std::reverse(this->mAlignmentSeqA.begin(), mAlignmentSeqA.end());
     std::reverse(this->mAlignmentSeqB.begin(), mAlignmentSeqB.end());
     //produceCigar();
 }

 // void LocalAlignerOrig::produceCigar(){
 //     std::string temp;
 //     temp = this->mAlignmentSeqA;
 //     this->mAlignmentSeqA = this->mAlignmentSeqB;
 //     this->mAlignmentSeqB = temp;
 //     int counter = 0, realCounter = 0;
 //     cigar = "";
 //     realCigar = "";
 //     int gapCount= 0, Icounter = 0, Dcounter = 0;
 //     std::ostringstream os;
 //     for (int i = 0 ; i < this->mAlignmentSeqA.size(); i++) {
 //         if (mAlignmentSeqA[i] == mAlignmentSeqB[i]){
 //             gapCount = 0;
 //             counter++;
 //             if (Icounter) {
 //                 os.str("");
 //                 os<<Icounter;
 //                 realCigar += os.str();
 //                 realCigar += "I";
 //             }
 //             if (Dcounter) {
 //                 os.str("");
 //                 os<<Dcounter;
 //                 realCigar += os.str();
 //                 realCigar += "D";
 //             }
 //             Icounter = 0;
 //             Dcounter = 0;
 //             realCounter++;
 //         }
 //         else if (mAlignmentSeqA[i] != mAlignmentSeqB[i] && mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] != '-'){
 //             gapCount = 0;
 //             if (Icounter) {
 //                 os.str("");
 //                 os<<Icounter;
 //                 realCigar += os.str();
 //                 realCigar += "I";
 //             }
 //             if (Dcounter) {
 //                 os.str("");
 //                 os<<Dcounter;
 //                 realCigar += os.str();
 //                 realCigar += "D";
 //             }
 //             if (counter) {
 //                 os.str("");
 //                 os<<counter;
 //                 cigar += os.str();
 //             }
 //             Icounter = 0;
 //             Dcounter = 0;
 //             counter = 0;
 //             cigar += mAlignmentSeqA[i];
 //             realCounter++;
 //         }
 //         else if (mAlignmentSeqA[i] == '-' && mAlignmentSeqB[i] != '-') {
 //             gapCount ++;
 //             Icounter ++;
 //             if (counter) {
 //                 os.str("");
 //                 os<<counter;
 //                 cigar += os.str();
 //             }
 //             if (gapCount == 1) cigar += "^";
 //             cigar += mAlignmentSeqB[i];
 //             counter = 0;

 //             if (realCounter) {
 //                 os.str("");
 //                 os<<realCounter;
 //                 realCigar += os.str();
 //                 realCigar += "M";
 //             }
 //             realCounter = 0;

 //             if (Dcounter) {
 //                 os.str("");
 //                 os<<Dcounter;
 //                 realCigar += os.str();
 //                 realCigar += "D";
 //             }
 //             Dcounter = 0;
 //         }
 //         else if (mAlignmentSeqA[i] != '-' && mAlignmentSeqB[i] == '-'){
 //             gapCount = 0;
 //             Dcounter ++;
 //             if (realCounter) {
 //                 os.str("");
 //                 os<<realCounter;
 //                 realCigar += os.str();
 //                 realCigar += "M";
 //             }
 //             realCounter = 0;
 //             if (Icounter) {
 //                 os.str("");
 //                 os<<Icounter;
 //                 realCigar += os.str();
 //                 realCigar += "I";
 //             }
 //             Icounter = 0;
 //         }

 //     }
 //     if (counter) {
 //         os.str("");
 //         os<<counter;
 //         cigar += os.str();
 //     }

 //     if (realCounter) {
 //         os.str("");
 //         os<<realCounter;
 //         realCigar += os.str();
 //         realCigar += "M";
 //     }
 //     if (Dcounter) {
 //                 os.str("");
 //                 os<<Dcounter;
 //                 realCigar += os.str();
 //                 realCigar += "D";
 //     }
 //     if (Icounter) {
 //                 os.str("");
 //                 os<<Icounter;
 //                 realCigar += os.str();
 //                 realCigar += "I";
 //     }
 // }


int LocalAlignerOrig::weight(size_t i, size_t j) {
    ++counter;
    if (this->mSeqA[i - 1] == this->mSeqB[j - 1])
        return MATCH;
    else
        return MISMATCH;
}
