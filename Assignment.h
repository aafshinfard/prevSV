#ifndef ASSIGNMENT_H
#define ASSIGNMENT_H

#include <time.h>
#include <string.h>
#include <ios>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "Path.h"
#include "LocalAligner.h"
#include <time.h>
#include <thread>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <stdio.h>
#include <tgmath.h>
#include <unistd.h>
#include <map>
#include "Aligning.h"
#include "Header.h"
#include "Tools.h"

#define GENE_LEN	    3094257500
#define NOISE_THREADS    	10
#define BOWTIE1         1
#define BOWTIE2         2
#define Pacbio  false

using namespace std;

struct Read_Vector{
    long size;
    double index;
};

struct thread_data
{
    int thread_num;
    int length;
    int numberOfTables;
};

struct Read {
    bool Histo;
    bool Checked;
    long long readName;
    long index;
    string mainHeader;
    vector<string> references;
    vector<long long> loc;
    vector<int> seqSizes;
    vector<int> print_M;
    vector<int> print_I;
    vector<int> flags;
    vector<int> directions;
    vector<double> scores;
    vector<double> T_M_vector;
    vector<double> T_I_vector;
};

struct Tiny_Read{
    long long readName;
    bool Passed;
    long index;
    int V_Length;
};

struct littleStruct{
    int flag;
    double score;
    long long position;
    string ref;
};

class Assignment {

private:
    vector <long long> indel_stat;
    std::map<string, long long> multiFasta;

    iRead* reads;
    vector<interval>* iDepth;
    vector<gene>* genes;
    vector<iReadNext>* readsNext;

    long long frag;
    int minScore;
    long long filledAll;
    long long counter = 0;

    double consecutiveThreshold;
    double K_factor;
    double IndelShift;
    int option;
    int d;
    int V;
    int error;
    int GAPPEN,MISPEN,MATCH;
    int L_MAX;
    int numberOfTables, numberOfThreads, myDepth;
    int myNumberOfThread;
    int readLength, partLength;
    int errorInIndex;
    int aligner;
    int N_bowtie,L_bowtie;
    int step;
    int characterPerRow;
    string outputDir;

    int shift;
    string firstOutputName;

    bool local_Flag=false;

    const std::string space = " ";
    const std::string enter = "\n";
    std::string alignerAddress, fastQAddress, genomeName;
    std::string indexAddress;
    std::string line2;

    // I/O pointers
    ofstream* myOSams;
    ifstream* mySams;
    ifstream* myFasts;
    ofstream* myOFasts;
    ofstream* O_Strings;
    thread* Threads;

    // definition of all the vectors
    vector<Read> allReads;
    vector<Tiny_Read> tinyReads;
    vector<Read> threadReads;
    Tiny_Read sampleRead;
    vector< vector<Read> > threadVectors;
    vector<long> errorCollector;


    // Prototype of every necessary method
    void Join2pathes(Path &, Path &);
    long long Calculate(long long, int);
    int checkNumOfTables(int, int);

    int flagDetection(long long);

    int numeralCounter(int);
    void emptyInterval_Detector(vector<bool> &, vector<long> *, long long, int);
    bool histo_Permission();
    double findInitMax();

    string indexMaker(int);
    string complement(char);

    string string2mD(string, string);
    std::string accessGenome(long long, int, long long, long long);


    string finalString(vector<string> &, vector<string> &, long long);
    string L_MAX_Gate(int, int, string);
    int get_mD_Length(string, int &, int &);
    void Save_Properties(vector<string> &,vector< vector<int> > &,vector< vector<string> > & ,vector< vector<long long> > &,vector< vector<string> > &,int, int);
    void gen_Paths(int, int, vector< vector <long long> > &, vector< vector <int> > &, vector< vector<string> > &, vector <vector <string> > &, vector<long long> &, vector<Path> &, vector <long long> &);
    void Filters(vector<Path> &, vector<Path> &, vector<Path> &);

    void fill_The_Gaps(vector<Path> &, int &, vector<long> &, string &, double &, int &, double &, int &, int &, bool &, vector<string> &, int &, int &);

    void Analysis(int, int);
    void readFasta(string,std::map<string ,long long> &);

    void align_Command(int,string, string, int, int, int, string, int, int, int, int);

    std::string NumToStr(long long);
    long long NextChar(std::string, long long, string);
    vector<string> split(const string &, const string &);
    void Split_Header(vector<string> &, string &,string &, long long &, long long &, int &, int &, bool, bool, bool, bool, bool);
    void Split_Token(vector<string> &, string &, string &, string &, string &);
    string mD2String(string);
    string mD2Cigar(string);
    string Glue_cigars(string, string);
    string currentPath();
    string slice(int, int, string);
    string string2mD(string);
    bool forward_Detection(int);
    bool reverse_Detection(int);
    void fill_The_Gaps(vector<Path> &, int &, vector<long> &, string &, double &, int &, double &, int &, int &, bool &, vector<string> &, double &, double &);
    void merge_Method(vector<int> &, int , int , vector<Path> &, vector<bool> &, vector<string> &, int &, int &);
    void reform_Header(vector<string> &, string &, long long &, long long &, int, bool, bool, bool, bool, int i);
    void main_Header(string &, string &);
    void system_sort(string, string);
    void system_remove(string);
    void remain_Local(int, Path &, int &, int &, int &, string &, string);

    tuple<bool,long long, long long, long long> calculateIntronOrExonAndDistance(long long);
    void updateDepth(long long, long long, bool, int);
    void addIntervalToDepth(long long, long long, long long);
    void addIntervalToDepth(long long, long long);
    void concatExons(long long);
    void AddToJenes(gene&);
    void mergeExonsInGenesAndTempGene(long long, long long, long long);
    string convertNumToStr(long long);

public:
    Assignment(iRead*, vector<interval>*, vector<gene>*, vector<iReadNext>*, int, int, int, string, int, int, int, int, string, string, int, string);
    void runForOnGenes();
    void createNewRead(int, string);
    ~Assignment();

};

#endif // ASSIGNMENT_H
