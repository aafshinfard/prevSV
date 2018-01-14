 //
 //  SVAnchor.cpp
 //  SV-Detection
 //
 //  Created by Ameerhosein Afshinfard - BRLab - Sharif Uiversity of Technology
 //  Copyright (c) 2016 Ameerhosein Afshinfard. All rights reserved.
 //

#include "SVAnchor.h"
#include "unistd.h"

//// CHANGES NEEDED YET:
/// 1- Yek moroor koli rooye code va hazfe khandan nevashtan haye ezafe ! (1:asoon 2:sakht o toolani)
/// 2- Regional local alignment, for incluster remained regions
///
// ================================|
// ================================|
//> Settings:                      |
// ================================|
// ==== Normal Settings:
#define MAX_BUFFER_SIZE 100     // each thread will read MAX_BUFFER_SIZE reads in each trial
bool BGLR = true;               // [true]: Left and Right Binary Genomes
                                // [false]: only 1 BinaryGenome (less accurate and slower But less Ram usage)
int anchoringIterations = 2;
int LocAlth = 130;
bool logTxt = false;
bool pairedEndedReads = false;  // if true, then reads must be in two files "readName"_1.fq and "readName"_2.fq
bool runPhase2 = false;

// ===========================
// ==== for multiple Files:
int readsFileCount = 2;

int multiFilesTogether = 50; // for readsFilesCount>1

// ===========================
// ==== for synteny Detection:
// No:
bool runInClusterExtension = true;
bool runExtensionLoop = true;
bool changeBinGenomeOnAnchoring = false;
// Yes:
//bool runInClusterExtension = false;
//bool runExtensionLoop = false;
//bool changeBinGenomeOnAnchoring = true/false;

// ===========================
// ==== old Versions :
int version = 2;// new multi threads
bool extensionIn1Thread = true; // works for version=1 only


// ================================|
//> End of Settings:               |
// ================================|

//> Global counters
bool onlyBG = false;
long long BGRchanges = 0, BGLchanges = 0, BGRchecks = 0;
long long changed=0, skipped=0;
long long cntCOM = 0, cntBG=0;

//> Global Variables
vector<int> shiftPos;
vector<int> prevShiftPos;

//> Binary Genomes
vector<bool> bGRight;
vector<bool> bGLeft;

SVAnchor::SVAnchor(string readName, string genomeName, string indexName, int chunkSize, int numberOfThread, int dashV, int numGap, int anchoringShift, int assignmentShift, bool isSlideNotShift, int readLength, string outputDir,long long numberOfReads1) {
    this->numberOfReads = numberOfReads1;
    this->genomeFileName = genomeName;
    //this->lenVar = 5;
    this->readName = readName;
    this->genomeName = genomeName;
    this->indexName = indexName;
    this->chunkSize = chunkSize;
    this->numberOfThread = numberOfThread;
    this->dashV = dashV;
    this->numGap = numGap;
    //this->anchoringShift = anchoringShift;
    this->readLength = readLength;
    this->outputDir = outputDir;
    this->extensionThreadIntervals = new int[numberOfThread];
    this->indelShift = 0.4;
    this->editDistance = 2;
    this->genomeLength = getGenomeLength();
    cout<<"\n Genome Size is: "<<this->genomeLength<<"\n";
    bGRight.resize(this->genomeLength);
    bGLeft.resize(this->genomeLength);
    for(int i=0 ; i < numberOfThread ; i++ ){
        shiftPos.push_back(0);
        prevShiftPos.push_back(0);
    }
}

SVAnchor::~SVAnchor() {
    delete [] this->reads;
    pthread_mutex_destroy(&lock);
    vector<interval>().swap(this->depth);
    vector<gene>().swap(this->genes);
    vector<iReadNext>().swap(this->readsNext);
}

void SVAnchor::writeScores(){
    string name = outputDir+"LocAliScores.txt";
    ofstream ofstr(name.c_str());
    for(vector<int>::iterator it= AlignmentScores.begin() ; it!=AlignmentScores.end() ; ++it )
        ofstr<< *it <<endl;
    ofstr.close();
    return;
}

void SVAnchor::automation(bool runAnchoring) {

    auto begin = std::chrono::high_resolution_clock::now();
    chdir(MAINADDRESS);

    if(logTxt){
        string command = "cat /dev/null > log.txt";
        system(command.c_str());
        command = "mkdir "+outputDir;
        system(command.c_str());
        freopen("log.txt", "w", stdout);
    }
    this->localAlThreshold = LocAlth;

    pFS(get_current_time());pFS(" Starting SV-Anchoring run (v1.4.1.2)");pFS("\n");
    pFS("\t\t Alignment Treshold: ");
    pFS(convertNumToStr((this->localAlThreshold))+"\n");
    pFS("\n===============================================\n");


    this->reads = new iRead[this->numberOfReads];
    this->characterPerRow = findCharacterPerRow(this->genomeName);

    /*=======================================================================*/
    /*============================   ANCHORING   ============================*/
    /*=======================================================================*/
    if (runAnchoring) {
        pFS(get_current_time());pFS(" SAMPLING file starts");pFS("\n");

        Sampling *sampling = new Sampling(this->reads, &this->depth, this->characterPerRow, this->genomeName, 5, this->outputDir);
        sampling->prepareSampling(this->readName , pairedEndedReads, readsFileCount );
        delete sampling;
        pFS(" Hi there ");pFS("\n");
        pFS(get_current_time());pFS(" Total number of reads is ");pFS(convertNumToStr(numberOfReads));pFS("\n");
        //pFS(get_current_time());pFS(" ANCHORING step:");pFS("\n");

        long long counter1 = 0;
        //int i = 0;
        //int shift = this->anchoringShift;
        //for (i = 0; i < floor(this->chunkSize/shift); i++) {

        pFS(get_current_time());pFS(" Aligning started....\n");

        if( extensionIn1Thread && version==1 )
            this->numberOfThread = 1;
        //
        pFS(get_current_time());
        pFS(" preparing files started....\n");//pFS(this->reads->);

        if(readsFileCount==1){

            // - // ************************************************************************************************************************************************** prepare
            cout<<"\n===============================================\n";
            pFS(" Anchoring Started.... ");
            cout<<"\n===============================================\n";
            int step = chunkSize / anchoringIterations ;
            for(int it = 0 ; it < anchoringIterations ; it++ ){
                cout<<"\n=====\n";
                pFS(" Anchoring  Loop #"+convertNumToStr(it+1));
                cout<<"\n=====\n";
                Anchoring *anchoring = new Anchoring(this->reads, &this->depth, this->indexName, "an", true, this->numGap, this->chunkSize, it*step, 5, this->dashV, this->numberOfThread, this->outputDir);
                anchoring->prepareFiles();
                //mkdirpFS(get_current_time());pFS(" files prepared");pFS("\n");

                // - // ************************************************************************************************************************************************** nextstep
                anchoring->nextStepAlignment();
                counter1 += anchoring->correspond;
                pFS(get_current_time());
                pFS("\n***********************************************\n this Loop Finished.... ");

                pFS(get_current_time());
                pFS("\n Number of anchoring in this step is ");pFS(convertNumToStr(anchoring->correspond));pFS("\n");
                delete anchoring;
            }
        }else{
            cout<<"\n===============================================\n";
            pFS(" Anchoring Started.... ");
            cout<<"\n===============================================\n";
            int i = 0;
            int step = chunkSize / anchoringIterations ;
            for(int it = 0 ; it < anchoringIterations ; it++ ){
                i = 0;
                cout<<"\n=====\n";
                pFS(" Anchoring  Loop #"+convertNumToStr(it+1));
                cout<<"\n=====\n";
                Anchoring *anchoring = new Anchoring(this->reads, &this->depth, this->indexName, "an", true, this->numGap, this->chunkSize, it*step, 5, this->dashV, this->numberOfThread, this->outputDir);
                do{
                    // use 3rd and 4th parameters to shift the anchoring, need some changes in function yet
                    // - // ************************************************************************************************************************************************** prepare
                    anchoring->prepareFiles(i*multiFilesTogether +1,(i+1)*multiFilesTogether,readsFileCount,0,0);
                    // - // ************************************************************************************************************************************************** nextstep
                    anchoring->nextStepAlignment(i*multiFilesTogether +1,(i+1)*multiFilesTogether,readsFileCount);

                    // use 3rd and 4th parameters to shift the anchoring, need some changes in function yet
                    i += 1;
                }while((i+1)*multiFilesTogether < (anchoring->getNumberFiles()/2) );
                this->writeReads( "beforExt_shift_"+(it*step) );
                counter1 += anchoring->correspond;
                pFS(get_current_time());
                pFS("\n***********************************************\n this Loop Finished.... ");

                pFS(get_current_time());
                pFS("\n Number of anchoring in this step is ");pFS(convertNumToStr(anchoring->correspond));pFS("\n");
                delete anchoring;
            }
        }

        //this->deleteAllJunk("an", shift*(i+1));

        //pFS(get_current_time());pFS(" Total anchoring is ");pFS(convertNumToStr(counter1));pFS("\n");
        // here we need to extend the reads
        this->writeReads("beforExt");
        if(changeBinGenomeOnAnchoring){
            cout<<"\n===============================================\n";
            pFS(" BinGenome modification using anchoring results...");
            cout<<"\n===============================================\n";
            changeBinGenomeByAnchoring();
        }
        if( version == 1 ){
            cout<<"\n===============================================\n";
            pFS(" Finding Intervals for out cluster Extension");
            cout<<"\n===============================================\n";
            this->keepReadsForExtension(); // now Intervals and read indices is stored
            this->breakReadsForExtension(true);
        }
        cout<<"\n===============================================\n";
        pFS(" Extension Starts...");
        cout<<"\n===============================================\n";
        pFS(convertNumToStr(this->numExtendables)) ;
        // MultiThread Extension:
        if(version==1){
            int rc;
            pthread_t threads[this->numberOfThread];
            pthread_attr_t attr;

            params_t params;
            pthread_mutex_init (&params.mutex , NULL);
            pthread_cond_init (&params.done, NULL);

            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

            for(int i=0; i < this->numberOfThread; i++ ){
                params.id = i;
                params.ptr = this;
                rc = pthread_create(&threads[i], NULL, &(SVAnchor::extensionStepWrapper), &params);
                if (rc)
                {
                    cout << "Error:unable to create thread," << rc << endl;
                    exit(-1);
                }else{
                    pthread_cond_wait (&params.done, &params.mutex);
                }
            }
            pthread_attr_destroy(&attr);
            void* status;
            for(int i=0; i < this->numberOfThread; i++ )
            {
                rc = pthread_join(threads[i], &status);
                if (rc)
                {
                    cout << "Error:unable to join," << rc << endl;
                    exit(-1);
                }
                cout << "Main: completed thread id :" << i ;
                cout << "  exiting with status :" << status << endl;
            }
            pthread_mutex_destroy (&params.mutex);
            pthread_cond_destroy (&params.done);
        }
        if(version == 2 && readsFileCount == 1){
            long long globalLineCounter;
            globalLineCounter = 0;
            int *ids ;
            string name =  this->readName;
            ifstream ifstr(name.c_str());
            ids = new int[numberOfThread];
            std::vector<std::thread> threads(numberOfThread);
            for (int i = 0; i < numberOfThread; i++) {
                ids[i] = i;
                threads[i] = std::thread(&SVAnchor::extensionStep2, this, &ifstr , &globalLineCounter , &ids[i], &readsFileCount);
            }
            // wait until all threads have finished
            for (int i = 0; i < numberOfThread; i++) {
                threads[i].join();
            }
            ifstr.close();
            delete ids;
        }
        if(version == 2 && readsFileCount > 1){
            long long globalLineCounter;
            globalLineCounter = 0;
            int *ids ;
            string name;
            for(int i = 0 ; i < readsFileCount ; i++){
                name =  this->readName +"_"+convertNumToStr(i)+".fq";
                ifstream ifstr(name.c_str());
                ids = new int[numberOfThread];
                std::vector<std::thread> threads(numberOfThread);
                for (int i = 0; i < numberOfThread; i++) {
                    ids[i] = i;
                    threads[i] = std::thread(&SVAnchor::extensionStep2, this, &ifstr , &globalLineCounter , &ids[i], &readsFileCount);
                }
                // wait until all threads have finished
                for (int i = 0; i < numberOfThread; i++) {
                    threads[i].join();
                }
                ifstr.close();
                delete ids;
            }
        }
        pFS("\n***********************************************\n SUCCESSFUL EXTENSION ");
        cout<<"\n Number of bases Remained after Extension 1 in BG(s) is :"<<endl;
        if(BGLR)
            cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
        else
            cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;
        pFS(get_current_time());
        cout<<"\n Extended by BG:"+convertNumToStr(cntBG)+" / by Alignment:"+convertNumToStr(cntCOM)+
              "\n***********************************************\n";
        this->writeReads("AfterExt");
        if(runInClusterExtension){
            cout<<"\n===============================================\n";
            cout<<" inClusterExtension() starts ... ";
            cout<<"\n===============================================\n";
            inClusterExtension(readsFileCount);// version will check inside
            pFS("\n***********************************************\n"
                " SUCCESSFUL inClusterExtension"
                "\n***********************************************\n" );
            cout<<"\n changed: "<<changed<<" Skipped: "<<skipped<<endl;
            this->writeReads("AfterInClusterExt");


            cout<<"\n Extended by BG:"+convertNumToStr(cntBG)+" / by Alignment:"+convertNumToStr(cntCOM)+
                  "\n***********************************************\n";
        }
        // //////// //////// ////// //////// //////// ////
        //          LOOP
        //        Extension
        // //////// //////// ////// //////// //////// ////
        if(version == 1){
            for(int h=0 ; h < (runExtensionLoop?3:0) /*scanCounts*/ ; h++)
            {
                cntCOM = 0, cntBG=0;
                onlyBG = true;
                this->keepReadsForExtension(); // now Intervals and read indices is stored
                this->breakReadsForExtension(true);
                cout<<"\n===============================================\n";
                cout<<" Rescan Extension "<<h+1<<" Starts...";
                cout<<"\n===============================================\n";
                pFS(convertNumToStr(this->numExtendables)) ;
                // MultiThread Extension:
                int rc2;
                pthread_t threads2[this->numberOfThread];
                pthread_attr_t attr2;

                params_t params2;
                pthread_mutex_init (&params2.mutex , NULL);
                pthread_cond_init (&params2.done, NULL);

                pthread_attr_init(&attr2);
                pthread_attr_setdetachstate(&attr2, PTHREAD_CREATE_JOINABLE);

                for(int i=0; i < this->numberOfThread; i++ ){
                    params2.id = i;
                    params2.ptr = this;
                    rc2 = pthread_create(&threads2[i], NULL, &(SVAnchor::extensionStepWrapper), &params2);
                    if (rc2)
                    {
                        cout << "Error:unable to create thread," << rc2 << endl;
                        exit(-1);
                    }else{
                        pthread_cond_wait (&params2.done, &params2.mutex);
                    }
                }
                pthread_attr_destroy(&attr2);
                void* status;
                for(int i=0; i < this->numberOfThread; i++ )
                {
                    rc2 = pthread_join(threads2[i], &status);
                    if (rc2)
                    {
                        cout << "Error:unable to join," << rc2 << endl;
                        exit(-1);
                    }
                    cout << "Main: completed thread id :" << i ;
                    cout << "  exiting with status :" << status << endl;
                }
                pthread_mutex_destroy (&params2.mutex);
                pthread_cond_destroy (&params2.done);
                pFS("\n***********************************************\n SUCCESSFUL EXTENSION ");
                cout<<"\n Number of bases Remained after Extension "<<h+2<<" in BG(s) is :"<<endl;
                if(BGLR)
                    cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
                else
                    cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;
                pFS(get_current_time());
                cout<<"\n Extended by BG:"+convertNumToStr(cntBG)+" / by Alignment:"+convertNumToStr(cntCOM)+
                      "\n***********************************************\n";

                //delete threads2;
                //delete params2;
                //delete attr2;
            }
        }
        if(version == 2){

            for(int h=0 ; h < (runExtensionLoop?3:0) /*scanCounts*/ ; h++)
            {
                cntCOM = 0, cntBG=0;
                onlyBG = true;
                this->keepReadsForExtension(); // now Intervals and read indices is stored
                this->breakReadsForExtension(true);
                cout<<"\n===============================================\n";
                cout<<" Rescan Extension "<<h+1<<" Starts...";
                cout<<"\n===============================================\n";
                pFS(convertNumToStr(this->numExtendables)) ;
                long long globalLineCounter;
                globalLineCounter = 0;
                int *ids ;
                string name =  this->readName;
                ifstream ifstr(name.c_str());
                ids = new int[numberOfThread];
                std::vector<std::thread> threads(numberOfThread);
                for (int i = 0; i < numberOfThread; i++) {
                    ids[i] = i;
                    threads[i] = std::thread(&SVAnchor::extensionStep2, this, &ifstr , &globalLineCounter , &ids[i], &readsFileCount);
                }
                // wait until all threads have finished
                for (int i = 0; i < numberOfThread; i++) {
                    threads[i].join();
                }
                ifstr.close();
                delete ids;
                pFS("\n***********************************************\n SUCCESSFUL reEXTENSION ");
                cout<<"\n Number of bases Remained after Extension "<<h+2<<" in BG(s) is :"<<endl;
                if(BGLR)
                    cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
                else
                    cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;
                pFS(get_current_time());
                cout<<"\n Extended by BG:"+convertNumToStr(cntBG)+" / by Alignment:"+convertNumToStr(cntCOM)+
                      "\n***********************************************\n";
                if(cntBG == 0 ){
                    cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
                    cout<<"\n No more Rescan will helps. so we quite loop of rescanning.";
                    cout<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
                    break;
                }
            }
        }
        if(runExtensionLoop)
            this->writeReads("AfterExtLoop");
        //// //////// //////// ////
        ////        End
        ////   Loop Extension
        //// //////// //////// ////

        cout<<"\n ==================\n start to write BinGenome and SingleUniques"<<endl;
        writeBinGenome();
        readSingleUniques();
        pruneSingleUniques();
        writeSingleUniques();

        //cout<<"start to write AfterExt"<<endl;

        writeScores();
        //this->writeDepth("an");
        pFS("-----------------------------------------------\n");

        // Signalling::Signalling(iRead* reads , string readName, string genomeName,
        // int chunkSize, int numberOfThread, int numGap, string outputDir,long long numberOfReads1)
        if(runPhase2){
            Signalling *signalling = new Signalling(this->reads,this->readName,this->genomeName,this->chunkSize,this->numberOfThread,this->numGap, this->outputDir, this->numberOfReads);
            signalling->automation();
            delete signalling;
        }
    }

    auto elapsed = std::chrono::high_resolution_clock::now() - begin;
    long long milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();

    int seconds = (int) (milliseconds / 1000) % 60 ;
    int minutes = (int) ((milliseconds / (1000*60)) % 60);
    int hours   = (int) ((milliseconds / (1000*60*60)) % 24);
    cout<<"\n # bases tried to changed: "<<BGRchanges<<" - and # bases tried to check: "<<BGRchecks<<endl;
    pFS("===============================================\n");
    pFS(get_current_time());pFS(" Run completed: ");
    pFS(convertNumToStr(hours)+":");
    pFS(convertNumToStr(minutes)+":");
    pFS(convertNumToStr(seconds)+ " elapsed\n");
    if(logTxt)
        freopen ("/dev/tty", "a", stdout);
    pFS(get_current_time());pFS(" Run completed: ");
    pFS(convertNumToStr(hours)+":");
    pFS(convertNumToStr(minutes)+":");
    pFS(convertNumToStr(seconds)+ " elapsed\n");

    //pthread_exit(NULL);
    return;
}


long long SVAnchor::basesRemained(bool isToRight){
    // Number of bases Remained in bin genomes (right OR left)
    long long bRemained = 0;
    //cout<<"\n Size of Genome"<<bGRight.size()<<endl;
    if(isToRight){
        for(long long i = 0 ; i < bGRight.size() ; i++ )
            if(bGRight[i]==false)
                bRemained++;
    }else{
        for(long long i = 0 ; i < bGRight.size() ; i++ )
            if(bGLeft[i]==false)
                bRemained++;
    }

    return bRemained;
}
inline bool SVAnchor::checkBinGenome(long long start , long long end, bool isOr, bool isToRight=true ){
    if( start > end){
        int temp = start;
        start = end;
        end = temp;
    }
    if( start < 0 || end > genomeLength-1 ){
        cout<<"\n Genome Ends Breaked.";
        return false;
    }
    if(!isOr) BGRchecks += abs(end - start)+1;
    if( isOr ){
        for(long long i = start; i < end+1 ; i++)
            if(bGRight[i] == false || bGLeft[i] == false )
                return false;

    }else if( isToRight ){
        for(long long i = start; i < end+1 ; i++)
            if(bGRight[i] == false)
                return false;
    }else{
        for(long long i = start; i < end+1 ; i++)
            if(bGLeft[i] == false)
                return false;
    }
    return true;
}

inline void SVAnchor::changeBinGenome(long long start , long long end, bool isToRight){

    if(start > end){
        long long temp = start;
        start = end;
        end = temp;
    }

    BGRchanges += abs(end - start)+1;
    if( isToRight ){
        bGRChangeMutex.lock();
        for(long long i = start; i < end+1 ; i++)
            bGRight[i] = true;
        bGRChangeMutex.unlock();
    }else{
        bGLChangeMutex.lock();
        for(long long i = start; i < end+1 ; i++)
            bGLeft[i] = true;
        bGLChangeMutex.unlock();
    }

    //pthread_mutex_unlock(&lock);// //////////////////////////////////////// UnLock

}
void SVAnchor::changeBinGenomeByAnchoring(){
    for (long long i = 0; i < this->numberOfReads; i++) {
        if( abs(reads[i].firstFragment - reads[i].lastFragment) > 2){
            changeBinGenome( min(reads[i].firstPosition,reads[i].lastPosition)+2*chunkSize , max(reads[i].firstPosition,reads[i].lastPosition)-chunkSize , false );
            changeBinGenome( min(reads[i].firstPosition,reads[i].lastPosition)+2*chunkSize , max(reads[i].firstPosition,reads[i].lastPosition)-chunkSize , true );
        }
        for( int j = 1 ; j < reads[i].firstFragClus2.size() ; j++ )
            if( abs(reads[i].firstFragClus2[j] - reads[i].lastFragClus2[j]) > 2){
                changeBinGenome( min(reads[i].firstPosClus2[j],reads[i].lastPosClus2[j])+2*chunkSize , max(reads[i].firstPosClus2[j],reads[i].lastPosClus2[j])-chunkSize  , false );
                changeBinGenome( min(reads[i].firstPosClus2[j],reads[i].lastPosClus2[j])+2*chunkSize , max(reads[i].firstPosClus2[j],reads[i].lastPosClus2[j])-chunkSize  , true );
            }
    }
}

inline string SVAnchor::revComplemACGT(string readSegment){
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
string SVAnchor::getGenome(long long first, int length) {
    struct stat sb;
    off_t len;
    char *p;
    int fd;
    string name = this->genomeFileName;// + ".fasta";

    getGenomeMutex.lock();

    ifstream ifstr(name.c_str());
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

    //close(fd);
    getGenomeMutex.unlock();
    return temp;
}
long long SVAnchor::getGenomeLength(){

    ifstream file((MAINADDRESS+genomeName).c_str(), ifstream::in | ifstream::binary);
    if(!file.is_open())
        return -1;

    string line1;
    getline(file,line1);
    int line1Size = line1.size();

    file.seekg(0, ios::end);
    long long fileSize = file.tellg();
    file.close();

    fileSize = fileSize - line1Size - 1 - 1 ;/*the last line has no \n*/
    long remainedLineCount = (fileSize)/81 ;
    long long GenomeSize = fileSize - remainedLineCount;
    return GenomeSize;
}

inline char SVAnchor::ifAlignable(long long startLoci, string readString, int flag, int threadNum, bool isToRight, int jumpsCount, bool onlyBGlocal, bool useBinGenome = true){
    int score = 0;
    // changing value of indelshift
    double indelShiftLocal = indelShift + (numGap*jumpsCount)/chunkSize;
    bool binGenome;
    if( useBinGenome ){
        binGenome = checkBinGenome(startLoci, startLoci + this->chunkSize, false, (BGLR ? isToRight : true));
        if(binGenome){
            cntBG++;
            return 1;
        }
        if( onlyBGlocal )
            return 0;
    }
    cntCOM++;
    LocalAligner *local;

    if( flag!=0 )
        readString = revComplemACGT(readString);

    int gapPen = -8, misPen = -5, matchPen = 10;
    if(startLoci - indelShiftLocal*readString.size() < 0 || startLoci + indelShiftLocal*readString.size()+this->chunkSize  > genomeLength-1 )
        return 0;
    string genomeString = this->getGenome(startLoci - indelShiftLocal*readString.size(), this->chunkSize + 2*readString.size()*indelShiftLocal );

    string seqConstantB = "";
    int readlen = readString.size();
    for(int i=0;i < (int)(indelShiftLocal*readlen);i++)
        seqConstantB += "A";
    readString += seqConstantB;

    local = new LocalAligner(readString, genomeString, gapPen, misPen, matchPen, (int)(indelShiftLocal*readlen), (int)(editDistance*readlen) );
    local->process();
    local->backtrack();
    if(startLoci == 7994 )
        local->Print();
    //cout//<<local->realCigar << endl;
    //cout<<local->cigar << endl;
    score = local->mScore;
    //AlignmentScores.push_back(score);
    char isAlignable = (score > localAlThreshold ? 2 : 0 );
    if(isAlignable && threadNum > -1 ){
        prevShiftPos.at(threadNum) = shiftPos.at(threadNum);
        shiftPos.at(threadNum) = local->totalGapr - local->totalGapfa;
    }
    delete local;
    return isAlignable;
}
inline int SVAnchor::ifAlignable(long long startLoci, string readString, int flag, bool isToRight, int numGapCount, bool useBinGenome, int *shift){

    int score = 0;
    // changing value of indelshift
    double indelShiftLocal = indelShift + (numGap*numGapCount)/chunkSize;
    bool binGenome;
    if( useBinGenome ){
        binGenome = checkBinGenome(startLoci, startLoci + this->chunkSize, false, (BGLR ? isToRight : true));
        if(binGenome){
            cntBG++;
            return 1;
        }
    }
    cntCOM++;
    LocalAligner *local;

    if( flag!=0 )
        readString = revComplemACGT(readString);

    int gapPen = -8, misPen = -5, matchPen = 10;
    if(startLoci - indelShiftLocal*readString.size() < 0 || startLoci + indelShiftLocal*readString.size()+this->chunkSize  > genomeLength-1 )
        return 0;
    string genomeString = this->getGenome(startLoci - indelShiftLocal*readString.size(), this->chunkSize + 2*readString.size()*indelShiftLocal );

    string seqConstantB = "";
    int readlen = readString.size();
    for(int i=0;i < (int)(indelShiftLocal*readlen);i++)
        seqConstantB += "A";
    readString += seqConstantB;

    local = new LocalAligner(readString, genomeString, gapPen, misPen, matchPen, (int)(indelShiftLocal*readlen), (int)(editDistance*readlen) );
    local->process();
    local->backtrack();
    //cout//<<local->realCigar << endl;
    //cout<<local->cigar << endl;
    score = local->mScore;
    //AlignmentScores.push_back(score);
    char isAlignable = (score > localAlThreshold ? 2 : 0 );
    if(isAlignable){
        *shift = local->totalGapr - local->totalGapfa;
    }
    delete local;
    return isAlignable;
}

inline int SVAnchor::numberOfExtendableFrags(vector<string> FragsOfRead, vector<ExtInterval>::iterator extensionInterval, int direction , int threadNum ){ //-/ AMEER NEWFUNC
    FragLoci extendFrom = (direction==1 ? extensionInterval->extFrom : extensionInterval->extTo );
    FragLoci extendTo = (direction==1 ? extensionInterval->extTo : extensionInterval->extFrom );
    bool isToRight = ((extendFrom - extendTo)*(extendFrom.Flag==0 ? 1 : -1) < 0 ? true : false );
    //if(isToRight != (direction==1 ? true:false ) )
       // cout<<"Error in SVAnchor 101";
    //-/ ////////////////////////////////////////////////////////////////////////////////////
    //-/ ////////////////////////////////////////////////////////////////////////////////////
    //-/ Genome binary model for Extension Prevention of previously extended repeat regions
    //    if(this region is known as repeat)
    //        then: find if any fragment in this intervals needs extension
    // we need to find a new FragLoci for extensio or an indicator tells that no extension needed here and we can update and close this gap
    // ********************************************** AND CLOSE THIS GAP ********************
    // +
    // a flag of previouslu bowtie alignemnt to check if its multi or not alignable at all (do not worse that much)
    //-/ ////////////////////////////////////////////////////////////////////////////////////
    // so:
    //    if(gapClosed){
    //    }

    char successFlag = 0;
    long long startLoci = extendFrom.Position-1;
    int counter = 0;
    int jumpsCount = 0;
    prevShiftPos.at(threadNum) = 0 ;
    shiftPos.at(threadNum) = 0 ;
    //-//changes
    //cout << "first: "<<startLoci <<endl;
    if( ((abs(extendTo.Fragment - extendFrom.Fragment )-1==0 && extendTo.Flag!=-1 )
         || ( abs(extendTo.Fragment - extendFrom.Fragment )==0 && extendTo.Flag==-1 ) ) )
        return 0;
    do{
        startLoci += (direction*(this->chunkSize + shiftPos.at(threadNum)))*( extendFrom.Flag== 0 ? 1 : -1 );//isForward(extendFrom.Flag);
        successFlag =  ifAlignable( startLoci ,FragsOfRead.at(extendFrom.Fragment + direction*(counter+1) -1 ) ,extendFrom.Flag ,threadNum ,isToRight ,jumpsCount ,(counter==0 && onlyBG ? true:false) );
        counter += ( successFlag ? 1:0 );
        jumpsCount += ( successFlag == 1 ? 1:0);
    }while( successFlag &&
            ((counter < abs(extendTo.Fragment - extendFrom.Fragment )-1 && extendTo.Flag!=-1 )
             || (counter < abs(extendTo.Fragment - extendFrom.Fragment ) && extendTo.Flag==-1 ) ) );//-//changes
    return counter;
}
void* SVAnchor::extensionStepWrapper(void* arg){
    return (((params*)arg)->ptr)->extensionStep(  ((params*)arg));
}

void* SVAnchor::extensionStep( void *threadarg  /*0 <= threadNum < numberofthreads*/ ){ //-/ AMEER NEWFUNC
    int threadNum;
    /* Lock.  */
    pthread_mutex_lock(&(*(params_t*)(threadarg)).mutex);
    /* Work.  */
    threadNum = (*(params_t*)(threadarg)).id;
    /* Unlock and signal completion.  */
    pthread_mutex_unlock(&(*(params_t*)(threadarg)).mutex);
    pthread_cond_signal (&(*(params_t*)(threadarg)).done);
    //pFS("\n -from threadNum: "+convertNumToStr(threadNum)+" +1 ;\n");

    string name =  outputDir+"BreakForExt_"+ convertNumToStr(threadNum + 1) + ".txt";
    ifstream ifstr(name.c_str());
    string line;
    string previous_line;
    string ReadSeq;
    int readLen = 0;
    int completeFrags = 0;

    vector< vector<ExtInterval> >::iterator row;
    vector< vector<ExtInterval> >::iterator rowEnd;
    vector<ExtInterval>::iterator col;
    vector<int>::iterator ind;
    vector<int>::iterator itemp;

    vector<string> FragsOfRead;

    long long lineCounter = 0;
    row = allExtIntervals.begin()+   (threadNum > 0 ?  extensionThreadIntervals[threadNum] : 0 ); // check -1
    rowEnd = (threadNum+1 != this->numberOfThread ? allExtIntervals.begin()+extensionThreadIntervals[threadNum+1]  : allExtIntervals.end() );// check +1
    ind = allExtIndices.begin()+(threadNum > 0 ? extensionThreadIntervals[threadNum] : 0 );// check -1
    cout<<" Thread: "<<threadNum<<" -First read: "<<reads[*(ind)-1].index<<" -Last read: "<< reads[*(ind+(rowEnd-row-1))-1 ].index<<endl;

    for ( ; row != rowEnd; ++row,++ind) {
        //-> Read Retreival
        while (getline(ifstr, line)) {
//            if(row >= rowEnd-2)
//                //pFS("thread "+convertNumToStr(threadNum)+" will be ends"+convertNumToStr(lineCounter));
//                pFS("\nreading read "+convertNumToStr((reads[*(ind)-1].index))+" in thread "+convertNumToStr(threadNum));
            if ( lineCounter % 2 == 0 /*&& (lineCounter/4)+1 == *ind*/ ){
                previous_line = line;
                getline(ifstr,ReadSeq);
                lineCounter+=2;
                break;
            }
            lineCounter++;
        }
        //<- end of Read Retreival

        //-> chunking Read:
        //cout<<"anyerror";
        readLen = (int) ReadSeq.size();
        FragsOfRead.clear();
        completeFrags = floor((readLen-1)/this->chunkSize);
        for( int i=0; i < completeFrags ; i++ )
            FragsOfRead.push_back( ReadSeq.substr(  (i*this->chunkSize)  , this->chunkSize ));
        FragsOfRead.push_back( ReadSeq.substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
        //<- end of chunking Read

        //-> Extending Intervals

        int newFrom = 0, newTo = 0 ;
        for (col = row->begin(); col != row->end(); ++col) {
//            if(*ind == 709)
//                cout<<endl;
            if(col->extFrom.Flag != -1 ){
                // extend from FROM to TO
                newFrom =  col->extFrom.Fragment + numberOfExtendableFrags(FragsOfRead, col , 1 , threadNum); // true: extend from FROM to TO
                //- FASTENING: store all this changes in a third party dataS and then apply all of them together
                if(newFrom != col->extFrom.Fragment ){
                    bool notReverse = col->extFrom.Flag == 0;

                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                    if( !BGLR && newFrom - col->extFrom.Fragment > 3 ){
                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : -1*this->chunkSize  ),
                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? -1 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                        true);
                    }
                    // V4-2BG : Expected (1or2)*chunksize distance with SV BP in 1 BG - extending toRight: -1chunk in right -2chink in left and wiseversa
                    if( BGLR && newFrom - col->extFrom.Fragment > 2 ){
                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : 0  ),
                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 0 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                        true);
                        changeBinGenome(col->extFrom.Position + (notReverse? 1*this->chunkSize : -1*this->chunkSize  ),
                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment +        -1           )*this->chunkSize + prevShiftPos.at(threadNum) ),
                                        false);
                    }
                    // V3 - wrong yet
//                    changeBinGenome(col->extFrom.Position + (notReverse? this->chunkSize : 0 ),
//                                    col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 0 : -1))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    true);

                    // V2 //prevShiftPos
//                    changeBinGenome(col->extFrom.Position + this->chunkSize,
//                                    col->extFrom.Position + this->chunkSize + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment)*this->chunkSize + shiftPos.at(threadNum) ),
//                                    true);
//                    changeBinGenome(col->extFrom.Position,
//                                    col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment)*this->chunkSize + shiftPos.at(threadNum) ),
//                                    false);
                    // V wrong
//                    changeBinGenome(col->extFrom.Position + (notReverse? 0 : this->chunkSize),
//                                    col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 1 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    notReverse);
//                    changeBinGenome(col->extFrom.Position + (notReverse? 0 : this->chunkSize),
//                                    col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 0 :-1))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    !notReverse);

                    if(reads[*(ind)-1].lastFragment != col->extFrom.Fragment ){
                        itemp = find(reads[*(ind)-1].lastFragClus2.begin(),reads[*(ind)-1].lastFragClus2.end() , col->extFrom.Fragment );
                        if(itemp == reads[*(ind)-1].lastFragClus2.end() ){
                            pFS("\n Err(11)- : fragment not found\n");// ERROR!
                            cout<<threadNum<<" - "<<reads[*(ind)-1].index<<endl;
                        }
                        else{
                            int tt = distance(reads[*(ind)-1].lastFragClus2.begin(),itemp);
                            reads[*(ind)-1].lastPosClus2.at(tt) +=
                                    (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[*(ind)-1].flagClus2.at(tt)==0?1:-1 );
                            *itemp = (itemp!=reads[*(ind)-1].lastFragClus2.end() ?  newFrom : *itemp );
                        }
                    }else{
                        reads[*(ind)-1].lastPosition += (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[*(ind)-1].flag==0?1:-1 )  ;
                        reads[*(ind)-1].lastFragment = newFrom;
                        //pFS("2");
                    }
                    col->extFrom.Fragment = newFrom;
                }
            }
            if(col->extTo.Flag != -1 ){
                // extend from TO to FROM
                newTo = col->extTo.Fragment - numberOfExtendableFrags(FragsOfRead, col ,-1 ,threadNum); // false: extend from TO to FROM
                //- FASTENING: store all this changes in a third party dataS and then apply all of them together
                if(newTo != col->extTo.Fragment){
                    bool notReverse = col->extTo.Flag == 0;

                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                    if( !BGLR && col->extTo.Fragment - newTo > 3 ){
                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 2*this->chunkSize ),
                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : -1))*this->chunkSize + shiftPos.at(threadNum) ),
                                        true);
                    }
                    if( BGLR && col->extTo.Fragment - newTo > 2 ){
                        changeBinGenome(col->extTo.Position + (notReverse? 0 : 2*this->chunkSize ),
                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
                                        true);
                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 1*this->chunkSize ),
                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  +            -1       )*this->chunkSize + shiftPos.at(threadNum) ),
                                        false);
                    }

                    // V3

//                    changeBinGenome(col->extTo.Position + (notReverse? 0 : this->chunkSize ),
//                                    col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -1 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    true);
                    // V2 //prevShiftPos
//                    changeBinGenome(col->extTo.Position + this->chunkSize,
//                                    col->extTo.Position + this->chunkSize - (notReverse? 1 : -1)*((col->extTo.Fragment - newTo)*this->chunkSize + shiftPos.at(threadNum) ),
//                                    true);
//                    changeBinGenome(col->extTo.Position ,
//                                    col->extTo.Position - (notReverse? 1 : -1)*((col->extTo.Fragment - newTo)*this->chunkSize + shiftPos.at(threadNum) ),
//                                    false);
                    // V wrong
//                    changeBinGenome(col->extTo.Position + (notReverse? this->chunkSize : 0 ),
//                                    col->extTo.Position - (notReverse? 1 : -1)*((col->extTo.Fragment - newTo + (notReverse? 0 : 1))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    !notReverse);

//                    changeBinGenome(col->extTo.Position + (notReverse? this->chunkSize : 0 ),
//                                    col->extTo.Position - (notReverse? 1 : -1)*((col->extTo.Fragment - newTo + (notReverse? -1 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
//                                    notReverse);
                    if(reads[*(ind)-1].firstFragment != col->extTo.Fragment ){
                        itemp = find(reads[*(ind)-1].firstFragClus2.begin(),reads[*(ind)-1].firstFragClus2.end() , col->extTo.Fragment );
                        if(itemp==reads[*(ind)-1].firstFragClus2.end()){ // ERROR!
                             pFS(" Err(21) : fragment not found\n");// ERROR!
                             cout<<threadNum<<" - "<<reads[*(ind)-1].index<<endl;
                        }
                        else{
                            int tt = distance(reads[*(ind)-1].firstFragClus2.begin(),itemp);
                            reads[*(ind)-1].firstPosClus2.at(tt) +=
                                    (newTo - col->extTo.Fragment)*this->chunkSize*(reads[*(ind)-1].flagClus2.at(tt)==0?1:-1 ) ;
                            *itemp = (itemp!=reads[*(ind)-1].firstFragClus2.end() ?  newTo : *itemp );
                        }

                    }else {
                        reads[*(ind)-1].firstPosition += (newTo - col->extTo.Fragment)*this->chunkSize*(reads[*(ind)-1].flag==0?1:-1 )  ;
                        reads[*(ind)-1].firstFragment = newTo;
                        //pFS("2");
                    }
                    col->extTo.Fragment = newTo;
                }
            }
        }
        //<- end of Extending Inetrvals
    }
    //cout<<"\n last read in thread "<<threadNum<< " proccessed with id :"<<reads[*(ind-1)-1].index<<endl;
}
void SVAnchor::extensionStep2(ifstream* file,long long* GLC,int* idt, int* readsFileCount){

    if(*readsFileCount == 1){
        //pFS("\nthread"+convertNumToStr((*idt))+"starts");
        int threadNum = *idt;
        string line;
        string line2;
        vector<long long> readIndex;
        vector<string> readSeq;
        //vector<string> readHeader;
        // extend-needed cluster numbers in each read.
        vector<ExtInterval> colExtInterval;
        vector<ExtInterval>::iterator col;
        vector<FragLoci> firsts;
        vector<FragLoci> lasts;
        bool flag = false;
        ExtInterval temp;
        vector<string> FragsOfRead;
        vector<int>::iterator itemp;
        int readLen = 0;

        int bufferCounter;
        do {
            //new-edit:
            readIndex.clear();
            readSeq.clear();

            bufferCounter = 0;
            // ////////////////// **********************************|
            // ////////////////// STEP 1 ---------------------------|
            // ////////////////// Threads compete to read from file |
            // ////////////////// **********************************|
            mainMutex.lock();
            while(bufferCounter < 4*MAX_BUFFER_SIZE && (*GLC) < 4*this->numberOfReads) {

                getline(*file,line);
                if( (*GLC) % 4 == 0 ){
                    //pFS("\n thread"+convertNumToStr(*idt)+" reading read "+convertNumToStr((*GLC)/4));
                    getline(*file,line2);
                    readIndex.push_back((*GLC)/4);
                    //readHeader.push_back(line);
                    readSeq.push_back(line2);
                    (*GLC)++;
                    bufferCounter++;
                }
                (*GLC)++;
                bufferCounter++;
            }
            mainMutex.unlock();

            // ////////////////// **********************************|
            // ////////////////// STEP 2 ---------------------------|
            // ////////////////// Threads seperated processing      |
            // ////////////////// **********************************|
            long long rI = 0;
            for(int k = 0 ; k < readIndex.size() ; k++){
                rI = readIndex[k];
                if((abs(reads[rI].firstFragment-reads[rI].lastFragment)+1)*this->chunkSize == this->reads[rI].length )
                    ;
                else{
                    firsts.clear();
                    lasts.clear();
                    colExtInterval.clear();
                    flag = false;

                    firsts.push_back({ reads[rI].firstFragment, reads[rI].flag, reads[rI].firstPosition });
                    lasts.push_back({reads[rI].lastFragment, reads[rI].flag, reads[rI].lastPosition});
                    for(int j = 1 ; j < this->reads[rI].firstPosClus2.size() ; j++ ){
                        firsts.push_back({reads[rI].firstFragClus2.at(j),reads[rI].flagClus2.at(j),reads[rI].firstPosClus2.at(j)});
                        lasts.push_back({reads[rI].lastFragClus2.at(j),reads[rI].flagClus2.at(j),reads[rI].lastPosClus2.at(j)});
                    }
                    sort(firsts.begin(),firsts.end());
                    sort(lasts.begin(),lasts.end());
                    if(firsts.front() > 1 ){
                        // pushback (1,firsts.front()) for this read
                        if(flag == false){
                            flag = true;
                            this->numExtendables++;
                            //allExtIndices.push_back(this->reads[rI].index);
                        }
                        temp.extFrom={1/*Frag*/,-1/*Flag*/,0/*Position*/} ; temp.extTo = firsts.front();
                        colExtInterval.push_back(temp);
                    }
                    if(lasts.back() < (this->reads[rI].length + this->chunkSize - 1 )/this->chunkSize ){
                        // pushback (lasts.back(),reads[i].length / this->chunkSize ) for this read
                        if(flag == false){
                            flag = true;
                            this->numExtendables++;
                            //allExtIndices.push_back(this->reads[rI].index);
                        }
                        temp.extFrom=lasts.back() ; temp.extTo= {(this->reads[rI].length + this->chunkSize - 1 )/this->chunkSize,-1,0};
                        colExtInterval.push_back(temp);
                    }

                    for (std::vector<FragLoci>::iterator it=firsts.begin()+1,it2=lasts.begin(); it!=firsts.end(); ++it,++it2)
                        if( *it - *it2 > 1 ){ //-/ **************** OR COULD BE 2 INSTEAD OF 1
                            // pushback (*it2,*it) for this read
                            if(flag == false){
                                flag = true;
                                this->numExtendables++;
                                //allExtIndices.push_back(this->reads[rI].index);
                            }
                            temp.extFrom=*it2 ; temp.extTo= *it;
                            colExtInterval.push_back(temp);
                        }/*else if(*it - *it2 < 0) { pFS(convertNumToStr(it->Fragment) );pFS("->");pFS(convertNumToStr(it2->Fragment) ); pFS("Error2");
                            for (std::vector<FragLoci>::iterator it3=firsts.begin(),it4=lasts.begin(); it3!=firsts.end(); ++it3,++it4)
                                pFS(convertNumToStr(this->reads[rI].index)+":"+convertNumToStr(it3->Fragment)+"+"+convertNumToStr(it4->Fragment)+" E2\n" );
                            pFS("\n---\n");
                        }*/
                    if( colExtInterval.size() ){
                        // //////// *************************
                        // //////// if Extension is probable:
                        // //////// *************************
                        int completeFrags = 0;
                        readLen = (int) readSeq[k].size();
                        FragsOfRead.clear();
                        completeFrags = floor((readLen-1)/this->chunkSize);
                        for( int i=0; i < completeFrags ; i++ )
                            FragsOfRead.push_back( readSeq[k].substr(  (i*this->chunkSize)  , this->chunkSize ));
                        FragsOfRead.push_back( readSeq[k].substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
                        int newFrom = 0, newTo = 0 ;
                        for (col = colExtInterval.begin(); col != colExtInterval.end(); ++col) {
                            if(col->extFrom.Flag != -1 ){
                                // extend from FROM to TO
                                if(rI == 45975 && col->extFrom.Fragment == 20 )
                                    cout<<"Err here\n";
                                newFrom =  col->extFrom.Fragment + numberOfExtendableFrags(FragsOfRead, col , 1 , threadNum); // true: extend from FROM to TO
                                if(newFrom != col->extFrom.Fragment ){
                                    bool notReverse = col->extFrom.Flag == 0;
                                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                                    if( !BGLR && newFrom - col->extFrom.Fragment > 3 ){
                                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : -1*this->chunkSize  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? -1 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        true);
                                    }
                                    // V4-2BGs : Expected (1or2)*chunksize distance with SV BP in 1 BG - extending toRight: -1chunk in right -2chink in left and wiseversa
                                    if( BGLR && newFrom - col->extFrom.Fragment > 2 ){
                                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : 0  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 0 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        true);
                                        changeBinGenome(col->extFrom.Position + (notReverse? 1*this->chunkSize : -1*this->chunkSize  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment +        -1           )*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        false);
                                    }

                                    if(reads[rI].lastFragment != col->extFrom.Fragment ){
                                        itemp = find(reads[rI].lastFragClus2.begin(),reads[rI].lastFragClus2.end() , col->extFrom.Fragment );
                                        if(itemp == reads[rI].lastFragClus2.end() ){
                                            pFS("\n Err(11)- : fragment not found\n");// ERROR!
                                            cout<<threadNum<<" - "<<reads[rI].index<<endl;
                                        }
                                        else{
                                            int tt = distance(reads[rI].lastFragClus2.begin(),itemp);
                                            extChangesMutex.lock();
                                            reads[rI].lastPosClus2.at(tt) +=
                                                    (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[rI].flagClus2.at(tt)==0?1:-1 );
                                            *itemp = (itemp!=reads[rI].lastFragClus2.end() ?  newFrom : *itemp );
                                            extChangesMutex.unlock();
                                        }
                                    }else{
                                        extChangesMutex.lock();
                                        reads[rI].lastPosition += (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[rI].flag==0?1:-1 )  ;
                                        reads[rI].lastFragment = newFrom;
                                        extChangesMutex.unlock();
                                        //pFS("2");
                                    }
                                    col->extFrom.Fragment = newFrom;
                                }
                            }
                            if(col->extTo.Flag != -1 ){
                                // extend from TO to FROM
                                newTo = col->extTo.Fragment - numberOfExtendableFrags(FragsOfRead, col ,-1 ,threadNum); // false: extend from TO to FROM
                                if(newTo != col->extTo.Fragment){
                                    bool notReverse = col->extTo.Flag == 0;

                                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                                    if( !BGLR && col->extTo.Fragment - newTo > 3 ){
                                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 2*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : -1))*this->chunkSize + shiftPos.at(threadNum) ),
                                                        true);
                                    }
                                    // V4-2BGs : Expected (1or2)*chunksize distance with SV BP in 1 BG - extending toRight: -1chunk in right -2chink in left and wiseversa
                                    if( BGLR && col->extTo.Fragment - newTo > 2 ){
                                        changeBinGenome(col->extTo.Position + (notReverse? 0 : 2*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
                                                        true);
                                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 1*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  +            -1       )*this->chunkSize + shiftPos.at(threadNum) ),
                                                        false);
                                    }

                                    if(reads[rI].firstFragment != col->extTo.Fragment ){
                                        itemp = find(reads[rI].firstFragClus2.begin(),reads[rI].firstFragClus2.end() , col->extTo.Fragment );
                                        if(itemp==reads[rI].firstFragClus2.end()){ // ERROR!
                                            pFS(" Err(21) : fragment not found\n");// ERROR!
                                            cout<<threadNum<<" - "<<reads[rI].index<<endl;
                                        }
                                        else{
                                            int tt = distance(reads[rI].firstFragClus2.begin(),itemp);
                                            extChangesMutex.lock();
                                            reads[rI].firstPosClus2.at(tt) +=
                                                    (newTo - col->extTo.Fragment)*this->chunkSize*(reads[rI].flagClus2.at(tt)==0?1:-1 ) ;
                                            *itemp = (itemp!=reads[rI].firstFragClus2.end() ?  newTo : *itemp );
                                            extChangesMutex.unlock();
                                        }

                                    }else {
                                        extChangesMutex.lock();
                                        reads[rI].firstPosition += (newTo - col->extTo.Fragment)*this->chunkSize*(reads[rI].flag==0?1:-1 )  ;
                                        reads[rI].firstFragment = newTo;
                                        extChangesMutex.unlock();
                                        //pFS("2");
                                    }
                                    col->extTo.Fragment = newTo;
                                }
                            }
                        }
                        ///
                        /// ////
                        /// ////

                    }
                }
            }

            //until EOF was not reached
        } while( (*GLC) < 4*this->numberOfReads );

        pFS(" thread "+ convertNumToStr(*idt) +" done.\n");

    }
    if(*readsFileCount > 1){
        //pFS("\nthread"+convertNumToStr((*idt))+"starts");
        int threadNum = *idt;
        string line;
        string line2;
        vector<long long> readIndex;
        vector<string> readSeq;
        //vector<string> readHeader;
        // extend-needed cluster numbers in each read.
        vector<ExtInterval> colExtInterval;
        vector<ExtInterval>::iterator col;
        vector<FragLoci> firsts;
        vector<FragLoci> lasts;
        bool flag = false;
        ExtInterval temp;
        vector<string> FragsOfRead;
        vector<int>::iterator itemp;
        int readLen = 0;
        bool notEOF = true;
        int bufferCounter;
        do {
            //new-edit:
            readIndex.clear();
            readSeq.clear();
            bufferCounter = 0;
            // ////////////////// **********************************|
            // ////////////////// STEP 1 ---------------------------|
            // ////////////////// Threads compete to read from file |
            // ////////////////// **********************************|
            mainMutex.lock();
            while( bufferCounter < 4*MAX_BUFFER_SIZE && notEOF ) {

                if(!getline(*file,line)){
                    notEOF = false;
                    break;
                }
                if( (*GLC) % 4 == 0 ){
                    //pFS("\n thread"+convertNumToStr(*idt)+" reading read "+convertNumToStr((*GLC)/4));
                    if(!getline(*file,line2)){
                        cout<<"\n WRONG READ FILE so skipped the remaining"<<endl;
                        notEOF = false;
                        break;
                    }
                    readIndex.push_back((*GLC)/4);
                    //readHeader.push_back(line);
                    readSeq.push_back(line2);
                    (*GLC)++;
                    bufferCounter++;
                }
                (*GLC)++;
                bufferCounter++;
            }
            mainMutex.unlock();

            // ////////////////// **********************************|
            // ////////////////// STEP 2 ---------------------------|
            // ////////////////// Threads seperated processing      |
            // ////////////////// **********************************|
            long long rI = 0;
            for(int k = 0 ; k < readIndex.size() ; k++){
                rI = readIndex[k];
                if((abs(reads[rI].firstFragment-reads[rI].lastFragment)+1)*this->chunkSize == this->reads[rI].length )
                    ;
                else{
                    firsts.clear();
                    lasts.clear();
                    colExtInterval.clear();
                    flag = false;

                    firsts.push_back({ reads[rI].firstFragment, reads[rI].flag, reads[rI].firstPosition });
                    lasts.push_back({reads[rI].lastFragment, reads[rI].flag, reads[rI].lastPosition});
                    for(int j = 1 ; j < this->reads[rI].firstPosClus2.size() ; j++ ){
                        firsts.push_back({reads[rI].firstFragClus2.at(j),reads[rI].flagClus2.at(j),reads[rI].firstPosClus2.at(j)});
                        lasts.push_back({reads[rI].lastFragClus2.at(j),reads[rI].flagClus2.at(j),reads[rI].lastPosClus2.at(j)});
                    }
                    sort(firsts.begin(),firsts.end());
                    sort(lasts.begin(),lasts.end());
                    if(firsts.front() > 1 ){
                        // pushback (1,firsts.front()) for this read
                        if(flag == false){
                            flag = true;
                            this->numExtendables++;
                            //allExtIndices.push_back(this->reads[rI].index);
                        }
                        temp.extFrom={1/*Frag*/,-1/*Flag*/,0/*Position*/} ; temp.extTo = firsts.front();
                        colExtInterval.push_back(temp);
                    }
                    if(lasts.back() < (this->reads[rI].length + this->chunkSize - 1 )/this->chunkSize ){
                        // pushback (lasts.back(),reads[i].length / this->chunkSize ) for this read
                        if(flag == false){
                            flag = true;
                            this->numExtendables++;
                            //allExtIndices.push_back(this->reads[rI].index);
                        }
                        temp.extFrom=lasts.back() ; temp.extTo= {(this->reads[rI].length + this->chunkSize - 1 )/this->chunkSize,-1,0};
                        colExtInterval.push_back(temp);
                    }

                    for (std::vector<FragLoci>::iterator it=firsts.begin()+1,it2=lasts.begin(); it!=firsts.end(); ++it,++it2)
                        if( *it - *it2 > 1 ){ //-/ **************** OR COULD BE 2 INSTEAD OF 1
                            // pushback (*it2,*it) for this read
                            if(flag == false){
                                flag = true;
                                this->numExtendables++;
                                //allExtIndices.push_back(this->reads[rI].index);
                            }
                            temp.extFrom=*it2 ; temp.extTo= *it;
                            colExtInterval.push_back(temp);
                        }/*else if(*it - *it2 < 0) { pFS(convertNumToStr(it->Fragment) );pFS("->");pFS(convertNumToStr(it2->Fragment) ); pFS("Error2");
                            for (std::vector<FragLoci>::iterator it3=firsts.begin(),it4=lasts.begin(); it3!=firsts.end(); ++it3,++it4)
                                pFS(convertNumToStr(this->reads[rI].index)+":"+convertNumToStr(it3->Fragment)+"+"+convertNumToStr(it4->Fragment)+" E2\n" );
                            pFS("\n---\n");
                        }*/
                    if( colExtInterval.size() ){
                        // //////// *************************
                        // //////// if Extension is probable:
                        // //////// *************************
                        int completeFrags = 0;
                        readLen = (int) readSeq[k].size();
                        FragsOfRead.clear();
                        completeFrags = floor((readLen-1)/this->chunkSize);
                        for( int i=0; i < completeFrags ; i++ )
                            FragsOfRead.push_back( readSeq[k].substr(  (i*this->chunkSize)  , this->chunkSize ));
                        FragsOfRead.push_back( readSeq[k].substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
                        int newFrom = 0, newTo = 0 ;
                        for (col = colExtInterval.begin(); col != colExtInterval.end(); ++col) {
                            if(col->extFrom.Flag != -1 ){
                                // extend from FROM to TO
                                if(rI == 45975 && col->extFrom.Fragment == 20 )
                                    cout<<"Err here\n";
                                newFrom =  col->extFrom.Fragment + numberOfExtendableFrags(FragsOfRead, col , 1 , threadNum); // true: extend from FROM to TO
                                if(newFrom != col->extFrom.Fragment ){
                                    bool notReverse = col->extFrom.Flag == 0;
                                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                                    if( !BGLR && newFrom - col->extFrom.Fragment > 3 ){
                                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : -1*this->chunkSize  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? -1 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        true);
                                    }
                                    // V4-2BGs : Expected (1or2)*chunksize distance with SV BP in 1 BG - extending toRight: -1chunk in right -2chink in left and wiseversa
                                    if( BGLR && newFrom - col->extFrom.Fragment > 2 ){
                                        changeBinGenome(col->extFrom.Position + (notReverse? 2*this->chunkSize : 0  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment + (notReverse? 0 : -2))*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        true);
                                        changeBinGenome(col->extFrom.Position + (notReverse? 1*this->chunkSize : -1*this->chunkSize  ),
                                                        col->extFrom.Position + (notReverse? 1 : -1)*((newFrom - col->extFrom.Fragment +        -1           )*this->chunkSize + prevShiftPos.at(threadNum) ),
                                                        false);
                                    }

                                    if(reads[rI].lastFragment != col->extFrom.Fragment ){
                                        itemp = find(reads[rI].lastFragClus2.begin(),reads[rI].lastFragClus2.end() , col->extFrom.Fragment );
                                        if(itemp == reads[rI].lastFragClus2.end() ){
                                            pFS("\n Err(11)- : fragment not found\n");// ERROR!
                                            cout<<threadNum<<" - "<<reads[rI].index<<endl;
                                        }
                                        else{
                                            int tt = distance(reads[rI].lastFragClus2.begin(),itemp);
                                            extChangesMutex.lock();
                                            reads[rI].lastPosClus2.at(tt) +=
                                                    (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[rI].flagClus2.at(tt)==0?1:-1 );
                                            *itemp = (itemp!=reads[rI].lastFragClus2.end() ?  newFrom : *itemp );
                                            extChangesMutex.unlock();
                                        }
                                    }else{
                                        extChangesMutex.lock();
                                        reads[rI].lastPosition += (newFrom - col->extFrom.Fragment)*this->chunkSize*(reads[rI].flag==0?1:-1 )  ;
                                        reads[rI].lastFragment = newFrom;
                                        extChangesMutex.unlock();
                                        //pFS("2");
                                    }
                                    col->extFrom.Fragment = newFrom;
                                }
                            }
                            if(col->extTo.Flag != -1 ){
                                // extend from TO to FROM
                                newTo = col->extTo.Fragment - numberOfExtendableFrags(FragsOfRead, col ,-1 ,threadNum); // false: extend from TO to FROM
                                if(newTo != col->extTo.Fragment){
                                    bool notReverse = col->extTo.Flag == 0;

                                    // V4-1BG : Expected 2*chunksize distance with maybeSV BP in 1 BG ->2chunksize distance in both ends
                                    if( !BGLR && col->extTo.Fragment - newTo > 3 ){
                                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 2*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : -1))*this->chunkSize + shiftPos.at(threadNum) ),
                                                        true);
                                    }
                                    // V4-2BGs : Expected (1or2)*chunksize distance with SV BP in 1 BG - extending toRight: -1chunk in right -2chink in left and wiseversa
                                    if( BGLR && col->extTo.Fragment - newTo > 2 ){
                                        changeBinGenome(col->extTo.Position + (notReverse? 0 : 2*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  + (notReverse? -2 : 0))*this->chunkSize + shiftPos.at(threadNum) ),
                                                        true);
                                        changeBinGenome(col->extTo.Position + (notReverse? -1*this->chunkSize  : 1*this->chunkSize ),
                                                        col->extTo.Position + (notReverse? -1 : 1)*((col->extTo.Fragment - newTo  +            -1       )*this->chunkSize + shiftPos.at(threadNum) ),
                                                        false);
                                    }

                                    if(reads[rI].firstFragment != col->extTo.Fragment ){
                                        itemp = find(reads[rI].firstFragClus2.begin(),reads[rI].firstFragClus2.end() , col->extTo.Fragment );
                                        if(itemp==reads[rI].firstFragClus2.end()){ // ERROR!
                                            pFS(" Err(21) : fragment not found\n");// ERROR!
                                            cout<<threadNum<<" - "<<reads[rI].index<<endl;
                                        }
                                        else{
                                            int tt = distance(reads[rI].firstFragClus2.begin(),itemp);
                                            extChangesMutex.lock();
                                            reads[rI].firstPosClus2.at(tt) +=
                                                    (newTo - col->extTo.Fragment)*this->chunkSize*(reads[rI].flagClus2.at(tt)==0?1:-1 ) ;
                                            *itemp = (itemp!=reads[rI].firstFragClus2.end() ?  newTo : *itemp );
                                            extChangesMutex.unlock();
                                        }

                                    }else {
                                        extChangesMutex.lock();
                                        reads[rI].firstPosition += (newTo - col->extTo.Fragment)*this->chunkSize*(reads[rI].flag==0?1:-1 )  ;
                                        reads[rI].firstFragment = newTo;
                                        extChangesMutex.unlock();
                                        //pFS("2");
                                    }
                                    col->extTo.Fragment = newTo;
                                }
                            }
                        }
                        ///
                        /// ////
                        /// ////

                    }
                }
            }

            //until EOF was not reached
        } while( notEOF );

        pFS(" thread "+ convertNumToStr(*idt) +" done.\n");

    }
}

void SVAnchor::breakReadsForExtension(bool allReads){

    int ReadCount = (allReads ? numberOfReads : allExtIndices.size());
    //cout<<"\nReadCount"<<ReadCount;
    int eachThread2 = ReadCount / this->numberOfThread;
    //cout<<"\neachThread2"<<eachThread2;
    int Rem = ReadCount % this->numberOfThread;
    //cout<<"\nRem "<<Rem ;
    int eachThread1 = (eachThread2+1);
    //cout<<"\neachThread1"<<eachThread1;
    int allThreads1 = eachThread1*Rem;
    //cout<<"\nallThreads1"<<allThreads1;
    long long lineCounter = 0;
    string name = readName;
    ifstream ifstr(name.c_str());

    string previous_line;
    string ReadSeq;
    string line;

    vector< vector<ExtInterval> >::iterator row;
    vector<int>::iterator ind;

    ofstream* myBreakReads;
    myBreakReads = new ofstream[this->numberOfThread];

    for (int thread = 0; thread < this->numberOfThread; thread++) {
        ostringstream filename;
        filename << outputDir<<"BreakForExt_"<< thread + 1 << ".txt";
        myBreakReads[thread].open(filename.str().c_str());
    }
    int k = -1;
    int index = 0,previousIndex = 0;
    extensionThreadIntervals =  new int[numberOfThread];
    extensionThreadIntervals[0]=0;
    for (row = allExtIntervals.begin(),ind = allExtIndices.begin() ; row != allExtIntervals.end(); ++row,++ind) {
        //-> Read Retreival
        while (getline(ifstr, line)) {
            if ( lineCounter % 4 == 0 && (lineCounter/4)+1 == *ind ){
                k++;
                previous_line = line;
                getline(ifstr,ReadSeq);
                lineCounter+=2;
                index = ( k < allThreads1 ? k/eachThread1 : (k-allThreads1)/eachThread2 + Rem );
                if(index!=previousIndex){
                    pFS("\n****\nThread "+convertNumToStr(index-1)+" : "+convertNumToStr(k)+" - "+convertNumToStr(allExtIndices.at(k)-1));
                    extensionThreadIntervals[index]=k;
                    //cout<<endl<< extensionThreadIntervals[index]<<endl ;
                }
                previousIndex = index;
                myBreakReads[index]<<previous_line<<endl<<ReadSeq<<endl;
                break;
            }
            lineCounter++;
        }
    }
    for (int thread = 0; thread < this->numberOfThread; thread++)
        myBreakReads[thread].close();
    //myOSams[threadCounter] << line << endl

}

void SVAnchor::keepReadsForExtension(/*string alter*/) { //-/ AMEER NEWFUNC
    allExtIndices.clear();
    allExtIntervals.clear();
    numExtendables = 0;


    vector<ExtInterval> colExtInterval;
    vector<FragLoci> firsts;
    vector<FragLoci> lasts;
    bool flag = false;
    ExtInterval temp;
    for (long long i = 0; i < this->numberOfReads; i++) {
        //-/ Ameerosein edits
        firsts.clear();
        lasts.clear();
        //if((abs(reads[i].firstFragment-reads[i].lastFragment)+1) == (this->reads[i].length + this->chunkSize - 1 )/this->chunkSize )
         if((abs(reads[i].firstFragment-reads[i].lastFragment)+1)*this->chunkSize == this->reads[i].length )
            ;
        else{
            firsts.push_back({ reads[i].firstFragment, reads[i].flag, reads[i].firstPosition });
            lasts.push_back({reads[i].lastFragment, reads[i].flag, reads[i].lastPosition});
            for(int j = 1 ; j < this->reads[i].firstPosClus2.size() ; j++ ){
                firsts.push_back({reads[i].firstFragClus2.at(j),reads[i].flagClus2.at(j),reads[i].firstPosClus2.at(j)});
                lasts.push_back({reads[i].lastFragClus2.at(j),reads[i].flagClus2.at(j),reads[i].lastPosClus2.at(j)});
            }
            sort(firsts.begin(),firsts.end());
            sort(lasts.begin(),lasts.end());
            //-/ ***************************************************************
            //-/ **************** OR COULD BE 2 INSTEAD OF 1
            //-/ **************** THE QUESTION IS: distance 2 or distance 1 ???
            //-/ ****************
            flag = false;
            colExtInterval.clear();
            if(firsts.front() > 1 ){
                // pushback (1,firsts.front()) for this read
                if(flag == false){
                    flag = true;
                    this->numExtendables++;
                    allExtIndices.push_back(this->reads[i].index);
                }
                temp.extFrom={1/*Frag*/,-1/*Flag*/,0/*Position*/} ; temp.extTo = firsts.front();
                colExtInterval.push_back(temp);
                // allExtIntervals.at( this->numExtendables-1 ).push_back( temp );

            }
            if(lasts.back() < (this->reads[i].length + this->chunkSize - 1 )/this->chunkSize ){
                // pushback (lasts.back(),reads[i].length / this->chunkSize ) for this read
                if(flag == false){
                    flag = true;
                    this->numExtendables++;
                    allExtIndices.push_back(this->reads[i].index);
                }
                temp.extFrom=lasts.back() ; temp.extTo= {(this->reads[i].length + this->chunkSize - 1 )/this->chunkSize,-1,0};
                colExtInterval.push_back(temp);
                //allExtIntervals.at(this->numExtendables-1).push_back(temp);
            }
            //std::vector<vector<ExtInterval>>::iterator numExtendables2 = allExtIndices.begin(); // not good because after every pushback a new allocation may results and ...
            for (std::vector<FragLoci>::iterator it=firsts.begin()+1,it2=lasts.begin(); it!=firsts.end(); ++it,++it2)
//                std::cout << ' ' << *it;
                if( *it - *it2 > 1 ){ //-/ **************** OR COULD BE 2 INSTEAD OF 1
                    // pushback (*it2,*it) for this read
                    if(flag == false){
                        flag = true;
                        this->numExtendables++;
                        allExtIndices.push_back(this->reads[i].index);
                    }
                    temp.extFrom=*it2 ; temp.extTo= *it;
                    colExtInterval.push_back(temp);
                    //allExtIntervals.at(this->numExtendables-1).push_back(temp);
                    //(numExtendables2++)
                }/*else if(*it - *it2 < 0) { pFS(convertNumToStr(it->Fragment) );pFS("->");pFS(convertNumToStr(it2->Fragment) ); pFS("Error2");
                    for (std::vector<FragLoci>::iterator it3=firsts.begin(),it4=lasts.begin(); it3!=firsts.end(); ++it3,++it4)
                        pFS(convertNumToStr(this->reads[i].index)+":"+convertNumToStr(it3->Fragment)+"+"+convertNumToStr(it4->Fragment)+" E2\n" );
                    pFS("\n---\n");
                }*/
            if( colExtInterval.size() )
                allExtIntervals.push_back(colExtInterval);
        }

    }
    return;
}
void SVAnchor::inClusterExtThread(ifstream* file, long long* GLC, int* idt, int *readsFileCount) {
    if( *readsFileCount == 1){
        string line;
        string line2;
        long long rI = 0 ; // temp read index
        vector<long long> readIndex;
        vector<string> readSeq;
        //vector<string> readHeader;
        vector<vector<int>> readExtClusters; // extend-needed cluster numbers in each read.
        int readLen = 0 ;

        //string outputDir = "/home/ameer/SVAnchoring/SV_out2/";
        //string name2 = outputDir+"read_thread"+convertNumToStr(*idt);
        //ofstream ofstr(name2.c_str());

        //    mutex->lock();
        //    cout<<"welcome to "<< (*idt) <<" thread."<<endl;
        //    mutex->unlock();

        int i;
        do {
            i = 0;
            // only 1 concurrent reader
            vector<int> extendingClusters; // extend-needed cluster numbers in each read.
            // ////////////////// **********************************|
            // ////////////////// STEP 1 ---------------------------|
            // ////////////////// Threads compete to read from file |
            // ////////////////// **********************************|
            mainMutex.lock();
            while(i < 4*MAX_BUFFER_SIZE && *GLC < 4*this->numberOfReads) {

                // ////////////////// ****************************************************************************************************
                // ////////////////// this Processes could be done in second step where there is no lock (if the bottleneck is first step)
                // ////////////////// just pushback the read index and then ...

                getline(*file,line);
                if( *GLC % 4 == 0 ){
                    rI = *GLC / 4;
                    extendingClusters.clear();
                    if( !checkBinGenome( min(reads[rI].firstPosition ,reads[rI].lastPosition)+chunkSize  , max(reads[rI].firstPosition ,reads[rI].lastPosition),true ) )
                        extendingClusters.push_back(0);
                    //cout<<"\nafter0\n";
                    for(int j = 1 ; j< reads[rI].firstFragClus2.size() ; j++)
                        if(!checkBinGenome( min(reads[rI].firstPosClus2.at(j),reads[rI].lastPosClus2.at(j))+chunkSize, max(reads[rI].firstPosClus2.at(j),reads[rI].lastPosClus2.at(j)) ,true ) )
                            extendingClusters.push_back(j);

                    getline(*file,line2);
                    if(extendingClusters.size()){
                        readIndex.push_back(rI);
                        //readHeader.push_back(line);
                        readSeq.push_back(line2);
                        readExtClusters.push_back(extendingClusters);
                    }
                    (*GLC)++;
                    i++;
                }

                (*GLC)++;
                i++;

            }
            mainMutex.unlock();
            // ////////////////// **********************************
            // ////////////////// STEP 2 ---------------------------
            // ////////////////// Threads seperated processing
            // ////////////////// **********************************
            // this process may takes no time in comparison with STEP1 because maybe there is no read to extend inClusterly,
            // so we cant manage STEP1 the way to collect MAX_BUFFER_SIZE process needed reads instead of MBS (un)/normal reads.
            // execute work for every line
            int completeFrags = 0;
            vector<string> FragsOfRead;
            for (int k = 0; k < readIndex.size(); k++) {
                // work with readSeq and readHeader and readExtClusters
                readLen = (int) readSeq[k].size();
                FragsOfRead.clear();
                completeFrags = floor((readLen-1)/this->chunkSize);
                for( int j = 0; j < completeFrags ; j++ )
                    FragsOfRead.push_back( readSeq[k].substr(  (j*this->chunkSize)  , this->chunkSize ));
                FragsOfRead.push_back( readSeq[k].substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
                int j = 0;
                if( readExtClusters[k].at(0) == 0 ){
                    inExtendability( FragsOfRead ,reads[readIndex[k]].firstPosition ,reads[readIndex[k]].lastPosition ,
                            reads[readIndex[k]].firstFragment ,reads[readIndex[k]].lastFragment ,reads[readIndex[k]].flag,
                            readIndex[k] ,j , *idt );
                    j = 1;
                }
                for(  ; j < readExtClusters[k].size() ; j++ )
                    //cout<<"\n2\n";
                    inExtendability( FragsOfRead ,reads[readIndex[k]].firstPosClus2.at(readExtClusters[k].at(j)) ,reads[readIndex[k]].lastPosClus2.at(readExtClusters[k].at(j)) ,
                            reads[readIndex[k]].firstFragClus2.at(readExtClusters[k].at(j)) ,reads[readIndex[k]].lastFragClus2.at(readExtClusters[k].at(j)),reads[readIndex[k]].flagClus2.at(readExtClusters[k].at(j))
                            ,readIndex[k] ,readExtClusters[k].at(j) , *idt );
                //cout<<"\nafter2\n";
            }
            // free old data
            readIndex.clear();
            //readHeader.clear();
            readSeq.clear();
            readExtClusters.clear();
            //until EOF was not reached
        } while( *GLC < 4*this->numberOfReads );

        pFS(" thread "+ convertNumToStr(*idt) +" done.\n");
    }

    if( *readsFileCount > 1){
        string line;
        string line2;
        long long rI = 0 ; // temp read index
        vector<long long> readIndex;
        vector<string> readSeq;
        //vector<string> readHeader;
        vector<vector<int>> readExtClusters; // extend-needed cluster numbers in each read.
        int readLen = 0 ;

        //string outputDir = "/home/ameer/SVAnchoring/SV_out2/";
        //string name2 = outputDir+"read_thread"+convertNumToStr(*idt);
        //ofstream ofstr(name2.c_str());

        //    mutex->lock();
        //    cout<<"welcome to "<< (*idt) <<" thread."<<endl;
        //    mutex->unlock();

        int i;
        bool notEOF = true;
        do {
            i = 0;
            // only 1 concurrent reader
            vector<int> extendingClusters; // extend-needed cluster numbers in each read.
            // ////////////////// **********************************|
            // ////////////////// STEP 1 ---------------------------|
            // ////////////////// Threads compete to read from file |
            // ////////////////// **********************************|
            mainMutex.lock();
            while( i < 4*MAX_BUFFER_SIZE && notEOF ) {

                // ////////////////// ****************************************************************************************************
                // ////////////////// this Processes could be done in second step where there is no lock (if the bottleneck is first step)
                // ////////////////// just pushback the read index and then ...

                if(!getline(*file,line) ){
                    notEOF = false;
                    break;
                }

                if( *GLC % 4 == 0 ){
                    rI = *GLC / 4;
                    extendingClusters.clear();
                    if( !checkBinGenome( min(reads[rI].firstPosition ,reads[rI].lastPosition)+chunkSize  , max(reads[rI].firstPosition ,reads[rI].lastPosition),true ) )
                        extendingClusters.push_back(0);
                    //cout<<"\nafter0\n";
                    for(int j = 1 ; j< reads[rI].firstFragClus2.size() ; j++)
                        if(!checkBinGenome( min(reads[rI].firstPosClus2.at(j),reads[rI].lastPosClus2.at(j))+chunkSize, max(reads[rI].firstPosClus2.at(j),reads[rI].lastPosClus2.at(j)) ,true ) )
                            extendingClusters.push_back(j);

                    if( !getline(*file,line2)){
                        notEOF = false;
                        cout<<"\n WRONG readFILE for inCluster Extension";
                        break;
                    }
                    if(extendingClusters.size()){
                        readIndex.push_back(rI);
                        //readHeader.push_back(line);
                        readSeq.push_back(line2);
                        readExtClusters.push_back(extendingClusters);
                    }
                    (*GLC)++;
                    i++;
                }

                (*GLC)++;
                i++;

            }
            mainMutex.unlock();
            // ////////////////// **********************************
            // ////////////////// STEP 2 ---------------------------
            // ////////////////// Threads seperated processing
            // ////////////////// **********************************
            // this process may takes no time in comparison with STEP1 because maybe there is no read to extend inClusterly,
            // so we cant manage STEP1 the way to collect MAX_BUFFER_SIZE process needed reads instead of MBS (un)/normal reads.
            // execute work for every line
            int completeFrags = 0;
            vector<string> FragsOfRead;
            for (int k = 0; k < readIndex.size(); k++) {
                // work with readSeq and readHeader and readExtClusters
                readLen = (int) readSeq[k].size();
                FragsOfRead.clear();
                completeFrags = floor((readLen-1)/this->chunkSize);
                for( int j = 0; j < completeFrags ; j++ )
                    FragsOfRead.push_back( readSeq[k].substr(  (j*this->chunkSize)  , this->chunkSize ));
                FragsOfRead.push_back( readSeq[k].substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
                int j = 0;
                if( readExtClusters[k].at(0) == 0 ){
                    inExtendability( FragsOfRead ,reads[readIndex[k]].firstPosition ,reads[readIndex[k]].lastPosition ,
                            reads[readIndex[k]].firstFragment ,reads[readIndex[k]].lastFragment ,reads[readIndex[k]].flag,
                            readIndex[k] ,j , *idt );
                    j = 1;
                }
                for(  ; j < readExtClusters[k].size() ; j++ )
                    //cout<<"\n2\n";
                    inExtendability( FragsOfRead ,reads[readIndex[k]].firstPosClus2.at(readExtClusters[k].at(j)) ,reads[readIndex[k]].lastPosClus2.at(readExtClusters[k].at(j)) ,
                            reads[readIndex[k]].firstFragClus2.at(readExtClusters[k].at(j)) ,reads[readIndex[k]].lastFragClus2.at(readExtClusters[k].at(j)),reads[readIndex[k]].flagClus2.at(readExtClusters[k].at(j))
                            ,readIndex[k] ,readExtClusters[k].at(j) , *idt );
                //cout<<"\nafter2\n";
            }
            // free old data
            readIndex.clear();
            //readHeader.clear();
            readSeq.clear();
            readExtClusters.clear();
            //until EOF was not reached
        } while( notEOF );
        pFS(" thread "+ convertNumToStr(*idt) +" done.\n");
    }
}
void SVAnchor::inClusterExtension(int readsFileCount=1){
    // Fastening: based on read coverage we can divide reads and inClusterExtending them in loops, this way in latter loops more reads will be skipped

    if(version == 1){
        string name =  this->readName;
        ifstream ifstr(name.c_str());
        string line;
        string previous_line;
        string ReadSeq;
        vector<string> FragsOfRead;
        int readLen = 0;
        int completeFrags = 0;
        long long lineCounter = 0;

        vector<int> extendingClusters; // extend-needed cluster numbers in each read.
        long long cnttemp = 0;
        for (long long i = 0; i < this->numberOfReads; i++) {
            extendingClusters.clear();
            //cout<<"\n0\n";
            //if( !checkBinGenome( reads[i].firstPosition ,reads[i].lastPosition  ,true ) )
            if( !checkBinGenome( min(reads[i].firstPosition ,reads[i].lastPosition)+chunkSize  , max(reads[i].firstPosition ,reads[i].lastPosition),true ) )
                extendingClusters.push_back(0);
            //cout<<"\nafter0\n";
            for(int j = 1 ; j< reads[i].firstFragClus2.size() ; j++)
                if(!checkBinGenome( min(reads[i].firstPosClus2.at(j),reads[i].lastPosClus2.at(j))+chunkSize, max(reads[i].firstPosClus2.at(j),reads[i].lastPosClus2.at(j)) ,true ) )
                    extendingClusters.push_back(j);
            //cout<<"\nafter02\n";
            //just read the readSeq for skipping or proccessing
            // // // // // // // // // // // // // // // // reading read info for skipping or proccessing
            while (getline(ifstr, line)) {
                if ( lineCounter % 4 == 0 /*&& (lineCounter/4)+1 == *ind*/ ){
                    previous_line = line;
                    getline(ifstr,ReadSeq);
                    lineCounter+=2;
                    break;
                }
                lineCounter++;
            }
            if( extendingClusters.size() ){
                cnttemp++;
                //if(extendingClusters.size() > 1)
                //cout<<"\n extension in this read for "<<extendingClusters.size()<<" Clusters, "<<cnttemp<<"\n";
                // /////////////////////////////////////////////////////
                // /////////////////////////////////////////////////////
                // /////////////////////////////////////////////////////
                if(true/*cnttemp < 20*/)
                {
                    readLen = (int) ReadSeq.size();
                    FragsOfRead.clear();
                    completeFrags = floor((readLen-1)/this->chunkSize);
                    for( int k=0; k < completeFrags ; k++ )
                        FragsOfRead.push_back( ReadSeq.substr(  (k*this->chunkSize)  , this->chunkSize ));
                    //cout<<((readLen-1)/this->chunkSize)*this->chunkSize <<endl;
                    FragsOfRead.push_back( ReadSeq.substr(((readLen-1)/this->chunkSize)*this->chunkSize , readLen-((readLen-1)/this->chunkSize)*this->chunkSize) );
                    //cout<<"\nafterreadret\n";
                    int j = 0;
                    if( extendingClusters.at(0) == 0 ){
                        //cout<<"\n1\n";
                        inExtendability( FragsOfRead ,reads[i].firstPosition ,reads[i].lastPosition ,reads[i].firstFragment ,reads[i].lastFragment ,reads[i].flag,i ,j ,0 );
                        //  cout<<"\nafter1\n";
                        j = 1;
                    }
                    for(  ; j < extendingClusters.size() ; j++ ){ // `on't Enter if extendingClusters.size() == 0
                        //cout<<"\n2\n";
                        inExtendability( FragsOfRead ,reads[i].firstPosClus2.at(extendingClusters.at(j)) ,reads[i].lastPosClus2.at(extendingClusters.at(j)) ,
                                         reads[i].firstFragClus2.at(extendingClusters.at(j)) ,reads[i].lastFragClus2.at(extendingClusters.at(j)),reads[i].flagClus2.at(extendingClusters.at(j))
                                         ,i ,extendingClusters.at(j), 0 );
                        //cout<<"\nafter2\n";
                    }
                }
            }
        }
        // /////////////////////////////////////////////////////
        // /////////////////////////////////////////////////////
        // /////////////////////////////////////////////////////
        cout<<"\n Number of bases Remained after inCluster Extension in BG(s) is :"<<endl;
        if(BGLR)
            cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
        else
            cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;
    }
    if(version == 2 && readsFileCount == 1){

        long long globalLineCounter;
        globalLineCounter = 0;
        // open file
        string name =  this->readName;
        ifstream ifstr(name.c_str());
        // mutex for synchronization
        // maybe a fair-lock would be a better solution
        // create threads and start them with thread_exec(&ifstr, &mutex, &global_line_counter);
        std::vector<std::thread> threads(numberOfThread);
        int *ids;
        ids = new int[numberOfThread];
        int rfc = readsFileCount ;
        for (int i = 0; i < numberOfThread; i++) {
            ids[i] = i;
            threads[i] = std::thread(&SVAnchor::inClusterExtThread, this, &ifstr, &globalLineCounter , &ids[i], &rfc);
        }
        // wait until all threads have finished
        for (int i = 0; i < numberOfThread; i++) {
            threads[i].join();
        }
        delete ids;
        ifstr.close();


        cout<<"\n Number of bases Remained after inCluster Extension in BG(s) is :"<<endl;
        if(BGLR)
            cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
        else
            cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;

    }
    if(version == 2 && readsFileCount > 1){

        long long globalLineCounter;
        globalLineCounter = 0;
        // open file
        string name;
        int rfc = readsFileCount ;
        for(int f = 0 ; f < readsFileCount ; f++){
            name =  this->readName +"_"+convertNumToStr(f)+".fq";
            ifstream ifstr(name.c_str());
            // mutex for synchronization
            // maybe a fair-lock would be a better solution
            // create threads and start them with thread_exec(&ifstr, &mutex, &global_line_counter);
            std::vector<std::thread> threads(numberOfThread);
            int *ids;
            ids = new int[numberOfThread];
            for (int i = 0; i < numberOfThread; i++) {
                ids[i] = i;
                threads[i] = std::thread(&SVAnchor::inClusterExtThread, this, &ifstr, &globalLineCounter , &ids[i], &rfc);
            }
            // wait until all threads have finished
            for (int i = 0; i < numberOfThread; i++) {
                threads[i].join();
            }
            delete ids;
            ifstr.close();
        }
        cout<<"\n Number of bases Remained after inCluster Extension in BG(s) is :"<<endl;
        if(BGLR)
            cout<<"   in BGLeft: "<<basesRemained(false)<<" - and in BGRight:"<<basesRemained(true)<<endl;
        else
            cout<<"   in BinaryGenome: "<<basesRemained(true)<<endl;
    }
}
inline bool SVAnchor::findNewClusters(vector<string> FragsOfRead, int firstFrag ,long long firstPos, int lastFrag,long long lastPos,int flag, long long readIndex){
    // check if you can find new clusters
    // frag1+1 if Inversion
    int shiftf1, shiftf2;
    int shiftl1, shiftl2;
    // ifAlignable(long long startLoci, string readString, int flag, bool isToRight, int numGapCount, bool useBinGenome, int *shift){
    // check if it is inversion (most probable event here)
    ////
    ////
    // ////////////////////////////////////
    // ////////////////////////////////////
    // ////////////////////////////////////
    // ////////////////////////////////////

    bool isToRight = false;
    int jumpsCount = 0;

    // ////////////////////////////////////
    // ////////////////////////////////////
    // ////////////////////////////////////



    bool first1 = ifAlignable(lastPos - ( chunkSize*(flag==0?1:-1) ) , FragsOfRead.at(firstFrag+1) , (flag==0?16:0), isToRight, jumpsCount, false, &shiftf1);
    bool first2 = ifAlignable(lastPos - (2*chunkSize*(flag==0?1:-1)) , FragsOfRead.at(firstFrag+2) , (flag==0?16:0), isToRight, jumpsCount, false, &shiftf2);

    bool last1 = ifAlignable(lastPos + ( chunkSize*(flag==0?1:-1) ) , FragsOfRead.at(lastFrag - 1) , (flag==0?16:0), isToRight, jumpsCount, false, &shiftl1);
    bool last2 = ifAlignable(lastPos + (2*chunkSize*(flag==0?1:-1)) , FragsOfRead.at(lastFrag - 2) , (flag==0?16:0), isToRight, jumpsCount, false, &shiftl2);


    // if not a full inversion so we have to do Local Alignments
    // or we can first mix this two parts together but if repeat?

}

inline bool SVAnchor::checkInBetween(vector<string> FragsOfRead, int firstFrag ,long long firstPos, int lastFrag,long long lastPos,int flag, long long readIndex){
    if( lastFrag - firstFrag < 2 && (lastPos-firstPos)*(flag==1?1:-1) > this->chunkSize + this->numGap){
    //if( abs(lastFrag - firstFrag) < 2 && abs((lastPos-firstPos)*(flag==1?1:-1)) > this->chunkSize + this->numGap){
        // Bridged Deletion
        // do nothing more
        return true;
    }
    if( lastFrag - firstFrag > 1 && (lastPos-firstPos)*(flag==1?1:-1) <= (this->chunkSize + this->numGap) ){
    //if( abs(lastFrag - firstFrag) > 1 && abs((lastPos-firstPos)*(flag==1?1:-1)) < 2 * this->chunkSize + this->numGap){
        // Bridged Insertion
        // need to choose between novel and mobile for further peak calling and type detecion
        // so add it to a vetor of Bridged unknown sequences to check if it's novel or not (treating as a read)
        return true;
    }
    return findNewClusters(FragsOfRead ,firstFrag ,firstPos ,lastFrag ,lastPos ,flag ,readIndex );

    // try to check if there is at least one more chunk corresponding to this region and flag, then it is noise,
    // otherwise maybeSV
}

inline void SVAnchor::inExtendability(vector<string> FragsOfRead ,long long firstPos ,long long lastPos ,int firstFrag ,int lastFrag,int flag, long long index , int clustCode, int threadNum){

    if(firstFrag == lastFrag){
        //cout<<"__\n baha: "<<firstFrag<<" "<<lastFrag<<" "<<clustCode<<endl;
        return ;
    }
    if( firstFrag > lastFrag){
        long long temp = firstPos;
        firstPos = lastPos;
        lastPos = temp;
        int temp2 = firstFrag;
        firstFrag = lastFrag;
        lastFrag = temp2;
    }
    int counter1 = 0;
    int counter2 = 0;
    int jumpsCount = 0;
    char successFlag = 0;
    long long startLoci = firstPos-1;
    prevShiftPos.at(threadNum) = 0 ;
    shiftPos.at(threadNum) = 0 ;
    int firstShiftPos = 0;
    int endtemp = abs(firstFrag-lastFrag)-1;
    bool notReverse = (flag == 0 ? true : false );
    long long lastStartLoci;

    // first to last
    //cout<<"inExtendability 001"<<endl;
    //cout<<firstFrag<<" "<<lastFrag<<endl;
    do{
        //cout<<"in 001 :"<<firstFrag + (counter1+1) -1<<endl;
        //startLoci += (direction*(this->chunkSize + shiftPos.at(threadNum)))*( extendFrom.Flag==0 ? 1 : -1 );//isForward(extendFrom.Flag);
        startLoci += ((this->chunkSize + shiftPos.at(threadNum))*(notReverse ? 1 : -1 ));
        successFlag =  ifAlignable(startLoci ,FragsOfRead.at(firstFrag + (counter1+1) -1 ) ,flag ,0  ,notReverse ,jumpsCount ,false );
        counter1 += ( successFlag ? 1:0 );
        jumpsCount += ( successFlag == 1 ? 1:0);
    }while( successFlag  && counter1 < endtemp );

    lastStartLoci = startLoci;
    if( !successFlag )
        lastStartLoci -= ((this->chunkSize + shiftPos.at(threadNum))*(notReverse ? 1 : -1 ));
    //cout<<"inExtendability 002"<<endl;
    // last to first
    firstShiftPos = shiftPos.at(threadNum);
    prevShiftPos.at(threadNum) = 0 ;
    shiftPos.at(threadNum) = 0;
    successFlag = 0;
    jumpsCount = 0;
    //cout<<"inExtendability 003"<<endl;
    startLoci = lastPos - 1 ;
    //endtemp = abs(firstFrag-lastFrag)-1-counter1;
    do{
        startLoci += (-1)*((this->chunkSize + shiftPos.at(threadNum))*(notReverse ? 1 : -1 ));//isForward(extendFrom.Flag);
        successFlag = ifAlignable( startLoci ,FragsOfRead.at(lastFrag + (-1)*(counter2+1) -1 ) ,flag ,0 ,!notReverse ,jumpsCount ,false  ) ;
        counter2 += (successFlag ? 1:0 );
        jumpsCount += ( successFlag == 1 ? 1:0);
    }while( successFlag && counter2 < endtemp );
    if( !successFlag )
        startLoci -= (-1)*((this->chunkSize + shiftPos.at(threadNum))*(notReverse ? 1 : -1 ));//isForward(extendFrom.Flag);
    int value = abs(firstFrag-lastFrag) - 1 - counter1 - counter2;
    if( value >= 0 &&
            !((value == 0) && abs(lastStartLoci-startLoci) <= chunkSize + numGap ) ){
        bool maybeSV = false;
        long long newFirstPos = firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);
        long long newLastPos = lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum));
        maybeSV = true;
        if(false)
            maybeSV = checkInBetween(FragsOfRead, firstFrag+counter1,newFirstPos,lastFrag-counter2,newLastPos,flag,index);
        // new Clusters
        //-/////////////////////////////////////////////////////////////////////////////////////////////////////
        //-////////////////////////////   BinGenome Changes Needed yet  ////////////////////////////////////////
        //-/////////////////////////////////////////////////////////////////////////////////////////////////////
        if(maybeSV){
            changed++;
            if(clustCode == 0){
                if(reads[index].lastFragment == lastFrag ){
                    cout<<" change1 in read Clusters with index: "<<reads[index].index<<endl;
                    inClustChangesMutex.lock();

                    reads[index].lastFragment = firstFrag+counter1;
                    //reads[index].lastPosition = firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);//newFirstPos
                    reads[index].lastPosition = newFirstPos;

                    reads[index].flagClus2.push_back(flag);
                    reads[index].firstFragClus2.push_back(lastFrag-counter2);
                    //reads[index].firstPosClus2.push_back(lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum)));
                    reads[index].firstPosClus2.push_back(newLastPos);
                    reads[index].lastFragClus2.push_back(lastFrag);
                    reads[index].lastPosClus2.push_back(lastPos);

                    inClustChangesMutex.unlock();

                    firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);
                    firstFrag+counter1;

                    lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum));
                    lastFrag-counter2;


                    // anyClustersInBetween();

                }else if(reads[index].lastFragment == firstFrag ){
                    cout<<" change2 in read Clusters with index: "<<reads[index].index<<endl;

                    inClustChangesMutex.lock();

                    reads[index].firstFragment = firstFrag;
                    reads[index].firstPosition = firstPos;
                    reads[index].lastFragment = firstFrag+counter1;
                    //reads[index].lastPosition = firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);// newFirstPos
                    reads[index].lastPosition = newFirstPos;
                    reads[index].flagClus2.push_back(flag);
                    reads[index].firstFragClus2.push_back(lastFrag-counter2);
                    //reads[index].firstPosClus2.push_back(lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum)));
                    reads[index].firstPosClus2.push_back(newLastPos);
                    reads[index].lastFragClus2.push_back(lastFrag);
                    reads[index].lastPosClus2.push_back(lastPos);

                    inClustChangesMutex.unlock();

                }else cout<<"\n Error in inExtend01\n";
            }else{
                if(reads[index].lastFragClus2.at(clustCode) == lastFrag ){
                    cout<<" change3 in read Clusters with index: "<<reads[index].index<<endl;

                    inClustChangesMutex.lock();

                    reads[index].lastFragClus2.at(clustCode) = firstFrag+counter1;
                    //reads[index].lastPosClus2.at(clustCode) = firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);//newFirstPos
                    reads[index].lastPosClus2.at(clustCode) = newFirstPos;

                    reads[index].flagClus2.push_back(flag);
                    reads[index].firstFragClus2.push_back(lastFrag-counter2);
                    //reads[index].firstPosClus2.push_back(lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum)));
                    reads[index].firstPosClus2.push_back(newLastPos);
                    reads[index].lastFragClus2.push_back(lastFrag);
                    reads[index].lastPosClus2.push_back(lastPos);

                    inClustChangesMutex.unlock();

                }else if(reads[index].lastFragClus2.at(clustCode) == firstFrag ){
                    cout<<" change4 in read Clusters with index: "<<reads[index].index<<endl;

                    inClustChangesMutex.lock();

                    reads[index].firstFragClus2.at(clustCode)  = firstFrag;
                    reads[index].firstPosClus2.at(clustCode)  = firstPos;
                    reads[index].lastFragClus2.at(clustCode)  = firstFrag+counter1;
                    //reads[index].lastPosClus2.at(clustCode)  = firstPos + (flag==0?1:-1)*((counter1*this->chunkSize) + firstShiftPos);//newFirstPos
                    reads[index].lastPosClus2.at(clustCode)  = newFirstPos;

                    reads[index].flagClus2.push_back(flag);
                    reads[index].firstFragClus2.push_back(lastFrag-counter2);
                    //reads[index].firstPosClus2.push_back(lastPos - (flag==0?1:-1)*((counter2*this->chunkSize) + shiftPos.at(threadNum)));//newLastPos
                    reads[index].firstPosClus2.push_back(newLastPos);//newLastPos
                    reads[index].lastFragClus2.push_back(lastFrag);
                    reads[index].lastPosClus2.push_back(lastPos);

                    inClustChangesMutex.unlock();

                }else cout<<"\n Error in inExtend02";
            }
        }
    }else { // comment this block to see changes in speed!!! and some changes in accuracy too!!!

        skipped++;

        // V4 - 1 BG
        if(!BGLR && abs(firstFrag-lastFrag) > 3 )
        changeBinGenome(firstPos + (notReverse ? 2*this->chunkSize : -1*this->chunkSize)
                        ,lastPos + (notReverse ? -1*this->chunkSize : 2*this->chunkSize)
                        ,true);
        // V4 - 2 BG
        if(BGLR && abs(firstFrag-lastFrag) > 2 ){

            changeBinGenome(firstPos + (notReverse ? 2*this->chunkSize : 0)
                            ,lastPos + (notReverse ? 0 : 2*this->chunkSize)
                            ,true);
            changeBinGenome(firstPos + (notReverse ? 1*this->chunkSize : -1*this->chunkSize)
                            ,lastPos + (notReverse ? -1*this->chunkSize : 1*this->chunkSize)
                            ,false);
        }
    }
}
void SVAnchor::writeBinGenome(){
    string name = outputDir+"BinGenomeRight.txt";
    string name2 = outputDir+"BinGenomeLeft.txt";

    ofstream ofstr1(name.c_str());
    for(long long i=0 ; i < bGRight.size() ;i++){
        ofstr1<< bGRight[i];
    }
    ofstr1.close();
    ofstream ofstr2(name2.c_str());
    for(long long i=0 ; i < bGLeft.size() ;i++){
        ofstr2<< bGLeft[i];
    }
    ofstr2.close();
}

void SVAnchor::writeReads(string alter) {

    string name = outputDir+alter+"_reads.txt";

    ofstream ofstr(name.c_str());
    // //////////////////
    // NORMAL READS :
    // //////////////////
    for (long long i = 0; i < this->numberOfReads; i++) {
        //-/ Ameerosein edits
        //ofstr << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" << this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition << "\t" << this->reads[i].flag << "\t" << this->reads[i].cigar << "\t";
        ofstr << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" <<
                 this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition <<"\t"<<reads[i].flag<<
                 "\t"<< (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].firstFragClus2.at(1):0 ) << "\t" << (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].lastFragClus2.at(1):0)<<
                 "\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].firstPosClus2.at(1):0)<<"\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].lastPosClus2.at(1):0)<<
                 "\t"<< (this->reads[i].firstFragClus2.size()>1?this->reads[i].flagClus2.at(1):0); //<< "\t" << this-   >reads[i].cigar << "\t";
        //this->reads[i].origin
        //this->reads[i].iRead
        //-/ end Ameerosein edits
#ifdef runInRealMode
#else
 //        for (int j = 0; j < 5; j++) {
 //            ofstr << this->reads[i].origin[j] << "\t";
 //        }
#endif
        ofstr << endl;
    }
    ofstr.close();
    // //////////////////
    // UN-NORMAL READS :
    // //////////////////
    string name2 = outputDir+alter+"_reads_unormal.txt";
    ofstream ofstr2(name2.c_str());

    for (long long i = 0; i < this->numberOfReads; i++) {
        //-/ Ameerosein edits
        //ofstr << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" << this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition << "\t" << this->reads[i].flag << "\t" << this->reads[i].cigar << "\t";
        if((abs(reads[i].firstFragment-reads[i].lastFragment)+1)*this->chunkSize <= this->reads[i].length - this->chunkSize ){
            ofstr2 << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" <<
                     this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition <<"\t"<<reads[i].flag<<
                     "\t"<< (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].firstFragClus2.at(1):0 ) << "\t" << (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].lastFragClus2.at(1):0)<<
                     "\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].firstPosClus2.at(1):0)<<"\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].lastPosClus2.at(1):0)<<
                     "\t"<< (this->reads[i].firstFragClus2.size()>1?this->reads[i].flagClus2.at(1):0); //<< "\t" << this->reads[i].cigar << "\t";
            //this->reads[i].origin
            //this->reads[i].iRead
            //-/ end Ameerosein edits
#ifdef runInRealMode
#else
            //        for (int j = 0; j < 5; j++) {
            //            ofstr << this->reads[i].origin[j] << "\t";
            //        }
#endif
            ofstr2 << endl;
        }
    }
    ofstr2.close();

    if(alter=="beforExt"){
        //pFS("\n ENTERED !!!!!!!!!!!!!!!!!!!\n");
        string name3 = outputDir+alter+"_Extindices.txt";
        ofstream ofstr3(name3.c_str());
        ofstr3<<allExtIndices.size()<<endl;
        for(vector<int>::iterator it = allExtIndices.begin() ; it != allExtIndices.end();++it){
            ofstr3<<*it<<endl;
        }
        ofstr3.close();
    }
    return;
}

void SVAnchor::readSingleUniques(){
    singleUniqes.resize( this->numberOfReads );
    string line = outputDir+"/SingleUniques.txt";
    ifstream ifstr(line.c_str());
    vector<string> tokens;
    while(getline(ifstr,line)){
        std::istringstream iss(line);
        std::string token;
        tokens.clear();
        while(std::getline(iss, token, '\t'))   // but we can specify a different one
            tokens.push_back(token);
        singleUniqes[ stoll(tokens[0])-1 ].push_back( singleBowtied(stoll(tokens[0]) , stoi(tokens[1]) , stoll(tokens[2]) , stoi(tokens[3]) ) );
    }
}

void SVAnchor::pruneSingleUniques(){
    long long readIndex;
    int fragNum;
    long long snap;
    int flag;

    vector<vector<singleBowtied>> newSingleUniques;
    newSingleUniques.resize(singleUniqes.size());
    bool seen = false;
    for( long long i = 0 ; i < this->numberOfReads ; i++ ){
        for( int j = 0 ; j < singleUniqes[ i ].size() ; j++ ){
            seen = false;
            if( singleUniqes[ i ][ j ].fragNum == reads[ i ].firstFragment )
                seen = true;
            else{
                if( singleUniqes[ i ][ j ].fragNum == reads[ i ].lastFragment )
                    seen = true;
                else{
                    for( int k = 0 ; k < reads[ i ].firstFragClus2.size() ; k++ )
                        if( singleUniqes[ i ][ j ].fragNum == reads[ i ].firstFragClus2[ k ] )
                            seen = true;
                        else{
                            for( int k = 0 ; k < reads[ i ].lastPosClus2.size() ; k++ )
                                if( singleUniqes[ i ][ j ].fragNum == reads[ i ].lastPosClus2[ k ] )
                                    seen = true;
                        }
                }
            }
            if( !seen )
                newSingleUniques[i].push_back( singleUniqes[ i ][ j ] );
        }
    }
    singleUniqes = newSingleUniques ;
}
void SVAnchor::writeSingleUniques(){
    string file = outputDir+"/SingleUniques.txt";
    ofstream ofstr(file.c_str());
    for( long long i = 0 ; i < this->numberOfReads ; i++ ){
        ofstr<<" readIndex :"<<i+1<<endl;
        for( int j = 0 ; j < singleUniqes[ i ].size() ; j++ ){
            ofstr<<singleUniqes[ i ][ j ].fragNum<<" - "<<singleUniqes[ i ][ j ].snap<<" - "<<singleUniqes[ i ][ j ].flag<<endl;
        }
    }

//long long readIndex;
//int fragNum;
//long long snap;
//int flag;
}
int SVAnchor::findCharacterPerRow(string genomeName) {
    int characterPerRow;
    ifstream ifstr(genomeName.c_str());
    string line;
    getline(ifstr, line);
    getline(ifstr, line);
    characterPerRow = (int)line.size();
    ifstr.close();
    return characterPerRow;
}


long long SVAnchor::findNumberOfRead(string fileName) {
    long long numberOfReads = 0;
    ifstream ifstr(fileName.c_str());
    string line;
    while (getline(ifstr, line)) {
         ++numberOfReads;
    }
    ifstr.close();
    return numberOfReads;
}




// // Commented codes:


/*=======================================================================*/
/*============================   ASSIGNMENT   ===========================*/
/*=======================================================================*/
//    if (runAssignment) {
//        pFS(get_current_time());pFS(" ASSIGNMENT step starts ");pFS("\n");

//        this->readReads("an");
//        this->readDepth("an");
//        this->printStatus();

//        if (this->isSlideNotShift) {
//            pFS(get_current_time());pFS(" Assigning with sliding starts");pFS("\n");
//            Assignment *assignment = new Assignment(this->reads, &this->depth, &this->genes, &this->readsNext, this->dashV, this->numberOfThread, this->chunkSize, "as", 0, this->readLength, this->characterPerRow, this->chunkSize-this->assignmentShift, this->genomeName, this->indexName, this->numGap, this->outputDir);
//            assignment->runForOnGenes();
//            assignment->createNewRead(0, "rb");
//            delete assignment;
//        } else {
//            for (int i = 0; i < floor(this->chunkSize/this->assignmentShift)-1; i++) {
//                pFS(get_current_time());pFS(" Assigning with d = ");pFS(convertNumToStr(this->assignmentShift*i));pFS(" starts");pFS("\n");
//                Assignment *assignment = new Assignment(this->reads, &this->depth, &this->genes, &this->readsNext, this->dashV, this->numberOfThread, this->chunkSize, "as", this->assignmentShift*i, this->readLength-this->assignmentShift*i, this->characterPerRow, 0, this->genomeName, this->indexName, this->numGap, this->outputDir);
//                if (i == floor(this->chunkSize/this->assignmentShift)-2) {
//                    assignment->runForOnGenes();
//                    assignment->createNewRead(0, "rb");
//                } else {
//                    assignment->createNewRead(this->assignmentShift*(i+1), "as");
//                }
//                delete assignment;
//            }
//        }


//        this->writeReads("as");
//        this->writeReadsNext("as");
//        this->writeDepth("as");
//        this->writeGenes("as");
//    }

//    if (runBowtie2) {
//        pFS(get_current_time());pFS(" RUNNING BOWTIE2 step starts");pFS("\n");

//        this->readReads("as");
//        this->readDepth("as");
//        this->printStatus();

//        RunBowtie2 *rb2 = new RunBowtie2(this->reads, &this->depth, this->outputDir);
//        rb2->run(0, "rb", this->numberOfThread, 1, 20, "CHR18");
//        delete rb2;

//        /*long long counter1 = 0;
//        int i = 0;
//        int shift = 10;
//        for (i = 0; i < floor(50/shift); i++) {
//            pFS(get_current_time());pFS(" Aligning with d = ");pFS(convertNumToStr(shift*i));pFS(" starts");pFS("\n");
//            Anchoring *anchoring = new Anchoring(this->reads, &this->depth, this->indexName, "rb", true, this->numGap, 50, shift*i, 5, this->dashV, this->numberOfThread, this->outputDir);
//            anchoring->prepareFiles();
//            anchoring->nextStepAlignment();
//            anchoring->findExonIntron();
//            anchoring->findTwoAdjacentExon();
//            anchoring->createNewReads(shift*(i+1));
//            counter1 += anchoring->correspond;
//            pFS(get_current_time());pFS(" Number of anchoring in this step is ");pFS(convertNumToStr(anchoring->correspond));pFS("\n");
//            delete anchoring;
//            this->deleteAllJunk("rb", shift*(i+1));
//        }
//        pFS(get_current_time());pFS(" Total anchoring is ");pFS(convertNumToStr(counter1));pFS("\n");*/

//        this->writeDepth("rb");
//    }

/*=======================================================================*/
/*=========================== LOCAL ALIGNMENT ===========================*/
/*=======================================================================*/
//    if (runLocalAlignment) {
//        pFS(get_current_time());pFS(" LOCAL ALIGNING starts");pFS("\n");

//        this->readReads("");
//        this->readDepth("");
//        this->printStatus();

//        int rc2;
//        int i2;
//        void *status2;

//        int threadsNO2 = 4;
//        pthread_t threads2[threadsNO2];
//        pthread_attr_t attr2;
//        struct wrapAligning td2[threadsNO2];

//        pthread_attr_init(&attr2);
//        pthread_attr_setdetachstate(&attr2, PTHREAD_CREATE_JOINABLE);

//        Aligning *alig = new Aligning(this->reads, 80, this->numberOfReads, "CHR19", 30);
//        alig->breakReadsForThreads(threadsNO2);

//        for (i2 = 0; i2 < threadsNO2; i2++) {
//            td2[i2].t_id = i2;
//            td2[i2].ins = alig;
//            //rc2 = pthread_create(&threads2[i2], NULL, callFuncAligning, (void *)&td2[i2]);
//            if (rc2) {
//                cout << "Error: unable to create thread," << rc2 << endl;
//                exit(-1);
//            }
//        }
//        pthread_attr_destroy(&attr2);
//        for (i2 = 0; i2 < threadsNO2; i2++) {
//            rc2 = pthread_join(threads2[i2], &status2);
//            if (rc2) {
//                cout << "Error: unable to join," << rc2 << endl;
//                exit(-1);
//            }
//        }

//        alig->prepareHeaderFile();
//        alig->combineSamForThreads(threadsNO2, false);
//        delete alig;

//        this->writeReads("");
//        this->printStatus();
//        this->deleteAllJunk("u", 0);
//    }

//void SVAnchor::printStatus() {
//    long long numberOfFiveChunkRemained = 0;
//    long long numberOfCompleted = 0;
//    long long numberOfOneChunkRemained = 0;
//    long long numberOfTwoChunkRemained = 0;
//    long long numberOfThreeChunkRemained = 0;

//    int temp1 = 0;
//    int temp2 = 0;
//    for (long long i = 0; i < this->numberOfReads; i++) {
//        if (this->reads[i].flag == -1) {
//            numberOfFiveChunkRemained++;
//        } else {
//            temp1 = (this->reads[i].firstFragment-1);
//            temp2 = (5-this->reads[i].lastFragment);
//            if (temp1+temp2 == 0) {
//                numberOfCompleted++;
//            } else if (temp1+temp2 == 1) {
//                numberOfOneChunkRemained++;
//            } else if (temp1+temp2 == 2) {
//                numberOfTwoChunkRemained++;
//            } else if (temp1+temp2 == 3) {
//                numberOfThreeChunkRemained++;
//            }
//        }
//    }

//    cout << "Number of unanchored = " << numberOfFiveChunkRemained << endl;
//    cout << "Number of 3 chunk remained = " << numberOfThreeChunkRemained << endl;
//    cout << "Number of 2 chunk remained = " << numberOfTwoChunkRemained << endl;
//    cout << "Number of 1 chunk remained = " << numberOfOneChunkRemained << endl;
//    cout << "Number of 0 chunk remained = " << numberOfCompleted << endl;
//    return;
//}
