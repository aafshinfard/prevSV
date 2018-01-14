///*
// * finalModule.cpp
// *
// *  Created on: Mar 27, 2015
// *      Author: amin
// */



//#include "Assignment.h"

//Assignment::Assignment(iRead* reads, vector<interval>* iDepth, vector<gene>* genes, vector<iReadNext>* readsNext, int V, int myNumberOfThread, int partLength, string firstOutputName, int shift, int readLength, int characterPerRow, int d, string genomeName, string indexAddress, int numGap, string outputDir) {
//    this->iDepth = iDepth;
//    this->reads = reads;
//    this->genes = genes;
//    this->readsNext = readsNext;

//    this->firstOutputName = firstOutputName;
//    this->shift = shift;
//    this->characterPerRow = characterPerRow;
//    this->outputDir = outputDir;

//    long long index;
//    long long number;
//    long long failNumber = 0;
//    long long current = 0,currentFast = 0,fastLineCounter = 0;
//    long long currentNumber,currentFastNumber;
//    long long readCounter = 0;
//    long long length;
//    vector<string> splitted;
//    bool isNew= true;
//    bool null = false;
//    int V_L = 0;
//    int row;
//    int threadCounter = 0;
//    int blockCounter = 0;
//    char sign = '@';
//    std::string tempString, sequencePrime, qualityPrime, Header, mainHeader;
//    std::string seq, line, quality;
//    std::string fastLine;
//    std::string token;
//    std::string command, concatenate, directory;
//    std::string sortedFinalReadsAddress = outputDir+"sortedReads.fq";
//    std::string semiSortedFinalReadsAddress = outputDir+"semiSortedReads.fq";
//    std::string firstReadsName = outputDir+this->firstOutputName+"_reads_d="+this->NumToStr(this->shift)+".fq";
//    std::string ReadsName = outputDir+"reformed_Reads.fq";
//    std::string outputName = outputDir+"StatisticResult.txt";
//    std::string AlignedAddress = outputDir+"Aligned.sam";
//    std::string ReOrderedAlignedAddress = outputDir+"Aligned_reordered.sam";
//    std::string failReads = outputDir+"Assign_Rem_Reads.fq";
//    fastQAddress = outputDir+"fastQAddress.fq";
//    this->genomeName = genomeName;
//    this->indexAddress = indexAddress;
//    this->numberOfThreads = 2;
//    this->myNumberOfThread = myNumberOfThread;
//    this->readLength = readLength;
//    this->partLength = partLength;
//    this->d = d;
//    this->consecutiveThreshold = 0.3;
//    this->error = numGap;
//    this->myDepth = 10;
//    IndelShift = 0.1;
//    errorInIndex = 2000;
//    GAPPEN = -8;
//    MISPEN = -6;
//    MATCH = 2;
//    this->V = V;
//    K_factor=0.5;
//    L_MAX = 2000;
//    L_bowtie = 20;
//    N_bowtie = 1;
//    step = 2;
//    aligner = 1;


//    //cout << "parameters" << endl;
//    //system(("mkdir " + directory + " 2>/dev/null").c_str());
//    //chdir(directory.c_str());
//    readFasta(genomeName, multiFasta);
//    //cout << "Rehearse the fasta file and sort it" << endl;

//    std::ifstream firstReads;
//    std::ofstream OfinalReads(semiSortedFinalReadsAddress.c_str());
//    std::ofstream Passed_reads((outputDir + "Assign_Rem_Reads.fq").c_str());
//    std::ofstream fastQfile(fastQAddress.c_str());
//    std::ofstream output(outputName.c_str());
//    myDepth++;
//    option = 40;
//    // In constant length mode the -L parameter stands for the constant read length
//    // In pacbio mode the -L parameter stands for L_MAX
//    // Defining the files for sam and fastQ
//    filledAll=0;
//    //COMPLETED
//    // Seperation of the Header considering the constant length Mode

//    if(Pacbio){
//        std::ifstream Reads(firstReadsName.c_str());
//        std::ofstream O_Reads(ReadsName.c_str());
//        readCounter=0;
//        while(getline(Reads,token)){
//            if(counter % 4 ==0){
//                Header = token;
//            }
//            if(counter % 4 ==1){
//                seq = token;
//            }
//            if(counter % 4 ==3){
//                readCounter++;
//                O_Reads << "@r" << readCounter << "_" << counter << "_" << seq.length() << "_" << Header << endl;
//                O_Reads << seq << endl;
//                O_Reads << "+" << endl;
//                O_Reads << token << endl;
//            }
//            counter++;
//        }
//        readCounter=0;
//        O_Reads.close();
//        firstReads.open(ReadsName.c_str());

//    }
//    else{
//        //firstReads.open(firstReadsName.c_str());
//        std::ifstream Reads(firstReadsName.c_str());
//        std::ofstream O_Reads(ReadsName.c_str());
//        readCounter=0;
//        while(getline(Reads,token)){
//            if(counter % 4 ==0){
//                Header = token;
//            }
//            if(counter % 4 ==1){
//                seq = token;
//            }
//            if(counter % 4 ==3){
//                readCounter++;
//                O_Reads << Header << "_" << seq.length() << endl;
//                O_Reads << seq << endl;
//                O_Reads << "+" << endl;
//                O_Reads << token << endl;
//            }
//            counter++;
//        }
//        readCounter=0;
//        O_Reads.close();
//        firstReads.open(ReadsName.c_str());

//    }
//    counter=0;
//    if(!firstReads.is_open()){
//        cerr << "the first read is not open" << endl;
//    }
//    if(!fastQfile.is_open()){
//        cerr << "fastQfile is not open" << endl;
//    }
//    while (getline(firstReads,token)){
//        if (counter % 4 == 0) {
//            Header = token;
//            //for constant length mode
//            Split_Header(splitted, Header, mainHeader, index, number, V_L, row, false, false, true, Pacbio,false);
//            //update the numberOfTables to calculate the indexMaker
//            numberOfTables = checkNumOfTables(partLength,V_L);
//        }
//        if (counter % 4 == 1) {
//            seq = token;
//            tempString = L_MAX_Gate(L_MAX,V_L,seq);
//            if(tempString.size()!=V_L){
//                reform_Header(splitted,Header,index,number,L_MAX, true,Pacbio,false,true,0);
//                numberOfTables = checkNumOfTables(partLength,L_MAX);
//            }
//            null = false;
//            size_t N = std::count(seq.begin(), seq.end(), 'N');
//            size_t n = std::count(seq.begin(), seq.end(), 'n');
//            if(n+N > 0.15*seq.length()){
//                null = true;
//                failNumber++;
//            }
//        }
//        if (counter % 4 == 3) {
//            if(!null){
//                reform_Header(splitted,Header,index,number,V_L,true,Pacbio, true,false,0);
//                OfinalReads << "@r" << Header <<"\t" << seq << "\t" << token << endl;
//                filledAll++;
//            }
//        }
//        counter++;
//    }
//    system_sort(semiSortedFinalReadsAddress, sortedFinalReadsAddress);
//    system_remove(semiSortedFinalReadsAddress);
//    std::ifstream finalReads(sortedFinalReadsAddress.c_str());
//    //cout << "First adjustment of input is done" << endl;
//    //cout << "All of the reads:" << "\t" <<filledAll << endl;
//    while (getline(finalReads,token)){
//        Split_Token(splitted, token, Header, seq, quality);
//        Split_Header(splitted, Header,mainHeader, index, number, V_L, row, true, false, true, Pacbio, false);
//        numberOfTables = checkNumOfTables(partLength,V_L);
//        if(numberOfTables>1){
//            for (int i=0;i<numberOfTables;i++){
//                sequencePrime = seq.substr(i * (partLength-d),
//                                           partLength);
//                qualityPrime = quality.substr(i * (partLength-d),
//                                              partLength);
//                reform_Header(splitted,Header,index,number,V_L,true,Pacbio, true,false,i);
//                fastQfile <<"@r" << Header << endl;
//                fastQfile << sequencePrime << endl;
//                fastQfile << "+"<<endl;
//                fastQfile << qualityPrime << endl;

//            }
//        }
//        else{
//            fastQfile <<"@r" << Header << endl;
//            fastQfile << sequencePrime << endl;
//            fastQfile << "+"<<endl;
//            fastQfile << qualityPrime << endl;
//        }
//    }
//    /*indel_stat.resize(NOISE_THREADS);

//     for (int i = 0; i<indel_stat.size() ; i++)
//     ifstr2 >> indel_stat[i];

//     for (int i=0;i<NOISE_THREADS;i++)
//     cout << indel_stat[i] << endl;

//     frag = (long long)floor(GENE_LEN/NOISE_THREADS);*/
//    //cout << "Arrangement of the reads is ended up" << endl;
//    //cout << "Creating the fragments in fastQ order to be aligned" << endl;
//    //cout << "Reads with Not Desired Qualification " << failNumber << endl;
//    finalReads.clear();
//    finalReads.seekg(0, ios::beg);
//    fastQfile.close();
//    // Reading the file contains noise of the reads
//    // splitting the fastqfile into smaller files for threads
//    // contain the read index and its sequence
//    length = Calculate(filledAll, numberOfThreads);
//    if(length==0)
//        numberOfThreads=1;
//    errorCollector.resize(numberOfThreads);
//    threadVectors.resize(numberOfThreads);
//    threadReads.resize(numberOfThreads);
//    myOSams = new ofstream[numberOfThreads];
//    myOFasts = new ofstream[numberOfThreads];
//    O_Strings = new ofstream[numberOfThreads];
//    myFasts = new ifstream[numberOfThreads];
//    mySams = new ifstream[numberOfThreads];
//    for (int thread = 0; thread < numberOfThreads; thread++) {
//        ostringstream filename,fastname,finalStrName;
//        filename << outputDir + "sam"<< thread + 1 << ".txt";
//        fastname << outputDir + "fast"<< thread + 1<< ".txt";
//        finalStrName << outputDir + "strings" << thread + 1 << ".txt";
//        myOSams[thread].open(filename.str().c_str());
//        myOFasts[thread].open(fastname.str().c_str());
//        O_Strings[thread].open(finalStrName.str().c_str());
//    }
//    // Creating thread fastQ files
//    //cout << "Assigning every Thread its quota...." << endl;
//    if(length!=0)
//        while (getline(finalReads,token) && blockCounter<numberOfThreads){
//            if((readCounter==length) && blockCounter!=(numberOfThreads-1)){
//                blockCounter++;
//                readCounter=0;
//            }
//            //Saving all of the properties of the read
//            Split_Token(splitted, token, Header, seq, quality);
//            Split_Header(splitted, Header, mainHeader, index, number, V_L, row, true, false, true, Pacbio, false);
//            if(!myOFasts[blockCounter].is_open())
//                cerr << "is not open" << endl;
//            else{
//                myOFasts[blockCounter] << "@r" << Header << "\t" << seq << "\t" << quality << endl;
//                readCounter++;
//            }
//        }
//    else{
//        //blockCounter=0;
//        while (getline(finalReads,token)){
//            //Saving all of the properties of the read
//            Split_Token(splitted, token, Header, seq, quality);
//            Split_Header(splitted, Header, mainHeader, index, number, V_L, row, true, false, true, Pacbio, false);
//            if(!myOFasts[0].is_open())
//                cerr << "is not open" << endl;
//            else{
//                myOFasts[0] << "@r" << Header << "\t" << seq << "\t" << quality << endl;
//                readCounter++;
//            }
//        }
//    }
//    for (int thread = 0; thread < numberOfThreads; thread++) {
//        ostringstream filename,fastname;
//        filename << outputDir + "sam"<< thread + 1 << ".txt";
//        fastname << outputDir + "fast"<< thread + 1<< ".txt";
//        mySams[thread].open(filename.str().c_str());
//        myFasts[thread].open(fastname.str().c_str());

//    }
//    //cout << "Ready for Aligning" << endl;
//    align_Command(aligner, indexAddress, fastQAddress, V, myDepth, numberOfThreads, AlignedAddress, option, minScore,N_bowtie, L_bowtie);
//    std::ifstream Aligned(AlignedAddress.c_str());
//    if(!Aligned.is_open()){
//        cerr << "can not open your sam file"<< endl;
//    }
//    command = "";
//    command.append("LC_COLLATE=C sort -k 1 ");
//    command.append(AlignedAddress.c_str());
//    command.append(" > ");
//    command.append(ReOrderedAlignedAddress.c_str());
//    system(command.c_str());
//    std::ifstream Aligned_reordered(ReOrderedAlignedAddress.c_str());
//    readCounter=0;
//    getline(myFasts[threadCounter],fastLine);
//    Split_Token(splitted, fastLine, Header, seq, quality);
//    Split_Header(splitted, Header,mainHeader, currentFast, currentFastNumber,V_L, row, false, false, true, Pacbio, false);
//    ++fastLineCounter;
//    readCounter++;
//    while(getline(Aligned_reordered,line)){
//        if ((int)line[0]!=(int)sign){
//            Split_Token(splitted,line,Header,seq,quality);
//            Split_Header(splitted, Header, mainHeader, current, currentNumber, V_L, row, false, true, true, Pacbio, false);
//            // the current read has splitted
//            // Detecting if the current line is new in sight of read Counting...
//            isNew=true;
//            if(currentFastNumber==currentNumber){
//                isNew=false;
//                myOSams[threadCounter] << line << endl;
//            }
//            // this condition is occured for all threads except the latest one
//            while(isNew){
//                if(readCounter==length && threadCounter!=(numberOfThreads-1)){
//                    threadCounter++;
//                    readCounter=0;
//                    fastLineCounter=0;
//                    getline(myFasts[threadCounter],fastLine);
//                    Split_Token(splitted, fastLine, Header, seq, quality);
//                    Split_Header(splitted, Header, mainHeader, currentFast, currentFastNumber,V_L, row, false, false, true, Pacbio,false);
//                    ++fastLineCounter;
//                    readCounter++;
//                }
//                else{
//                    getline(myFasts[threadCounter],fastLine);
//                    Split_Token(splitted, fastLine, Header, seq, quality);
//                    Split_Header(splitted, Header,mainHeader, currentFast, currentFastNumber,V_L,row, false, false, true, Pacbio, false);
//                    ++fastLineCounter;
//                    readCounter++;
//                }
//                if(currentNumber==currentFastNumber){
//                    isNew=false;
//                    myOSams[threadCounter] << line << endl;
//                }
//            }
//        }
//    }
//    for (int thread = 0; thread < numberOfThreads; thread++) {
//        myOSams[thread].close();
//        myOFasts[thread].close();
//    }
//    //cout << "All the reads have to be processed " << filledAll << endl;
//    if(length==0){
//        Threads = new std::thread[numberOfThreads];
//        //Threads[numberOfThreads-1]=std::thread(&Assignment::Analysis, this,numberOfThreads,filledAll);
//        Analysis(0,(int)filledAll);
//        //cout << "the only thread has been created" << endl;
//    }
//    else{
//        Threads = new std::thread[numberOfThreads];
//        for (int threads = 0; threads < numberOfThreads; threads++) {
//            if (threads!=numberOfThreads-1){
//                Threads[threads]=std::thread(&Assignment::Analysis, this,threads,length);
//                //Analysis(threads,length);
//                //cout << threads <<"\t"<<length<< endl;
//            }
//            else{
//                Threads[threads]=std::thread(&Assignment::Analysis, this,threads,filledAll-(numberOfThreads-1)*length);
//                //Analysis(threads, filledAll-(numberOfThreads-1)*length);
//                //cout << threads<<"\t"<< filledAll-(numberOfThreads-1)*length<< endl;
//            }
//        }
//    }
//    //cout << "threads are created" << endl;
//    for (int T=0;T<numberOfThreads;T++){
//        Threads[T].join();
//        //cout << "THREAD "<<T <<" HAS ENDED"<< endl;
//    }

//    bool trueFlag=false;
//    int histo_filtering_error=0;
//    int histo_processed=0;
//    //int filter_Processed=0;
//    //int filtering_error=0;
//    allReads = threadVectors[0];
//    concatenate.append("cat ");
//    for(int i=1;i<numberOfThreads;i++){
//        allReads.insert(allReads.end(), threadVectors[i].begin(), threadVectors[i].end());
//    }
//    for(int i=0;i<numberOfThreads;i++){
//        ostringstream name;
//        name << outputDir + "strings" << i+1 << ".txt ";
//        O_Strings[i].close();
//        concatenate.append(name.str().c_str());
//        {
//            vector<Read>().swap(threadVectors[i]);
//        }
//    }
//    concatenate.append(("> " + outputDir + "finalStrings.txt"));
//    system(concatenate.c_str());
//    std::ifstream finalStrings((outputDir + "finalStrings.txt").c_str());
//    if(!finalStrings.is_open())
//        cerr << "is not open" << endl;
//    finalReads.clear();
//    finalReads.seekg(0, ios::beg);
//    sampleRead.readName=0;
//    sampleRead.index=0;
//    sampleRead.V_Length=0;
//    sampleRead.Passed=false;
//    while (getline(finalReads,token)){
//        Split_Token(splitted, token, Header, seq, quality);
//        Split_Header(splitted, Header, mainHeader, index, number, V_L, row, true, false, true, Pacbio, false);
//        sampleRead.index = number;
//        sampleRead.readName = index;
//        sampleRead.V_Length = V_L;
//        tinyReads.push_back(sampleRead);
//    }
//    for(int i=0;i<allReads.size();i++){
//        if(!allReads[i].Histo){
//            histo_processed++;
//            for(int k=0;k<tinyReads.size();k++){
//                if((tinyReads[k].index==allReads[i].index)&&(tinyReads[k].readName==allReads[i].readName)){
//                    allReads[i].Checked=true;
//                    tinyReads[k].Passed=true;
//                }
//            }
//            for(int j=0;j<allReads[i].scores.size();j++){
//                if(allReads[i].flags[j]==3){
//                    trueFlag = true;
//                }
//            }
//            if(!trueFlag)
//                histo_filtering_error++;
//        }
//    }

//    littleStruct tempStruct;
//    vector<littleStruct> tempArray;
//    vector<littleStruct> tempArray2;
//    vector<littleStruct> chosenArray;

//    //Saving all of the properties of the read
//    int additionalScore=0;
//    int k=0;
//    string Line;
//    while(getline(finalStrings, Line)){
//        if(allReads[k].Checked){
//            splitted = split(Line, "\t");
//            seq = splitted[splitted.size()-2];
//            quality = splitted[splitted.size()-1];
//            for(int j=0;j<allReads[k].flags.size();j++){
//                if(allReads[k].flags[j]!=-1){
//                    //additionalScore=0;
//                    additionalScore = (int)(MATCH*(seq.length())+(MISPEN+MATCH)*(allReads[k].print_M[j])+(GAPPEN+MATCH)*(allReads[k].print_I[j]));
//                    output << allReads[k].mainHeader << "\t";
//                    output << allReads[k].flags[j] << "\t";
//                    output << allReads[k].references[j] << "\t";
//                    output << allReads[k].loc[j] << "\t";
//                    output << "255 \t";
//                    output << splitted[2*j+1] << "\t";
//                    output << "* \t 0 \t 0 \t";
//                    output << seq << "\t";
//                    output << quality << "\t";
//                    output << "XA:i:" << additionalScore << "\t";
//                    output << splitted[2*j] << "\t";
//                    output << "NM:i:" <<allReads[k].print_M[j] << "\t";
//                    output << endl;
//                }
//            }
//        }
//        k++;
//    }
//    {
//        vector<Read>().swap(allReads);
//    }
//    for(int i=0;i<numberOfThreads;i++){
//        {
//            std::string remove="";
//            std::string remove2="";
//            std::string remove3="";
//            remove.append("rm ");
//            remove2.append("rm ");
//            remove3.append("rm ");
//            {
//                ostringstream filename;
//                ostringstream filename2;
//                ostringstream filename3;
//                filename << outputDir + "sam" << NumToStr(i+1) << ".txt";
//                filename2 << outputDir + "fast" << NumToStr(i+1) << ".txt";
//                filename3 << outputDir + "strings" << NumToStr(i+1) << ".txt";
//                remove.append(filename.str().c_str());
//                remove2.append(filename2.str().c_str());
//                remove3.append(filename3.str().c_str());
//            }
//            system(remove.c_str());
//            system(remove2.c_str());
//            system(remove3.c_str());
//        }
//    }
//    long failCounter=0;
//    finalReads.clear();
//    finalReads.seekg(0, ios::beg);
//    long tinyCounter=0;
//    while (getline(finalReads,token)){
//        Split_Token(splitted, token, Header, seq, quality);
//        Split_Header(splitted, Header,mainHeader, index, number, V_L, row, true, false, true, Pacbio, false);
//        if(!tinyReads[tinyCounter].Passed){
//            if(step==3 && Pacbio){
//                Passed_reads << mainHeader << endl;
//                Passed_reads << seq << endl;
//                Passed_reads << "+" << endl;
//                Passed_reads << quality << endl;
//            }else if(step==3 && !Pacbio){
//                Passed_reads << "@r" << index << "_" << number << "_" << V_L << endl;
//                Passed_reads << seq << endl;
//                Passed_reads << "+" << endl;
//                Passed_reads << quality << endl;
//            }else if(step==2){
//                Passed_reads << mainHeader << "4 \t * \t 0 \t 0 \t * \t * \t 0 \t 0 \t" << seq << "\t" << quality << "XA:i:0" << endl;
//            }
//            failCounter++;
//        }
//        tinyCounter++;
//    }
//    {
//        vector<Tiny_Read>().swap(tinyReads);
//    }
//    if(step==2){
//        system(("cat "+ outputDir +"Assign_Rem_Reads.fq >> " + outputName).c_str());
//    }
//    //cout << "End of Concatenation" << endl;
//    long finalError = 0;
//    for(int i=0;i<numberOfThreads;i++){
//        finalError = finalError + errorCollector[i];
//    }
//    system_remove(sortedFinalReadsAddress);
//    system_remove(ReadsName);
//    system_remove(failReads);
//    system_remove(outputName);
//    system_remove(AlignedAddress);
//    system_remove(ReOrderedAlignedAddress);
//    system_remove(fastQAddress);
//    system_remove(outputDir + "finalStrings.txt");
//    system_remove(outputDir + "headerFile.txt");
//    //cout << "All failures for passing to the Rem_reads file and Bowtie2: " << failCounter << endl;
//    //cout << "Error Among All The Processed Reads before histo-filtering: " << finalError << endl;
//    //cout << "All reads processed in filtering " << filter_Processed << endl;
//    //cout << "Error among the processed reads in filtering" << filtering_error << endl;
//}

//Assignment::~Assignment() {

//}

//void Assignment::Analysis(int thread, int length) {
//    /*void* returnValue;
//     int thread,length;
//     int numberOfTables;
//     struct thread_data *args;
//     args = (struct thread_data *)arguments;
//     thread = args->thread_num;
//     length = args->length;
//     numberOfTables = args->numberOfTables;*/
//    myFasts[thread].clear();
//    myFasts[thread].seekg(0, ios::beg);
//    string line, fastLine, currentSeq, currentHeader, currentQuality, samHeader, currentMainHeader, samMainHeader;
//    string tempStr, localStr;
//    bool isNew;
//    bool trueFlag=false;
//    bool hasRemained=false;
//    long long fastLineCounter=0, currentFast=0,readCounter=0, current=0;
//    long long currentFastNumber, currentNumber,lastNumber=0;
//    int row=0,currentLength,samLength, currentNumOfTables;
//    int fake_Row,right=0;
//    int seqSize;
//    int localCalling=0;
//    vector <long long> faileds;
//    vector <Path> localPaths;
//    vector <Path> chosenArray, orderedchosenarray,filteredChosenArray;//=================================
//    vector <Path> tempArray;
//    vector <string> final_empty_Cigars;
//    vector <string> final_Cigars;
//    vector <long long> truePath;
//    vector <string> splitted;
//    vector <string> samSplitted;
//    vector <string> samVector;
//    vector < vector<long long> > myTables;
//    vector < vector<string> > mySeq;
//    vector < vector<string> > myRef;
//    vector < vector<int> > myFlags;
//    vector<string> Refs;
//    vector<double> scores;
//    vector<double> Path_M;
//    vector<double> Path_I;
//    vector<double> Print_M;
//    vector<double> Print_I;
//    vector<int> scoreFlags;
//    vector<int> Flags;
//    vector<int> seqSizes;
//    vector<long long> positions;
//    vector<long> emptyIntervals;
//    vector<int> mergeParts;
//    vector<bool> quantArray;
//    double totalScore=0;
//    int total_M=0, total_I=0;
//    int remain_Score=0, remain_I=0, remain_M=0;
//    int M_cigar=0,I_cigar=0;
//    double NScore;
//    double N_Path_M;
//    double N_Path_I;

//    myTables.resize(2*checkNumOfTables(partLength, L_MAX));
//    mySeq.resize(checkNumOfTables(partLength, L_MAX));
//    myRef.resize(checkNumOfTables(partLength, L_MAX));
//    myFlags.resize(checkNumOfTables(partLength, L_MAX));
//    errorCollector[thread]=0;
//    // All we need has prepared
//    // while loop for identifying the reads and acquisition of their tables
//    // /*******************************************************************
//    getline(myFasts[thread],fastLine);
//    Split_Token(splitted, fastLine, currentHeader, currentSeq, currentQuality);
//    Split_Header(splitted,currentHeader, currentMainHeader,currentFast, currentFastNumber, currentLength,fake_Row, false, false, true, Pacbio, false);
//    fastLineCounter++;
//    readCounter++;
//    currentNumOfTables = checkNumOfTables(partLength, currentLength);
//    truePath.resize(currentNumOfTables);
//    for (int k=0;k<currentNumOfTables;k++){
//        myTables[2*k+1].push_back(currentFast+k*(partLength-d));
//    }
//    // /*******************************************************************
//    while(getline(mySams[thread],line)){
//        {
//            vector<string>().swap(samSplitted);
//            vector<string>().swap(samVector);
//        }
//        samVector = split(line, "\t");
//        samHeader = samVector[0];
//        Split_Header(samSplitted, samHeader, samMainHeader, current, currentNumber,samLength,row,false,true, false, Pacbio, true);
//        current = current + row*(partLength-d);
//        isNew=true;
//        if(currentNumber==currentFastNumber){
//            isNew=false;
//            lastNumber = currentNumber;
//        }
//        if(!isNew){
//            // Saving the read properties
//            Save_Properties(samVector,myFlags, mySeq, myTables, myRef, row, BOWTIE1);
//        }
//        // /*****************************************************************************
//        if(isNew){
//            // /*****************************************************************************
//            // start analysis of the tables
//            for(int l=0;l<currentNumOfTables;l++){
//                if(myTables[2*l+1].size()==(myDepth+1)){
//                    //vector<long long>().swap(myTables[2*l+1]);
//                    vector<long long>(myTables[2*l+1].begin(), myTables[2*l+1].begin()+1).swap(myTables[2*l+1]);//========================
//                    vector<long long>().swap(myTables[2*l]);
//                    vector<string>().swap(mySeq[l]);
//                    vector<int>().swap(myFlags[l]);
//                    vector<string>().swap(myRef[l]);
//                }
//            }
//            for(int j=0;j<currentNumOfTables;j++){
//                if(myTables[2*j+1].size()!=1){
//                    int size=(int)myTables[2*j+1].size()-1;
//                    for(int k=0;k<myTables[2*j].size();k++){
//                        myTables[2*j][k]=(exp(-0.5*(size-1)))*myTables[2*j][k];
//                    }
//                }
//            }
//            {
//                vector<Path>().swap(tempArray);
//                vector<Path>().swap(chosenArray);
//                vector<Path>().swap(filteredChosenArray);
//                vector<Path>().swap(orderedchosenarray);
//            }
//            gen_Paths(currentNumOfTables,currentLength,myTables, myFlags, myRef, mySeq, truePath, localPaths,faileds);
//            if(!localPaths.empty()){
//                Filters(localPaths,tempArray, chosenArray);

//                // Finding paths with chosen scores is done
//                //Writing on output file
//                // Completing the sequence of current Read by LocalAligner class//
//                {
//                    vector<string>().swap(final_Cigars);
//                    vector<string>().swap(final_empty_Cigars);
//                }
//                trueFlag=false;
//                right=0;
//                for (int completing=0;completing<chosenArray.size();completing++){
//                    emptyIntervals.clear();
//                    quantArray.resize(currentLength);
//                    if(chosenArray[completing].clearance()){
//                        //tempScore = chosenArray[completing].getScore();
//                        if(chosenArray[completing].hasMore(1)){
//                            filteredChosenArray.push_back(chosenArray[completing]);
//                            totalScore = 0;
//                            localCalling=0;
//                            seqSize=0;
//                            total_I=0;
//                            total_M=0;
//                            M_cigar=0;
//                            I_cigar=0;
//                            remain_I=0;
//                            remain_M=0;
//                            remain_Score=0;
//                            localStr.clear();
//                            // this path must be reported
//                            int tempFlag = chosenArray[completing].getFlag();
//                            for(int j=0;j<currentLength;j++){
//                                quantArray[j]=false;
//                            }
//                            merge_Method(mergeParts,currentNumOfTables,completing, chosenArray, quantArray, final_Cigars,M_cigar,I_cigar);
//                            emptyInterval_Detector(quantArray, &emptyIntervals, currentFast, currentLength);
//                            if(emptyIntervals.size()%2==0){
//                                hasRemained=false;
//                            }
//                            else{
//                                hasRemained=true;
//                            }
//                            //**********************************************************
//                            if(!emptyIntervals.empty()){
//                                fill_The_Gaps(chosenArray, completing, emptyIntervals,currentSeq, IndelShift, tempFlag, totalScore, seqSize, localCalling,hasRemained, final_empty_Cigars,total_M, total_I);
//                                chosenArray[completing].setTotalCigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0]));
//                            }
//                            else{
//                                chosenArray[completing].setTotalCigar(final_Cigars[0]);
//                            }

//                            int additionflag = 0;
//                            for (int chosenarrayelement = 0;chosenarrayelement<orderedchosenarray.size();chosenarrayelement++){
//                                if (orderedchosenarray[chosenarrayelement].getBirthRow() > chosenArray[completing].getBirthRow()){
//                                    orderedchosenarray.insert(orderedchosenarray.begin()+chosenarrayelement, chosenArray[completing]);
//                                    additionflag = 1;
//                                    break;
//                                }

//                            }
//                            if (additionflag == 0){
//                                orderedchosenarray.push_back(chosenArray[completing]);
//                            }
//                            //**********************************************************

//                            if(!emptyIntervals.empty() && local_Flag){
//                                fill_The_Gaps(chosenArray, completing, emptyIntervals,currentSeq, IndelShift, tempFlag, totalScore, seqSize, localCalling,hasRemained, final_empty_Cigars,total_M, total_I);
//                                NScore = totalScore/seqSize;
//                                N_Path_I = total_I;
//                                N_Path_M = total_M;
//                                if(seqSize > (int)(0.1*currentLength)){
//                                    right++;
//                                    if(chosenArray[completing].isTruePath()){
//                                        trueFlag =true;
//                                    }
//                                    tempStr = L_MAX_Gate(L_MAX,(int)currentSeq.length(),currentSeq);
//                                    if(tempStr.size()!=currentSeq.length()){
//                                        remain_Local(currentLength, chosenArray[completing],remain_M, remain_I, remain_Score,localStr, currentSeq);
//                                    }
//                                    scores.push_back(NScore);
//                                    Path_I.push_back(N_Path_I);
//                                    Path_M.push_back(N_Path_M);
//                                    Print_M.push_back(M_cigar+total_M+remain_M);
//                                    Print_I.push_back(I_cigar+total_I+remain_I);
//                                    seqSizes.push_back(seqSize);
//                                    scoreFlags.push_back(chosenArray[completing].isTruePath());
//                                    positions.push_back(chosenArray[completing].getMain());
//                                    Flags.push_back(chosenArray[completing].getFlag());
//                                    Refs.push_back(chosenArray[completing].getRef());
//                                    O_Strings[thread] << finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr << "\t";
//                                    O_Strings[thread] << mD2Cigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr) << "\t";
//                                }
//                            }
//                            else{
//                                right++;
//                                if(chosenArray[completing].isTruePath()){
//                                    trueFlag =true;
//                                }
//                                tempStr = L_MAX_Gate(L_MAX,(int)currentSeq.length(),currentSeq);
//                                if(tempStr.size()!=currentSeq.length() && local_Flag){
//                                    remain_Local(currentLength, chosenArray[completing],remain_M, remain_I, remain_Score,localStr, currentSeq);
//                                }
//                                scores.push_back((double)chosenArray[completing].getScore());
//                                Print_M.push_back(M_cigar+remain_M);
//                                Print_I.push_back(I_cigar+remain_I);
//                                scoreFlags.push_back((int)(chosenArray[completing].isTruePath()+2));
//                                positions.push_back(chosenArray[completing].getMain());
//                                Flags.push_back(chosenArray[completing].getFlag());
//                                Refs.push_back(chosenArray[completing].getRef());
//                                /*if(true){
//                                 if(final_empty_Cigars.size()>18)
//                                 out << "@#$%" << endl;

//                                 for(int j=0;j<quantArray.size();j++){
//                                 out << quantArray[j];
//                                 }
//                                 out << endl;
//                                 out << "remain " << hasRemained  << "\t" << emptyIntervals.size() << endl;

//                                 out << final_Cigars.size() << "\t" << final_empty_Cigars.size() << "\t" << currentNumOfTables << endl;

//                                 for(int i=0;i<currentNumOfTables;i++){
//                                 out << chosenArray[completing].getIndices()[i] << "\t";
//                                 }
//                                 }
//                                 out << endl;*/

//                                O_Strings[thread] << finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr << "\t";
//                                O_Strings[thread] << mD2Cigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr) << "\t";
//                            }
//                            {
//                                vector<bool>().swap(quantArray);
//                                vector<long>().swap(emptyIntervals);
//                                vector<string>().swap(final_Cigars);
//                                vector<string>().swap(final_empty_Cigars);
//                            }
//                        }
//                    }
//                }

//                for (int chosenarrayelement = 0;chosenarrayelement<orderedchosenarray.size();chosenarrayelement++){
//                    for (int chosenarrayelement1 = chosenarrayelement+1;chosenarrayelement1<orderedchosenarray.size();chosenarrayelement1++){
//                        if (orderedchosenarray[chosenarrayelement].getLastRow() < orderedchosenarray[chosenarrayelement1].getBirthRow()){
//                            if (forward_Detection(orderedchosenarray[chosenarrayelement].getFlag()) == forward_Detection(orderedchosenarray[chosenarrayelement1].getFlag())){
//                                if (forward_Detection(orderedchosenarray[chosenarrayelement].getFlag())){
//                                    if((orderedchosenarray[chosenarrayelement1].getBirth()-orderedchosenarray[chosenarrayelement].getLastIndex()<10000) && (orderedchosenarray[chosenarrayelement1].getBirth()-orderedchosenarray[chosenarrayelement].getLastIndex()>0)){
//                                        Join2pathes(orderedchosenarray[chosenarrayelement],orderedchosenarray[chosenarrayelement1]);
//                                        orderedchosenarray.erase(orderedchosenarray.begin()+chosenarrayelement1);
//                                        chosenarrayelement--;
//                                        break;
//                                    }

//                                }
//                                if (reverse_Detection(orderedchosenarray[chosenarrayelement].getFlag())){
//                                    if((orderedchosenarray[chosenarrayelement].getLastIndex()-orderedchosenarray[chosenarrayelement1].getBirth()<10000) && (orderedchosenarray[chosenarrayelement].getLastIndex()-orderedchosenarray[chosenarrayelement1].getBirth()>0)){
//                                        Join2pathes(orderedchosenarray[chosenarrayelement],orderedchosenarray[chosenarrayelement1]);
//                                        orderedchosenarray.erase(orderedchosenarray.begin()+chosenarrayelement1);
//                                        chosenarrayelement--;
//                                        break;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                /*for (int chosenarrayelement = 0; chosenarrayelement < orderedchosenarray.size(); chosenarrayelement++) {
//                    cout << lastNumber << "\t" << orderedchosenarray[chosenarrayelement].getBirth() << "\t";
//                    for (int cosenArray = 0; cosenArray < orderedchosenarray[chosenarrayelement].getIndices().size(); cosenArray++) {
//                        cout << orderedchosenarray[chosenarrayelement].getIndices().at(cosenArray) << "\t";
//                    }
//                    cout << "\n";
//                }*/

//                if (orderedchosenarray.size() == 1) {
//                    if (orderedchosenarray[0].getSpliced()) {
//                        if (this->reads[lastNumber-1].flag == -1) {
//                            if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                                this->reads[lastNumber-1].flag = 5;
//                            } else {
//                                this->reads[lastNumber-1].flag = 21;
//                            }
//                            this->reads[lastNumber-1].cigar = orderedchosenarray[0].getTotalCigar();
//                            this->reads[lastNumber-1].firstFragment = orderedchosenarray[0].getBirthRow();
//                            this->reads[lastNumber-1].lastFragment = orderedchosenarray[0].getLastRow();
//                            this->reads[lastNumber-1].firstPosition = orderedchosenarray[0].getBirth();
//                            this->reads[lastNumber-1].lastPosition  = orderedchosenarray[0].getLastIndex();

//                        } else {
//                            this->readsNext->push_back(iReadNext());
//                            if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                                this->readsNext->back().flag = 5;
//                            } else {
//                                this->readsNext->back().flag = 21;
//                            }
//                            this->readsNext->back().cigar = orderedchosenarray[0].getTotalCigar();
//                            this->readsNext->back().firstFragment = orderedchosenarray[0].getBirthRow();
//                            this->readsNext->back().lastFragment = orderedchosenarray[0].getLastRow();
//                            this->readsNext->back().firstPosition = orderedchosenarray[0].getBirth();
//                            this->readsNext->back().lastPosition  = orderedchosenarray[0].getLastIndex();
//                            this->readsNext->back().previousIndex = lastNumber;
//                        }


//                        if (reverse_Detection(orderedchosenarray[0].getFlag())) {
//                            bool firstNonZero = false;
//                            long long startOfExon = 0;
//                            int startOfExonFrag = 0;
//                            long long endOfExon = 0;
//                            int endOfExonFrag = 0;
//                            for (int partition = 0; partition < orderedchosenarray[0].getIndices().size(); partition++) {
//                                if (orderedchosenarray[0].getIndices()[partition]>0 && firstNonZero == false) {
//                                    endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                                    endOfExonFrag = partition;
//                                    firstNonZero = true;
//                                }
//                                else if(orderedchosenarray[0].getIndices()[partition]>0 && abs(endOfExon-partLength- (partition-endOfExonFrag)*(partLength-d) - orderedchosenarray[0].getIndices()[partition]) <= ((partition-startOfExonFrag)*error)) {
//                                    startOfExon = orderedchosenarray[0].getIndices()[partition];
//                                    startOfExonFrag = partition;
//                                    if (partition == orderedchosenarray[0].getIndices().size()-1) {
//                                        this->addIntervalToDepth(startOfExon, endOfExon);
//                                    }
//                                }
//                                else if (orderedchosenarray[0].getIndices()[partition]>0){
//                                    this->addIntervalToDepth(startOfExon, endOfExon);
//                                    endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                                    endOfExonFrag = partition;
//                                }

//                            }
//                        } else {
//                            bool firstNonZero = false;
//                            long long startOfExon = 0;
//                            int startOfExonFrag = 0;
//                            long long endOfExon = 0;
//                            for (int partition = 0; partition < orderedchosenarray[0].getIndices().size(); partition++) {
//                                if (orderedchosenarray[0].getIndices()[partition] > 0 && firstNonZero == false) {
//                                    startOfExon = orderedchosenarray[0].getIndices()[partition];
//                                    startOfExonFrag = partition;
//                                    firstNonZero = true;
//                                }
//                                else if(orderedchosenarray[0].getIndices()[partition] > 0 && abs(orderedchosenarray[0].getIndices()[partition] - (partition-startOfExonFrag)*(partLength-d) -startOfExon) <= ((partition-startOfExonFrag)*error)) {
//                                    endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                                    //endOfExonFrag = partition;
//                                    if (partition == orderedchosenarray[0].getIndices().size()-1) {
//                                        this->addIntervalToDepth(startOfExon, endOfExon);
//                                    }
//                                }
//                                else if (orderedchosenarray[0].getIndices()[partition] > 0){
//                                    this->addIntervalToDepth(startOfExon, endOfExon);
//                                    startOfExon = orderedchosenarray[0].getIndices()[partition];
//                                    startOfExonFrag = partition;
//                                }

//                            }
//                        }
//                        //add to spliced database
//                        //update depth
//                    }
//                    else{
//                        if (this->reads[lastNumber-1].flag == -1) {
//                            if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                                this->reads[lastNumber-1].flag = 4;
//                            } else {
//                                this->reads[lastNumber-1].flag = 20;
//                            }
//                            this->reads[lastNumber-1].firstFragment = orderedchosenarray[0].getBirthRow();
//                            this->reads[lastNumber-1].lastFragment = orderedchosenarray[0].getLastRow();
//                            this->reads[lastNumber-1].firstPosition = orderedchosenarray[0].getBirth();
//                            this->reads[lastNumber-1].lastPosition  = orderedchosenarray[0].getLastIndex();
//                            updateDepth(this->reads[lastNumber-1].firstPosition, this->reads[lastNumber-1].lastPosition, forward_Detection(orderedchosenarray[0].getFlag()), partLength);
//                        } else {
//                            this->readsNext->push_back(iReadNext());
//                            if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                                this->readsNext->back().flag = 4;
//                            } else {
//                                this->readsNext->back().flag = 20;
//                            }
//                            this->readsNext->back().firstFragment = orderedchosenarray[0].getBirthRow();
//                            this->readsNext->back().lastFragment = orderedchosenarray[0].getLastRow();
//                            this->readsNext->back().firstPosition = orderedchosenarray[0].getBirth();
//                            this->readsNext->back().lastPosition  = orderedchosenarray[0].getLastIndex();

//                            this->readsNext->back().previousIndex = lastNumber;
//                            updateDepth(this->readsNext->back().firstPosition, this->readsNext->back().lastPosition, forward_Detection(orderedchosenarray[0].getFlag()), partLength);
//                        }
//                        //add to unspliced database
//                        //update depth
//                    }
//                }
//                else{
//                    for (int element = 0; element < orderedchosenarray.size(); element++) {
//                        gene tempGene = gene();
//                        tempGene.inRead.push_back(splicedRead());
//                        tempGene.inRead.back().index = lastNumber;
//                        tempGene.inRead.back().cigar = orderedchosenarray[element].getTotalCigar();
//                        if (reverse_Detection(orderedchosenarray[element].getFlag())) {
//                            tempGene.start = orderedchosenarray.at(element).getLastIndex();
//                            tempGene.end = orderedchosenarray.at(element).getBirth();
//                            if (orderedchosenarray[element].getSpliced()) {
//                                bool firstNonZero = false;
//                                long long startOfExon = 0;
//                                int startOfExonFrag = 0;
//                                long long endOfExon = 0;
//                                int endOfExonFrag = 0;
//                                for (int partition = 0; partition < orderedchosenarray[element].getIndices().size(); partition++) {
//                                    if (orderedchosenarray[element].getIndices()[partition]>0 && firstNonZero == false) {
//                                        endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                        endOfExonFrag = partition;
//                                        firstNonZero = true;
//                                    }
//                                    else if(orderedchosenarray[element].getIndices()[partition]>0 && abs(endOfExon-partLength- (partition-endOfExonFrag)*(partLength-d) - orderedchosenarray[element].getIndices()[partition]) <= ((partition-startOfExonFrag)*error)) {
//                                        startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                        startOfExonFrag = partition;
//                                        if (partition == orderedchosenarray[element].getIndices().size()-1) {
//                                            tempGene.inExon.push_back(interval());
//                                            tempGene.inExon.back().start = startOfExon;
//                                            tempGene.inExon.back().end = endOfExon;
//                                            tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                            tempGene.inRead.back().fragments.back().firstFragment = endOfExonFrag;
//                                            tempGene.inRead.back().fragments.back().lastFragment = startOfExonFrag;
//                                            tempGene.inRead.back().fragments.back().firstPosition = endOfExon - partLength;
//                                            tempGene.inRead.back().fragments.back().lastPosition = startOfExon;
//                                        }
//                                    }
//                                    else if (orderedchosenarray[element].getIndices()[partition]>0){
//                                        tempGene.inExon.push_back(interval());
//                                        tempGene.inExon.back().start = startOfExon;
//                                        tempGene.inExon.back().end = endOfExon;
//                                        tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                        tempGene.inRead.back().fragments.back().firstFragment = endOfExonFrag;
//                                        tempGene.inRead.back().fragments.back().lastFragment = startOfExonFrag;
//                                        tempGene.inRead.back().fragments.back().firstPosition = endOfExon - partLength;
//                                        tempGene.inRead.back().fragments.back().lastPosition = startOfExon;
//                                        endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                        endOfExonFrag = partition;
//                                    }
//                                }
//                            }
//                            else{
//                                tempGene.inExon.push_back(interval());
//                                tempGene.inExon.back().start = orderedchosenarray.at(element).getLastIndex();
//                                tempGene.inExon.back().end = orderedchosenarray.at(element).getBirth() + partLength;
//                                tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                tempGene.inRead.back().fragments.back().firstFragment = orderedchosenarray.at(element).getBirthRow();
//                                tempGene.inRead.back().fragments.back().lastFragment = orderedchosenarray.at(element).getLastRow();
//                                tempGene.inRead.back().fragments.back().firstPosition = orderedchosenarray.at(element).getBirth();
//                                tempGene.inRead.back().fragments.back().lastPosition = orderedchosenarray.at(element).getLastIndex();
//                            }
//                        } else {
//                            tempGene.start = orderedchosenarray.at(element).getBirth();
//                            tempGene.end = orderedchosenarray.at(element).getLastIndex();
//                            if (orderedchosenarray[element].getSpliced()) {
//                                bool firstNonZero = false;
//                                long long startOfExon = 0;
//                                int startOfExonFrag = 0;
//                                long long endOfExon = 0;
//                                int endOfExonFrag = 0;
//                                for (int partition = 0; partition < orderedchosenarray[element].getIndices().size(); partition++) {
//                                    if (orderedchosenarray[element].getIndices()[partition] > 0 && firstNonZero == false) {
//                                        startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                        startOfExonFrag = partition;
//                                        firstNonZero = true;
//                                    }
//                                    else if(orderedchosenarray[element].getIndices()[partition] > 0 && abs(orderedchosenarray[element].getIndices()[partition] - (partition-startOfExonFrag)*(partLength-d) -startOfExon) <= ((partition-startOfExonFrag)*error)) {
//                                        endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                        endOfExonFrag = partition;
//                                        if (partition == orderedchosenarray[element].getIndices().size()-1) {
//                                            tempGene.inExon.push_back(interval());
//                                            tempGene.inExon.back().start = startOfExon;
//                                            tempGene.inExon.back().end = endOfExon;
//                                            tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                            tempGene.inRead.back().fragments.back().firstFragment = startOfExonFrag;
//                                            tempGene.inRead.back().fragments.back().lastFragment = endOfExonFrag;
//                                            tempGene.inRead.back().fragments.back().firstPosition = startOfExon;
//                                            tempGene.inRead.back().fragments.back().lastPosition = endOfExon - partLength;
//                                        }
//                                    }
//                                    else if (orderedchosenarray[element].getIndices()[partition] > 0){
//                                        tempGene.inExon.push_back(interval());
//                                        tempGene.inExon.back().start = startOfExon;
//                                        tempGene.inExon.back().end = endOfExon;
//                                        tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                        tempGene.inRead.back().fragments.back().firstFragment = startOfExonFrag;
//                                        tempGene.inRead.back().fragments.back().lastFragment = endOfExonFrag;
//                                        tempGene.inRead.back().fragments.back().firstPosition = startOfExon;
//                                        tempGene.inRead.back().fragments.back().lastPosition = endOfExon - partLength;
//                                        startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                        startOfExonFrag = partition;
//                                    }
//                                }
//                            }
//                            else{
//                                tempGene.inExon.push_back(interval());
//                                tempGene.inExon.back().start = orderedchosenarray.at(element).getBirth();
//                                tempGene.inExon.back().end = orderedchosenarray.at(element).getLastIndex() + partLength;
//                                tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                tempGene.inRead.back().fragments.back().firstFragment = orderedchosenarray.at(element).getBirthRow();
//                                tempGene.inRead.back().fragments.back().lastFragment = orderedchosenarray.at(element).getLastRow();
//                                tempGene.inRead.back().fragments.back().firstPosition = orderedchosenarray.at(element).getBirth();
//                                tempGene.inRead.back().fragments.back().lastPosition = orderedchosenarray.at(element).getLastIndex();
//                            }
//                        }
//                        tempGene.isForward = forward_Detection(orderedchosenarray[element].getFlag());
//                        bool isCompatibale = false;
//                        for (int counterOnInExon = 0; counterOnInExon < tempGene.inExon.size() && !isCompatibale; counterOnInExon++) {
//                            long long startOnInExon = tempGene.inExon.at(counterOnInExon).start;
//                            long long endOnInExon = tempGene.inExon.at(counterOnInExon).end;
//                            for (long long counterOnRealExon = 0; counterOnRealExon < this->iDepth->size(); counterOnRealExon++) {
//                                if ((this->iDepth->at(counterOnRealExon).start < startOnInExon &&
//                                     startOnInExon < this->iDepth->at(counterOnRealExon).end) ||
//                                    (this->iDepth->at(counterOnRealExon).start < endOnInExon &&
//                                     endOnInExon < this->iDepth->at(counterOnRealExon).end)) {
//                                        isCompatibale = true;
//                                        break;
//                                    }
//                            }
//                        }

//                        if (isCompatibale && !(tempGene.start == tempGene.end && tempGene.start == 0)) {
//                            //add to spliced database and update depth
//                            for (int counterOnGeneTempForCigar = 0; counterOnGeneTempForCigar < tempGene.inRead.size(); counterOnGeneTempForCigar++) {

//                                long long index = tempGene.inRead.at(counterOnGeneTempForCigar).index-1;
//                                /*string cigarForSplice = "";
//                                 int lengthOfMatch = 0;
//                                 int fragmentAlignedBefore = 0;
//                                 long long firstPositionBefore = 0;

//                                 for (int counterFragsOnInRead = 0; counterFragsOnInRead < tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size(); counterFragsOnInRead++) {

//                                 long long firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).firstPosition;
//                                 int firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).firstFragment;
//                                 int lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).lastFragment;

//                                 if (counterFragsOnInRead == 0) {
//                                 lengthOfMatch = (lastFragment-firstFragment+1)*partLength-(lastFragment-firstFragment)*d;
//                                 cigarForSplice += this->convertNumToStr(lengthOfMatch) + "M";
//                                 fragmentAlignedBefore = firstFragment;
//                                 } else {
//                                 cigarForSplice += this->convertNumToStr(firstPosition-firstPositionBefore-lengthOfMatch) + "N";
//                                 cigarForSplice += this->convertNumToStr((firstFragment-fragmentAlignedBefore)*(partLength-d)) + "R";
//                                 fragmentAlignedBefore = firstFragment;
//                                 lengthOfMatch = (lastFragment-firstFragment+1)*partLength-(lastFragment-firstFragment)*d;
//                                 cigarForSplice += this->convertNumToStr(lengthOfMatch) + "M";
//                                 }
//                                 firstPositionBefore = firstPosition;
//                                 }*/

//                                long size = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size();

//                                if (this->reads[index].flag == -1) {
//                                    this->reads[index].cigar = tempGene.inRead.at(counterOnGeneTempForCigar).cigar;
//                                    this->reads[index].firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstPosition;
//                                    this->reads[index].lastPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastPosition;
//                                    this->reads[index].firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstFragment;
//                                    this->reads[index].lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastFragment;
//                                    if (tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size() > 1) {
//                                        if (tempGene.isForward) {
//                                            this->reads[index].flag = 7;
//                                        } else {
//                                            this->reads[index].flag = 23;
//                                        }
//                                    }
//                                    else{
//                                        if (tempGene.isForward) {
//                                            this->reads[index].flag = 6;
//                                        } else {
//                                            this->reads[index].flag = 22;
//                                        }
//                                    }

//                                } else {
//                                    this->readsNext->push_back(iReadNext());
//                                    this->readsNext->back().cigar = tempGene.inRead.at(counterOnGeneTempForCigar).cigar;
//                                    this->readsNext->back().firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstPosition;
//                                    this->readsNext->back().lastPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastPosition;
//                                    this->readsNext->back().firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstFragment;
//                                    this->readsNext->back().lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastFragment;
//                                    if (tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size() > 1) {
//                                        if (tempGene.isForward) {
//                                            this->readsNext->back().flag = 7;
//                                        } else {
//                                            this->readsNext->back().flag = 23;
//                                        }
//                                    }
//                                    else{
//                                        if (tempGene.isForward) {
//                                            this->readsNext->back().flag = 6;
//                                        } else {
//                                            this->readsNext->back().flag = 22;
//                                        }
//                                    }
//                                    this->readsNext->back().previousIndex = lastNumber;
//                                }

//                                for (int counterOnExonsGenesForAddingToDepth = 0; counterOnExonsGenesForAddingToDepth < tempGene.inExon.size(); counterOnExonsGenesForAddingToDepth++) {
//                                    this->addIntervalToDepth(tempGene.inExon.at(counterOnExonsGenesForAddingToDepth).start, tempGene.inExon.at(counterOnExonsGenesForAddingToDepth).end);
//                                }
//                            }
//                        }
//                        else if (!(tempGene.start == tempGene.end && tempGene.start == 0)){
//                            AddToJenes(tempGene);
//                        }

//                    }

//                }


//                if(right!=0){
//                    {
//                        vector<double>().swap(threadReads[thread].scores);
//                        vector<double>().swap(threadReads[thread].T_I_vector);
//                        vector<double>().swap(threadReads[thread].T_M_vector);
//                        vector<int>().swap(threadReads[thread].print_M);
//                        vector<int>().swap(threadReads[thread].print_I);
//                        vector<int>().swap(threadReads[thread].flags);
//                        vector<int>().swap(threadReads[thread].directions);
//                        vector<int>().swap(threadReads[thread].seqSizes);                                                                                               vector<long long>().swap(threadReads[thread].loc);
//                        vector<string>().swap(threadReads[thread].references);
//                    }
//                    threadReads[thread].readName=currentFast;
//                    threadReads[thread].index = currentFastNumber;
//                    threadReads[thread].mainHeader = currentMainHeader;
//                    threadReads[thread].Histo = false;
//                    threadReads[thread].Checked = false;
//                    O_Strings[thread] << currentSeq << "\t" << currentQuality << "\t" << endl;
//                    if(!trueFlag)
//                        errorCollector[thread]++;
//                    for(int j=0;j<scores.size();j++){
//                        threadReads[thread].scores.push_back(scores[j]);
//                        threadReads[thread].flags.push_back(scoreFlags[j]);
//                        if(!Path_I.empty()){
//                            threadReads[thread].T_I_vector.push_back(Path_I[j]);
//                            threadReads[thread].T_M_vector.push_back(Path_M[j]);
//                            threadReads[thread].seqSizes.push_back(seqSizes[j]);
//                        }
//                        threadReads[thread].print_I.push_back(Print_I[j]);
//                        threadReads[thread].print_M.push_back(Print_M[j]);
//                        threadReads[thread].loc.push_back(positions[j]);
//                        threadReads[thread].references.push_back(Refs[j]);
//                        threadReads[thread].directions.push_back(Flags[j]);
//                        if(scoreFlags[j]==1 || scoreFlags[j]==0)
//                            threadReads[thread].Histo=true;
//                    }
//                    if(!threadReads[thread].scores.empty()){
//                        threadVectors[thread].push_back(threadReads[thread]);
//                    }
//                    scores.clear();
//                    Path_I.clear();
//                    Path_M.clear();
//                    Print_M.clear();
//                    Print_I.clear();
//                    scoreFlags.clear();
//                    seqSizes.clear();
//                    positions.clear();
//                    Refs.clear();
//                    Flags.clear();
//                }
//            }
//            // after analysis the new read we should clear all the vectors
//            for(int i=0;i<currentNumOfTables;i++){
//                vector<string>().swap(mySeq[i]);
//                vector<int>().swap(myFlags[i]);
//                vector<string>().swap(myRef[i]);
//                vector<long long>().swap(myTables[2*i]);
//                vector<long long>().swap(myTables[2*i+1]);
//            }
//            {
//                vector<Path>().swap(localPaths);
//            }
//        }
//        // /*****************************************************************************
//        while(isNew){
//            if(readCounter==length){
//                return;
//            }
//            else{
//                getline(myFasts[thread],fastLine);
//                Split_Token(splitted, fastLine, currentHeader, currentSeq, currentQuality);
//                Split_Header(splitted,currentHeader, currentMainHeader,currentFast, currentFastNumber, currentLength,fake_Row, false, false, true, Pacbio, false);
//                fastLineCounter++;
//                readCounter++;
//                currentNumOfTables = checkNumOfTables(partLength, currentLength);
//                truePath.resize(currentNumOfTables);
//                for (int k=0;k<currentNumOfTables;k++){
//                    myTables[2*k+1].push_back(currentFast+k*(partLength-d));
//                }
//            }
//            if(currentNumber==currentFastNumber){
//                isNew=false;
//                lastNumber = currentNumber;
//            }
//            if(!isNew){
//                // Saving the read properties for the first time
//                Save_Properties(samVector,myFlags, mySeq, myTables, myRef, row, BOWTIE1);
//            }
//        }
//    }
//    for(int l=0;l<currentNumOfTables;l++){
//        if(myTables[2*l+1].size()==(myDepth+1)){
//            vector<long long>().swap(myTables[2*l+1]);
//            vector<long long>().swap(myTables[2*l]);
//            vector<string>().swap(mySeq[l]);
//            vector<int>().swap(myFlags[l]);
//            vector<string>().swap(myRef[l]);
//        }
//    }
//    for(int j=0;j<currentNumOfTables;j++){
//        if(myTables[2*j+1].size()!=1){
//            int size=(int)myTables[2*j+1].size()-1;
//            for(int k=0;k<myTables[2*j].size();k++){
//                myTables[2*j][k]=(exp(-0.5*(size-1)))*myTables[2*j][k];
//            }
//        }
//    }
//    {
//        vector<Path>().swap(tempArray);
//        vector<Path>().swap(chosenArray);
//    }
//    gen_Paths(currentNumOfTables,currentLength,myTables, myFlags, myRef, mySeq, truePath, localPaths,faileds);
//    // /*************************************
//    // /************************************
//    if(!localPaths.empty()){
//        Filters(localPaths,tempArray, chosenArray);
//        // Finding paths with chosen scores is done
//        //Writing on output file
//        // Completing the sequence of current Read by LocalAligner class//
//        {
//            vector<string>().swap(final_Cigars);
//            vector<string>().swap(final_empty_Cigars);
//        }
//        trueFlag=false;
//        right=0;
//        for (int completing=0;completing<chosenArray.size();completing++){
//            emptyIntervals.clear();
//            quantArray.resize(currentLength);
//            if(chosenArray[completing].clearance()){
//                //tempScore = chosenArray[completing].getScore();
//                if(chosenArray[completing].hasMore(1)){
//                    filteredChosenArray.push_back(chosenArray[completing]);
//                    totalScore = 0;
//                    localCalling=0;
//                    seqSize=0;
//                    total_I=0;
//                    total_M=0;
//                    M_cigar=0;
//                    I_cigar=0;
//                    remain_I=0;
//                    remain_M=0;
//                    remain_Score=0;
//                    localStr.clear();
//                    // this path must be reported
//                    int tempFlag = chosenArray[completing].getFlag();
//                    for(int j=0;j<currentLength;j++){
//                        quantArray[j]=false;
//                    }
//                    merge_Method(mergeParts,currentNumOfTables,completing, chosenArray, quantArray, final_Cigars,M_cigar,I_cigar);
//                    emptyInterval_Detector(quantArray, &emptyIntervals, currentFast, currentLength);
//                    if(emptyIntervals.size()%2==0){
//                        hasRemained=false;
//                    }
//                    else{
//                        hasRemained=true;
//                    }
//                    //**********************************************************
//                    if(final_Cigars.size() > 1){
//                        fill_The_Gaps(chosenArray, completing, emptyIntervals,currentSeq, IndelShift, tempFlag, totalScore, seqSize, localCalling,hasRemained, final_empty_Cigars,total_M, total_I);
//                        chosenArray[completing].setTotalCigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0]));
//                    }
//                    else{
//                        chosenArray[completing].setTotalCigar(final_Cigars[0]);
//                    }

//                    int additionflag = 0;
//                    for (int chosenarrayelement = 0;chosenarrayelement<orderedchosenarray.size();chosenarrayelement++){
//                        if (orderedchosenarray[chosenarrayelement].getBirthRow() > chosenArray[completing].getBirthRow()){
//                            orderedchosenarray.insert(orderedchosenarray.begin()+chosenarrayelement, chosenArray[completing]);
//                            additionflag = 1;
//                            break;
//                        }

//                    }
//                    if (additionflag == 0){
//                        orderedchosenarray.push_back(chosenArray[completing]);
//                    }
//                    //**********************************************************

//                    if(!emptyIntervals.empty() && local_Flag){
//                        fill_The_Gaps(chosenArray, completing, emptyIntervals,currentSeq, IndelShift, tempFlag, totalScore, seqSize, localCalling,hasRemained, final_empty_Cigars,total_M, total_I);
//                        NScore = totalScore/seqSize;
//                        N_Path_I = total_I;
//                        N_Path_M = total_M;
//                        if(seqSize > (int)(0.1*currentLength)){
//                            right++;
//                            if(chosenArray[completing].isTruePath()){
//                                trueFlag =true;
//                            }
//                            tempStr = L_MAX_Gate(L_MAX,(int)currentSeq.length(),currentSeq);
//                            if(tempStr.size()!=currentSeq.length()){
//                                remain_Local(currentLength, chosenArray[completing],remain_M, remain_I, remain_Score,localStr, currentSeq);
//                            }
//                            scores.push_back(NScore);
//                            Path_I.push_back(N_Path_I);
//                            Path_M.push_back(N_Path_M);
//                            Print_M.push_back(M_cigar+total_M+remain_M);
//                            Print_I.push_back(I_cigar+total_I+remain_I);
//                            seqSizes.push_back(seqSize);
//                            scoreFlags.push_back(chosenArray[completing].isTruePath());
//                            positions.push_back(chosenArray[completing].getMain());
//                            Flags.push_back(chosenArray[completing].getFlag());
//                            Refs.push_back(chosenArray[completing].getRef());
//                            O_Strings[thread] << finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr << "\t";
//                            O_Strings[thread] << mD2Cigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr) << "\t";
//                        }
//                    }
//                    else{
//                        right++;
//                        if(chosenArray[completing].isTruePath()){
//                            trueFlag =true;
//                        }
//                        tempStr = L_MAX_Gate(L_MAX,(int)currentSeq.length(),currentSeq);
//                        if(tempStr.size()!=currentSeq.length() && local_Flag){
//                            remain_Local(currentLength, chosenArray[completing],remain_M, remain_I, remain_Score,localStr, currentSeq);
//                        }
//                        scores.push_back((double)chosenArray[completing].getScore());
//                        Print_M.push_back(M_cigar+remain_M);
//                        Print_I.push_back(I_cigar+remain_I);
//                        scoreFlags.push_back((int)(chosenArray[completing].isTruePath()+2));
//                        positions.push_back(chosenArray[completing].getMain());
//                        Flags.push_back(chosenArray[completing].getFlag());
//                        Refs.push_back(chosenArray[completing].getRef());
//                        O_Strings[thread] << finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr << "\t";
//                        O_Strings[thread] << mD2Cigar(finalString(final_Cigars,final_empty_Cigars,chosenArray[completing].getIndices()[0])+localStr) << "\t";
//                    }
//                    {
//                        vector<bool>().swap(quantArray);
//                        vector<long>().swap(emptyIntervals);
//                        vector<string>().swap(final_Cigars);
//                        vector<string>().swap(final_empty_Cigars);
//                    }
//                }
//            }
//        }

//        for (int chosenarrayelement = 0;chosenarrayelement<orderedchosenarray.size();chosenarrayelement++){
//            for (int chosenarrayelement1 = chosenarrayelement+1;chosenarrayelement1<orderedchosenarray.size();chosenarrayelement1++){
//                if (orderedchosenarray[chosenarrayelement].getLastRow() < orderedchosenarray[chosenarrayelement1].getBirthRow()){
//                    if (forward_Detection(orderedchosenarray[chosenarrayelement].getFlag()) == forward_Detection(orderedchosenarray[chosenarrayelement1].getFlag())){
//                        if (forward_Detection(orderedchosenarray[chosenarrayelement].getFlag())){
//                            if((orderedchosenarray[chosenarrayelement1].getBirth()-orderedchosenarray[chosenarrayelement].getLastIndex()<10000) && (orderedchosenarray[chosenarrayelement1].getBirth()-orderedchosenarray[chosenarrayelement].getLastIndex()>0)){
//                                Join2pathes(orderedchosenarray[chosenarrayelement],orderedchosenarray[chosenarrayelement1]);
//                                orderedchosenarray.erase(orderedchosenarray.begin()+chosenarrayelement1);
//                                chosenarrayelement--;
//                                break;
//                            }

//                        }
//                        if (reverse_Detection(orderedchosenarray[chosenarrayelement].getFlag())){
//                            if((orderedchosenarray[chosenarrayelement].getLastIndex()-orderedchosenarray[chosenarrayelement1].getBirth()<10000) && (orderedchosenarray[chosenarrayelement].getLastIndex()-orderedchosenarray[chosenarrayelement1].getBirth()>0)){
//                                Join2pathes(orderedchosenarray[chosenarrayelement],orderedchosenarray[chosenarrayelement1]);
//                                orderedchosenarray.erase(orderedchosenarray.begin()+chosenarrayelement1);
//                                chosenarrayelement--;
//                                break;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//        /*for (int chosenarrayelement = 0; chosenarrayelement < orderedchosenarray.size(); chosenarrayelement++) {
//            cout << lastNumber << "\t" << orderedchosenarray[chosenarrayelement].getBirth() << "\t";
//            for (int cosenArray = 0; cosenArray < orderedchosenarray[chosenarrayelement].getIndices().size(); cosenArray++) {
//                cout << orderedchosenarray[chosenarrayelement].getIndices().at(cosenArray) << "\t";
//            }
//            cout << "\n";
//        }*/

//        if (orderedchosenarray.size() == 1) {
//            if (orderedchosenarray[0].getSpliced()) {
//                if (this->reads[lastNumber-1].flag == -1) {
//                    if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                        this->reads[lastNumber-1].flag = 5;
//                    } else {
//                        this->reads[lastNumber-1].flag = 21;
//                    }
//                    this->reads[lastNumber-1].cigar = orderedchosenarray[0].getTotalCigar();
//                    this->reads[lastNumber-1].firstFragment = orderedchosenarray[0].getBirthRow();
//                    this->reads[lastNumber-1].lastFragment = orderedchosenarray[0].getLastRow();
//                    this->reads[lastNumber-1].firstPosition = orderedchosenarray[0].getBirth();
//                    this->reads[lastNumber-1].lastPosition  = orderedchosenarray[0].getLastIndex();

//                } else {
//                    this->readsNext->push_back(iReadNext());
//                    if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                        this->readsNext->back().flag = 5;
//                    } else {
//                        this->readsNext->back().flag = 21;
//                    }
//                    this->readsNext->back().cigar = orderedchosenarray[0].getTotalCigar();
//                    this->readsNext->back().firstFragment = orderedchosenarray[0].getBirthRow();
//                    this->readsNext->back().lastFragment = orderedchosenarray[0].getLastRow();
//                    this->readsNext->back().firstPosition = orderedchosenarray[0].getBirth();
//                    this->readsNext->back().lastPosition  = orderedchosenarray[0].getLastIndex();
//                    this->readsNext->back().previousIndex = lastNumber;
//                }


//                if (reverse_Detection(orderedchosenarray[0].getFlag())) {
//                    bool firstNonZero = false;
//                    long long startOfExon = 0;
//                    int startOfExonFrag = 0;
//                    long long endOfExon = 0;
//                    int endOfExonFrag = 0;
//                    for (int partition = 0; partition < orderedchosenarray[0].getIndices().size(); partition++) {
//                        if (orderedchosenarray[0].getIndices()[partition]>0 && firstNonZero == false) {
//                            endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                            endOfExonFrag = partition;
//                            firstNonZero = true;
//                        }
//                        else if(orderedchosenarray[0].getIndices()[partition]>0 && abs(endOfExon-partLength- (partition-endOfExonFrag)*(partLength-d) - orderedchosenarray[0].getIndices()[partition]) <= ((partition-startOfExonFrag)*error)) {
//                            startOfExon = orderedchosenarray[0].getIndices()[partition];
//                            startOfExonFrag = partition;
//                            if (partition == orderedchosenarray[0].getIndices().size()-1) {
//                                this->addIntervalToDepth(startOfExon, endOfExon);
//                            }
//                        }
//                        else if (orderedchosenarray[0].getIndices()[partition]>0){
//                            this->addIntervalToDepth(startOfExon, endOfExon);
//                            endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                            endOfExonFrag = partition;
//                        }

//                    }
//                } else {
//                    bool firstNonZero = false;
//                    long long startOfExon = 0;
//                    int startOfExonFrag = 0;
//                    long long endOfExon = 0;
//                    for (int partition = 0; partition < orderedchosenarray[0].getIndices().size(); partition++) {
//                        if (orderedchosenarray[0].getIndices()[partition] > 0 && firstNonZero == false) {
//                            startOfExon = orderedchosenarray[0].getIndices()[partition];
//                            startOfExonFrag = partition;
//                            firstNonZero = true;
//                        }
//                        else if(orderedchosenarray[0].getIndices()[partition] > 0 && abs(orderedchosenarray[0].getIndices()[partition] - (partition-startOfExonFrag)*(partLength-d) -startOfExon) <= ((partition-startOfExonFrag)*error)) {
//                            endOfExon = orderedchosenarray[0].getIndices()[partition] + partLength;
//                            //endOfExonFrag = partition;
//                            if (partition == orderedchosenarray[0].getIndices().size()-1) {
//                                this->addIntervalToDepth(startOfExon, endOfExon);
//                            }
//                        }
//                        else if (orderedchosenarray[0].getIndices()[partition] > 0){
//                            this->addIntervalToDepth(startOfExon, endOfExon);
//                            startOfExon = orderedchosenarray[0].getIndices()[partition];
//                            startOfExonFrag = partition;
//                        }

//                    }
//                }
//                //add to spliced database
//                //update depth
//            }
//            else{
//                if (this->reads[lastNumber-1].flag == -1) {
//                    if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                        this->reads[lastNumber-1].flag = 4;
//                    } else {
//                        this->reads[lastNumber-1].flag = 20;
//                    }
//                    this->reads[lastNumber-1].firstFragment = orderedchosenarray[0].getBirthRow();
//                    this->reads[lastNumber-1].lastFragment = orderedchosenarray[0].getLastRow();
//                    this->reads[lastNumber-1].firstPosition = orderedchosenarray[0].getBirth();
//                    this->reads[lastNumber-1].lastPosition  = orderedchosenarray[0].getLastIndex();
//                    updateDepth(this->reads[lastNumber-1].firstPosition, this->reads[lastNumber-1].lastPosition, forward_Detection(orderedchosenarray[0].getFlag()), partLength);
//                } else {
//                    this->readsNext->push_back(iReadNext());
//                    if (forward_Detection(orderedchosenarray[0].getFlag()) == true) {
//                        this->readsNext->back().flag = 4;
//                    } else {
//                        this->readsNext->back().flag = 20;
//                    }
//                    this->readsNext->back().firstFragment = orderedchosenarray[0].getBirthRow();
//                    this->readsNext->back().lastFragment = orderedchosenarray[0].getLastRow();
//                    this->readsNext->back().firstPosition = orderedchosenarray[0].getBirth();
//                    this->readsNext->back().lastPosition  = orderedchosenarray[0].getLastIndex();

//                    this->readsNext->back().previousIndex = lastNumber;
//                    updateDepth(this->readsNext->back().firstPosition, this->readsNext->back().lastPosition, forward_Detection(orderedchosenarray[0].getFlag()), partLength);
//                }
//                //add to unspliced database
//                //update depth
//            }
//        }
//        else{
//            for (int element = 0; element < orderedchosenarray.size(); element++) {
//                gene tempGene = gene();
//                tempGene.inRead.push_back(splicedRead());
//                tempGene.inRead.back().index = lastNumber;
//                tempGene.inRead.back().cigar = orderedchosenarray[element].getTotalCigar();
//                if (reverse_Detection(orderedchosenarray[element].getFlag())) {
//                    tempGene.start = orderedchosenarray.at(element).getLastIndex();
//                    tempGene.end = orderedchosenarray.at(element).getBirth();
//                    if (orderedchosenarray[element].getSpliced()) {
//                        bool firstNonZero = false;
//                        long long startOfExon = 0;
//                        int startOfExonFrag = 0;
//                        long long endOfExon = 0;
//                        int endOfExonFrag = 0;
//                        for (int partition = 0; partition < orderedchosenarray[element].getIndices().size(); partition++) {
//                            if (orderedchosenarray[element].getIndices()[partition]>0 && firstNonZero == false) {
//                                endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                endOfExonFrag = partition;
//                                firstNonZero = true;
//                            }
//                            else if(orderedchosenarray[element].getIndices()[partition]>0 && abs(endOfExon-partLength- (partition-endOfExonFrag)*(partLength-d) - orderedchosenarray[element].getIndices()[partition]) <= ((partition-startOfExonFrag)*error)) {
//                                startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                startOfExonFrag = partition;
//                                if (partition == orderedchosenarray[element].getIndices().size()-1) {
//                                    tempGene.inExon.push_back(interval());
//                                    tempGene.inExon.back().start = startOfExon;
//                                    tempGene.inExon.back().end = endOfExon;
//                                    tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                    tempGene.inRead.back().fragments.back().firstFragment = endOfExonFrag;
//                                    tempGene.inRead.back().fragments.back().lastFragment = startOfExonFrag;
//                                    tempGene.inRead.back().fragments.back().firstPosition = endOfExon - partLength;
//                                    tempGene.inRead.back().fragments.back().lastPosition = startOfExon;
//                                }
//                            }
//                            else if (orderedchosenarray[element].getIndices()[partition]>0){
//                                tempGene.inExon.push_back(interval());
//                                tempGene.inExon.back().start = startOfExon;
//                                tempGene.inExon.back().end = endOfExon;
//                                tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                tempGene.inRead.back().fragments.back().firstFragment = endOfExonFrag;
//                                tempGene.inRead.back().fragments.back().lastFragment = startOfExonFrag;
//                                tempGene.inRead.back().fragments.back().firstPosition = endOfExon - partLength;
//                                tempGene.inRead.back().fragments.back().lastPosition = startOfExon;
//                                endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                endOfExonFrag = partition;
//                            }
//                        }
//                    }
//                    else{
//                        tempGene.inExon.push_back(interval());
//                        tempGene.inExon.back().start = orderedchosenarray.at(element).getLastIndex();
//                        tempGene.inExon.back().end = orderedchosenarray.at(element).getBirth() + partLength;
//                        tempGene.inRead.back().fragments.push_back(splicedFrag());
//                        tempGene.inRead.back().fragments.back().firstFragment = orderedchosenarray.at(element).getBirthRow();
//                        tempGene.inRead.back().fragments.back().lastFragment = orderedchosenarray.at(element).getLastRow();
//                        tempGene.inRead.back().fragments.back().firstPosition = orderedchosenarray.at(element).getBirth();
//                        tempGene.inRead.back().fragments.back().lastPosition = orderedchosenarray.at(element).getLastIndex();
//                    }
//                } else {
//                    tempGene.start = orderedchosenarray.at(element).getBirth();
//                    tempGene.end = orderedchosenarray.at(element).getLastIndex();
//                    if (orderedchosenarray[element].getSpliced()) {
//                        bool firstNonZero = false;
//                        long long startOfExon = 0;
//                        int startOfExonFrag = 0;
//                        long long endOfExon = 0;
//                        int endOfExonFrag = 0;
//                        for (int partition = 0; partition < orderedchosenarray[element].getIndices().size(); partition++) {
//                            if (orderedchosenarray[element].getIndices()[partition] > 0 && firstNonZero == false) {
//                                startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                startOfExonFrag = partition;
//                                firstNonZero = true;
//                            }
//                            else if(orderedchosenarray[element].getIndices()[partition] > 0 && abs(orderedchosenarray[element].getIndices()[partition] - (partition-startOfExonFrag)*(partLength-d) -startOfExon) <= ((partition-startOfExonFrag)*error)) {
//                                endOfExon = orderedchosenarray[element].getIndices()[partition] + partLength;
//                                endOfExonFrag = partition;
//                                if (partition == orderedchosenarray[element].getIndices().size()-1) {
//                                    tempGene.inExon.push_back(interval());
//                                    tempGene.inExon.back().start = startOfExon;
//                                    tempGene.inExon.back().end = endOfExon;
//                                    tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                    tempGene.inRead.back().fragments.back().firstFragment = startOfExonFrag;
//                                    tempGene.inRead.back().fragments.back().lastFragment = endOfExonFrag;
//                                    tempGene.inRead.back().fragments.back().firstPosition = startOfExon;
//                                    tempGene.inRead.back().fragments.back().lastPosition = endOfExon - partLength;
//                                }
//                            }
//                            else if (orderedchosenarray[element].getIndices()[partition] > 0){
//                                tempGene.inExon.push_back(interval());
//                                tempGene.inExon.back().start = startOfExon;
//                                tempGene.inExon.back().end = endOfExon;
//                                tempGene.inRead.back().fragments.push_back(splicedFrag());
//                                tempGene.inRead.back().fragments.back().firstFragment = startOfExonFrag;
//                                tempGene.inRead.back().fragments.back().lastFragment = endOfExonFrag;
//                                tempGene.inRead.back().fragments.back().firstPosition = startOfExon;
//                                tempGene.inRead.back().fragments.back().lastPosition = endOfExon - partLength;
//                                startOfExon = orderedchosenarray[element].getIndices()[partition];
//                                startOfExonFrag = partition;
//                            }
//                        }
//                    }
//                    else{
//                        tempGene.inExon.push_back(interval());
//                        tempGene.inExon.back().start = orderedchosenarray.at(element).getBirth();
//                        tempGene.inExon.back().end = orderedchosenarray.at(element).getLastIndex() + partLength;
//                        tempGene.inRead.back().fragments.push_back(splicedFrag());
//                        tempGene.inRead.back().fragments.back().firstFragment = orderedchosenarray.at(element).getBirthRow();
//                        tempGene.inRead.back().fragments.back().lastFragment = orderedchosenarray.at(element).getLastRow();
//                        tempGene.inRead.back().fragments.back().firstPosition = orderedchosenarray.at(element).getBirth();
//                        tempGene.inRead.back().fragments.back().lastPosition = orderedchosenarray.at(element).getLastIndex();
//                    }
//                }
//                tempGene.isForward = forward_Detection(orderedchosenarray[element].getFlag());
//                bool isCompatibale = false;
//                for (int counterOnInExon = 0; counterOnInExon < tempGene.inExon.size() && !isCompatibale; counterOnInExon++) {
//                    long long startOnInExon = tempGene.inExon.at(counterOnInExon).start;
//                    long long endOnInExon = tempGene.inExon.at(counterOnInExon).end;
//                    for (long long counterOnRealExon = 0; counterOnRealExon < this->iDepth->size(); counterOnRealExon++) {
//                        if ((this->iDepth->at(counterOnRealExon).start < startOnInExon &&
//                             startOnInExon < this->iDepth->at(counterOnRealExon).end) ||
//                            (this->iDepth->at(counterOnRealExon).start < endOnInExon &&
//                             endOnInExon < this->iDepth->at(counterOnRealExon).end)) {
//                                isCompatibale = true;
//                                break;
//                            }
//                    }
//                }

//                if (isCompatibale && !(tempGene.start == tempGene.end && tempGene.start == 0)) {
//                    //add to spliced database and update depth
//                    for (int counterOnGeneTempForCigar = 0; counterOnGeneTempForCigar < tempGene.inRead.size(); counterOnGeneTempForCigar++) {

//                        long long index = tempGene.inRead.at(counterOnGeneTempForCigar).index-1;
//                        /*string cigarForSplice = "";
//                         int lengthOfMatch = 0;
//                         int fragmentAlignedBefore = 0;
//                         long long firstPositionBefore = 0;

//                         for (int counterFragsOnInRead = 0; counterFragsOnInRead < tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size(); counterFragsOnInRead++) {

//                         long long firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).firstPosition;
//                         int firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).firstFragment;
//                         int lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(counterFragsOnInRead).lastFragment;

//                         if (counterFragsOnInRead == 0) {
//                         lengthOfMatch = (lastFragment-firstFragment+1)*partLength-(lastFragment-firstFragment)*d;
//                         cigarForSplice += this->convertNumToStr(lengthOfMatch) + "M";
//                         fragmentAlignedBefore = firstFragment;
//                         } else {
//                         cigarForSplice += this->convertNumToStr(firstPosition-firstPositionBefore-lengthOfMatch) + "N";
//                         cigarForSplice += this->convertNumToStr((firstFragment-fragmentAlignedBefore)*(partLength-d)) + "R";
//                         fragmentAlignedBefore = firstFragment;
//                         lengthOfMatch = (lastFragment-firstFragment+1)*partLength-(lastFragment-firstFragment)*d;
//                         cigarForSplice += this->convertNumToStr(lengthOfMatch) + "M";
//                         }
//                         firstPositionBefore = firstPosition;
//                         }*/

//                        long size = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size();

//                        if (this->reads[index].flag == -1) {
//                            this->reads[index].cigar = tempGene.inRead.at(counterOnGeneTempForCigar).cigar;
//                            this->reads[index].firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstPosition;
//                            this->reads[index].lastPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastPosition;
//                            this->reads[index].firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstFragment;
//                            this->reads[index].lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastFragment;
//                            if (tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size() > 1) {
//                                if (tempGene.isForward) {
//                                    this->reads[index].flag = 7;
//                                } else {
//                                    this->reads[index].flag = 23;
//                                }
//                            }
//                            else{
//                                if (tempGene.isForward) {
//                                    this->reads[index].flag = 6;
//                                } else {
//                                    this->reads[index].flag = 22;
//                                }
//                            }

//                        } else {
//                            this->readsNext->push_back(iReadNext());
//                            this->readsNext->back().cigar = tempGene.inRead.at(counterOnGeneTempForCigar).cigar;
//                            this->readsNext->back().firstPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstPosition;
//                            this->readsNext->back().lastPosition = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastPosition;
//                            this->readsNext->back().firstFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(0).firstFragment;
//                            this->readsNext->back().lastFragment = tempGene.inRead.at(counterOnGeneTempForCigar).fragments.at(size-1).lastFragment;
//                            if (tempGene.inRead.at(counterOnGeneTempForCigar).fragments.size() > 1) {
//                                if (tempGene.isForward) {
//                                    this->readsNext->back().flag = 7;
//                                } else {
//                                    this->readsNext->back().flag = 23;
//                                }
//                            }
//                            else{
//                                if (tempGene.isForward) {
//                                    this->readsNext->back().flag = 6;
//                                } else {
//                                    this->readsNext->back().flag = 22;
//                                }
//                            }
//                            this->readsNext->back().previousIndex = lastNumber;
//                        }

//                        for (int counterOnExonsGenesForAddingToDepth = 0; counterOnExonsGenesForAddingToDepth < tempGene.inExon.size(); counterOnExonsGenesForAddingToDepth++) {
//                            this->addIntervalToDepth(tempGene.inExon.at(counterOnExonsGenesForAddingToDepth).start, tempGene.inExon.at(counterOnExonsGenesForAddingToDepth).end);
//                        }
//                    }
//                }
//                else if (!(tempGene.start == tempGene.end && tempGene.start == 0)){
//                    AddToJenes(tempGene);
//                }

//            }

//        }

//        if(right!=0){
//            {
//                vector<double>().swap(threadReads[thread].scores);
//                vector<double>().swap(threadReads[thread].T_I_vector);
//                vector<double>().swap(threadReads[thread].T_M_vector);
//                vector<int>().swap(threadReads[thread].print_M);
//                vector<int>().swap(threadReads[thread].print_I);
//                vector<int>().swap(threadReads[thread].flags);
//                vector<int>().swap(threadReads[thread].directions);
//                vector<int>().swap(threadReads[thread].seqSizes);                                                                                               vector<long long>().swap(threadReads[thread].loc);
//                vector<string>().swap(threadReads[thread].references);
//            }
//            threadReads[thread].readName=currentFast;
//            threadReads[thread].index = currentFastNumber;
//            threadReads[thread].mainHeader = currentMainHeader;
//            threadReads[thread].Histo = false;
//            threadReads[thread].Checked = false;
//            O_Strings[thread] << currentSeq << "\t" << currentQuality << "\t" << endl;
//            if(!trueFlag)
//                errorCollector[thread]++;
//            for(int j=0;j<scores.size();j++){
//                threadReads[thread].scores.push_back(scores[j]);
//                threadReads[thread].flags.push_back(scoreFlags[j]);
//                if(!Path_I.empty()){
//                    threadReads[thread].seqSizes.push_back(seqSizes[j]);
//                    threadReads[thread].T_I_vector.push_back(Path_I[j]);
//                    threadReads[thread].T_M_vector.push_back(Path_M[j]);
//                }
//                threadReads[thread].print_I.push_back(Print_I[j]);
//                threadReads[thread].print_M.push_back(Print_M[j]);
//                threadReads[thread].loc.push_back(positions[j]);
//                threadReads[thread].references.push_back(Refs[j]);
//                threadReads[thread].directions.push_back(Flags[j]);
//                if(scoreFlags[j]==1 || scoreFlags[j]==0)
//                    threadReads[thread].Histo=true;
//            }
//            if(!threadReads[thread].scores.empty()){
//                threadVectors[thread].push_back(threadReads[thread]);
//            }
//            scores.clear();
//            Path_I.clear();
//            Path_M.clear();
//            Print_M.clear();
//            Print_I.clear();
//            scoreFlags.clear();
//            seqSizes.clear();
//            positions.clear();
//            Refs.clear();
//            Flags.clear();
//        }
//    }
//    // after analysis the new read we should clear all the vectors
//    for(int i=0;i<currentNumOfTables;i++){
//        vector<string>().swap(mySeq[i]);
//        vector<int>().swap(myFlags[i]);
//        vector<string>().swap(myRef[i]);
//        vector<long long>().swap(myTables[2*i]);
//        vector<long long>().swap(myTables[2*i+1]);
//    }
//    {
//        vector<Path>().swap(localPaths);
//    }
//}
//void Assignment::Join2pathes(Path &path1,Path &path2){

//    if (path2.isTruePath()){
//        path1.setTruePath();
//    }
//    for(int i=0;i<path1.getTable().size();i++){
//        if (path2.getTable()[i]!=-1){
//            path1.setTable(i,path2.getTable()[i]);
//        }
//    }
//    for(int i=0;i<path1.getIndices().size();i++){
//        if (path2.getIndices()[i]!=0){
//            path1.setIndexTable(i,path2.getIndices()[i]);
//        }
//    }
//    for (int i=0; i<path1.getCigars().size(); i++) {
//        if (path2.getCigars()[i].size()>0){
//            path1.addToCigars(path2.getCigars()[i], i);
//        }
//    }
//    //***********************************************
//    if (forward_Detection(path1.getFlag())) {
//        long distanceInGenome = path2.getBirth() - path1.getLastIndex();
//        int distanceInRead = (partLength-d)*(path2.getBirthRow() - path1.getLastRow());
//        path1.setTotalCigar(path1.getTotalCigar()+"N"+convertNumToStr(distanceInGenome)+"N"+convertNumToStr(distanceInRead)+"R"+path2.getTotalCigar());
//    }
//    else{
//        long distanceInGenome = path1.getLastIndex() - path2.getBirth();
//        int distanceInRead = (partLength-d)*(path2.getBirthRow() - path1.getLastRow());
//        path1.setTotalCigar(path1.getTotalCigar()+"N"+convertNumToStr(distanceInGenome)+"N"+convertNumToStr(distanceInRead)+"R"+path2.getTotalCigar());
//    }
//    //***********************************************

//    path1.addScore(path2.getScore());
//    path1.joining();
//}

//string Assignment::Glue_cigars(string first, string second) {
//    string collect;
//    string cut_first, cut_second;
//    if(isdigit(first[first.length()-1]) && isdigit(second[0])){
//        string first_digits;
//        string second_digits;
//        int i,j;
//        for( i=(int)first.length()-1;i>=0;i--){
//            if(isdigit(first[i])){
//                first_digits += first[i];
//            }
//            else{
//                break;
//            }
//        }
//        for( j=0;j<second.length();j++){
//            if(isdigit(second[j])){
//                second_digits += second[j];
//            }
//            else{
//                break;
//            }
//        }
//        cut_first = first.substr(0,i+1);
//        cut_second = second.substr(j);
//        collect = cut_first + NumToStr(atoi(first_digits.c_str())+atoi(second_digits.c_str()));
//        collect += cut_second;
//    }
//    else{
//        collect = first+second;
//    }
//    return collect;
//}

//long long Assignment::Calculate(long long filledAll, int numberOfThreads) {
//    long long step;
//    step = floor(filledAll / numberOfThreads);
//    return step;
//}

//int Assignment::checkNumOfTables(int partLength, int V_L) {
//    int numOfTables;
//    if(V_L < partLength){
//        numOfTables=1;
//        return numOfTables;
//    }
//    else{
//        numOfTables=(int)(floor((V_L-partLength)/(partLength-d)))+1;
//        return numOfTables;
//    }
//}

//int Assignment::flagDetection(long long ReportedFlag) {
//    int flag=0;
//    if (ReportedFlag==4){
//        // this flag means the read is not aligned
//        flag=4;
//    }
//    else if (ReportedFlag==0){
//        // this flag means the read is forward
//        flag = 0;
//    }
//    else if(ReportedFlag%256==0){
//        // this flag means this read is aligned forwardly
//        // and has more than one option for alignment
//        flag = 0;
//    }
//    else if(ReportedFlag==16){
//        // this flag means the read is reversed
//        flag = 16;
//    }
//    else if((ReportedFlag-256)!=0){
//        // this flag means this read is reversed
//        // and has more than one option for alignment
//        flag = 16;
//    }
//    return flag;
//}

//int Assignment::numeralCounter(int number) {
//    return floor(log10(number))+1;
//}

//void Assignment::emptyInterval_Detector(vector<bool> &quantArray, vector<long> *emptyIntervals, long long current, int currentLength) {
//    bool flag=false,flag1=false;//**********************************************************
//    for(int i=0;i<currentLength;i++){
//        if(quantArray[i]){
//            if(flag){
//                emptyIntervals->push_back(i-1);
//            }
//            flag = false;
//            flag1 = true;//**********************************************************
//        }
//        else{
//            if(!flag && flag1){//**********************************************************
//                flag = true;
//                emptyIntervals->push_back(i);
//            }
//        }
//    }
//}

//bool Assignment::histo_Permission() {
//    long counter=0;
//    for(int i=0;i<allReads.size();i++){
//        for(int j=0;j<allReads[i].scores.size();j++){
//            if(allReads[i].flags[j]==1 || allReads[i].flags[j]==0)
//                counter++;
//        }
//    }
//    if(counter>=100){
//        return true;
//    }
//    else{
//        return false;
//    }
//}

//double Assignment::findInitMax() {
//    bool flag=false;
//    double max = 0.0;
//    for(int i=0;i<allReads.size();i++){
//        for(int j=0;j<allReads[i].scores.size();j++){
//            if(allReads[i].flags[j]==0 || allReads[i].flags[j]==1){
//                max = allReads[i].scores[j];
//                flag=true;
//                break;
//            }
//        }
//        if(flag)
//            break;
//    }
//    return max;
//}

//string Assignment::indexMaker(int index) {
//    int numberOfDigits;
//    int numberOfIndexDigits;
//    std::ostringstream str;
//    if(index==0){
//        for(int i=0;i<numeralCounter(numberOfTables);i++){
//            str << "0";
//        }
//    }
//    else{
//        numberOfDigits = numeralCounter(numberOfTables);
//        numberOfIndexDigits = numeralCounter(index);
//        for(int i=0;i<(numberOfDigits-numberOfIndexDigits);i++){
//            str << "0";
//        }
//        str << index;
//    }
//    return str.str();
//}

//string Assignment::complement(char character) {
//    string rev;
//    if(character == 'A'){
//        rev="T";
//    }
//    if(character == 'C'){
//        rev="G";
//    }
//    if(character=='G'){
//        rev="C";
//    }
//    if(character=='T'){
//        rev="A";
//    }
//    if(character=='a'){
//        rev="t";
//    }
//    if(character=='c'){
//        rev="g";
//    }
//    if(character=='g'){
//        rev="c";
//    }
//    if(character=='t'){
//        rev="a";
//    }
//    if(character=='N' || character=='n'){
//        rev="N";
//    }
//    return rev;
//}

//std::string Assignment::accessGenome(long long index, int flag, long long length, long long refPos) {
//    struct stat sb;
//    off_t len;
//    char *p;
//    int fd;
//    ifstream ifstr(genomeName.c_str());
//    string firstLine;
//    getline(ifstr, firstLine);
//    ifstr.close();
//    index += refPos+(index/characterPerRow);
//    fd = open(genomeName.c_str(), O_RDONLY);
//    if (fd == -1)
//        perror ("open");

//    if (fstat (fd, &sb) == -1)
//        perror ("fstat");

//    p = (char *)mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
//    if (close (fd) == -1)
//        perror ("close");

//    string temp = "";
//    //index--;
//    if(forward_Detection(flag)){
//        for (len = index; temp.size() < length; ){
//            if (p[len] != '\n')
//                temp += p[len];
//            ++len;
//        }
//    }

//    else if(reverse_Detection(flag)){
//        for (len = (index+length); temp.size() < (length); ){
//            temp += complement(p[len]);
//            --len;
//        }
//    }
//    if (munmap (p, sb.st_size) == -1)
//        perror ("munmap");
//    return temp;
//}

//bool myComparison(Path first, Path second) {
//    return first.getScore() < second.getScore();
//}

//bool vectorComparison(Read_Vector first ,Read_Vector second) {
//    return first.size < second.size;
//}

//bool littleSort(littleStruct first, littleStruct second) {
//    return first.score < second.score;
//}

//string Assignment::finalString(vector<string> &final_Cigars, vector<string> &final_empty_Cigars, long long first_Part){
//    string outStr;
//    int empty_pointer=0, filled_pointer=0;
//    bool readFilled=false, readEmpty=false;
//    if(local_Flag){
//        if(first_Part==0){
//            outStr = final_empty_Cigars[empty_pointer];
//            empty_pointer++;
//            readFilled=true;
//        }
//        else{
//            outStr = final_Cigars[filled_pointer];
//            filled_pointer++;
//            readEmpty=true;
//        }
//        while(empty_pointer<final_empty_Cigars.size() || filled_pointer<final_Cigars.size()){
//            if(readFilled){
//                if(final_Cigars.empty() || filled_pointer==(final_Cigars.size())){
//                    readEmpty=true;
//                    readFilled=false;
//                }
//                else{
//                    outStr = Glue_cigars(outStr,final_Cigars[filled_pointer]);
//                    readFilled=false;
//                    readEmpty=true;
//                    filled_pointer++;
//                }
//            }
//            if(readEmpty){
//                if(final_empty_Cigars.empty() || empty_pointer==(final_empty_Cigars.size())){
//                    readEmpty=false;
//                    readFilled=true;
//                }
//                else{
//                    outStr = Glue_cigars(outStr,final_empty_Cigars[empty_pointer]);
//                    readEmpty=false;
//                    readFilled=true;
//                    empty_pointer++;
//                }
//            }
//        }
//        return outStr;
//    }

//    else{
//        for(int i=0;i<final_Cigars.size();i++){
//            outStr = Glue_cigars(outStr,final_Cigars[i]);
//        }
//        return outStr;
//    }
//}

//string Assignment::L_MAX_Gate(int L_max, int currentLength, string seq) {
//    if(currentLength<=L_max){
//        return seq;
//    }
//    else{
//        string temp=seq.substr(0,L_max);
//        return temp;
//    }
//}

//int Assignment::get_mD_Length(string mD, int &M_cigar, int &I_cigar) {
//    bool misMatch=false,indel=false,match=false;
//    int temp=0;
//    string str;
//    for(int i=0;i<mD.length();i++){
//        if(isalpha(mD[i])){
//            if(indel){
//                temp++;
//                I_cigar++;
//                continue;
//            }
//            else{
//                if(match){
//                    int matchNumber = atoi(str.c_str());
//                    temp+=matchNumber;
//                    str.clear();
//                    match=false;
//                }
//                misMatch=true;
//                M_cigar++;
//                temp++;
//                continue;
//            }
//        }
//        else if(!isalnum(mD[i])){
//            indel=true;
//            if(match){
//                int matchNumber = atoi(str.c_str());
//                temp+=matchNumber;
//                str.clear();
//                match=false;
//            }
//            else if(misMatch){
//                misMatch=false;
//            }
//            continue;
//        }
//        else{
//            if(misMatch){
//                misMatch=false;
//            }
//            else if(indel){
//                indel=false;
//            }
//            match=true;
//            str+=mD[i];
//            continue;
//        }
//    }
//    if(match){
//        int matchNumber = atoi(str.c_str());
//        temp+=matchNumber;
//        str.clear();
//    }
//    return temp;
//}

//void Assignment::Save_Properties(vector<string> &samVector,vector< vector<int> > &myFlags,vector< vector<string> > &mySeq ,vector< vector<long long> > &myTables,vector< vector<string> > &myRef,int row, int mode) {
//    string scoreFlag,cigarFlag, currentCigar, currentRef, tempString;
//    int tempFlag, currentScore=0;
//    long long currentAlign=0;
//    bool NullFlag=false;
//    cigarFlag = "MD:Z:";
//    if(mode==BOWTIE2)
//        scoreFlag = "AS:i:";
//    else if(mode==BOWTIE1)
//        scoreFlag = "NM:i:";
//    for (int i=0;i<samVector.size();i++) {
//        if(i==1){
//            //Reading the flag of the Read
//            // flags 0 and 256 mean reverse alignment
//            // flags 16 and 272 mean forward alignment
//            // flag 4 means non-aligned read
//            // other flags will be added
//            if(myFlags[row].size() < (myDepth)){
//                tempString = samVector[i];
//                tempFlag = flagDetection(atoll(tempString.c_str()));
//                if(tempFlag==4){
//                    NullFlag = true;
//                    currentAlign=0;
//                    currentCigar="null";
//                    currentScore=0;
//                }
//            }
//        }
//        if(!NullFlag){
//            if(i==3){
//                //Saving the position of the Aligned Read
//                if(myTables[2*row+1].size()<(myDepth+1)){
//                    tempString = samVector[i];
//                    currentAlign = 0;
//                    if (tempString != "*" && tempString != "0") {
//                        currentAlign = atoll(tempString.c_str());
//                    }
//                }
//            }
//            if(i==2){
//                //Saving the reference
//                if(myRef[row].size()<(myDepth)){
//                    tempString = samVector[i];
//                    if (tempString != "*" && tempString != "0") {
//                        currentRef = tempString.c_str();
//                    }
//                }
//            }
//            if(samVector[i].find(scoreFlag)!=std::string::npos){
//                //Saving the Score of the Read
//                if(myTables[2*row].size()<(myDepth)){
//                    tempString = samVector[i];
//                    currentScore = 0;
//                    if(tempString !="*" && tempString != "0"){
//                        tempString  = tempString.substr(scoreFlag.length());
//                        currentScore = atoi(tempString.c_str());
//                        if(mode==BOWTIE1)
//                            currentScore = MATCH*partLength - (MATCH-MISPEN)*currentScore;
//                    }
//                }
//            }
//            if(samVector[i].find(cigarFlag)!=std::string::npos){
//                //Saving cigar string if it exists
//                if(mySeq[row].size()<(myDepth)){
//                    tempString = samVector[i];
//                    tempString  = tempString.substr(cigarFlag.length());
//                    currentCigar = tempString;
//                }
//            }
//        }
//    }
//    if(!NullFlag){
//        if(myTables[2*row+1].size()<(myDepth+1)){
//            myFlags[row].push_back(tempFlag);
//            myTables[2*row+1].push_back(currentAlign);
//            myTables[2*row].push_back(currentScore);
//            mySeq[row].push_back(currentCigar);
//            myRef[row].push_back(currentRef);
//        }
//    }
//    return;
//}

//void Assignment::gen_Paths(int currentNumOfTables, int currentLength, vector< vector <long long> > &myTables, vector< vector <int> > &myFlags, vector< vector<string> > &myRef, vector <vector <string> > &mySeq, vector<long long> &truePath, vector<Path> &localPaths, vector <long long> &faileds) {
//    vector <Path> minorLocalPaths;
//    Path currentPath;
//    for (int j = 0; j < currentNumOfTables; j++){
//        long long counter = 0;
//        //Initializing all paths by creating Path objects in first row
//        if (!myTables[0].empty() && j == 0 && myTables[1][counter+1]!=0 && myTables[1][counter+1]!=-1){
//            truePath[0] = myTables[1][counter];
//            while ((counter+1)<myTables[1].size() && myTables[1][counter + 1] != 0 && myTables[1][counter+1]!=-1){
//                localPaths.push_back(Path(0,myTables[1][counter + 1], currentNumOfTables,myFlags[0][counter],partLength, d, myRef[0][counter], currentLength));
//                localPaths[localPaths.size()-1].setTable(0, 0);
//                localPaths[localPaths.size()-1].addScore(myTables[0][counter]);
//                localPaths[localPaths.size()-1].setIndexTable(0,myTables[1][counter + 1]);
//                localPaths[localPaths.size()-1].addToCigars(mySeq[0][counter],0);
//                long long pos_dam1=localPaths[localPaths.size()-1].getBirth();

//                if(reverse_Detection(localPaths[localPaths.size()-1].getFlag())){
//                    pos_dam1=pos_dam1-currentLength;
//                }
//                if (std::abs(truePath[0]-pos_dam1) <= errorInIndex) {
//                    localPaths[localPaths.size()-1].setTruePath();
//                }
//                counter++;
//            }
//        }
//        else if((j==0 && myTables[0].empty()) || (j==0 && myTables[1][counter+1]==0)||(j==0 && myTables[1][counter+1]==-1)){
//            // make no new paths
//            continue;
//        }
//        else if(j!=0){

//            if(j!=currentNumOfTables-1 &&(myTables[2*j].empty() || myTables[2*j+1][counter+1]==0 || myTables[2*j+1][counter+1]==-1)){
//                // make no new paths in next lines
//                continue;
//            }
//            else if(!myTables[2*j].empty() && myTables[2*j+1][counter+1]!=0 && myTables[2*j+1][counter+1]!=-1){
//                truePath[j] = myTables[2 * j + 1][0];
//                while (myTables[2 * j + 1][counter + 1] != -1 && myTables[2 * j + 1][counter + 1] != 0 && (counter+1)<myTables[2*j+1].size()) {
//                    int flag = 0;
//                    for (int p = 0; p < localPaths.size(); p++) {
//                        long long equivalent = 0;
//                        int PathFlag = 0;
//                        int PartFlag = 0;
//                        string PartRef, PathRef;
//                        currentPath = localPaths[p];
//                        PathFlag = currentPath.getFlag();
//                        PathRef = currentPath.getRef();
//                        PartFlag = myFlags[j][counter];
//                        PartRef = myRef[j][counter];
//                        if(forward_Detection(PathFlag)){
//                            // this means the read is aligned forwardly
//                            equivalent = currentPath.getBirth() + (j-currentPath.getBirthRow())*(partLength-d);
//                        }
//                        else if(reverse_Detection(PathFlag)){
//                            // this means the read is aligned reservely
//                            equivalent = currentPath.getBirth()-(j - currentPath.getBirthRow())*(partLength-d);
//                        }
//                        if ((std::abs(myTables[2*j+1][counter+1]-equivalent)<=(j*error))&&(PartFlag==PathFlag)&&(PartRef==PathRef)) {
//                            localPaths[p].addScore(myTables[2 * j][counter]);
//                            localPaths[p].setTable(j, j);
//                            localPaths[p].setIndexTable(j,myTables[2 * j + 1][counter + 1]);
//                            localPaths[p].addToCigars(mySeq[j][counter],j);
//                            flag = 1;
//                        }
//                    }
//                    if (flag == 0) {
//                        minorLocalPaths.push_back(Path(j,myTables[2 * j + 1][counter + 1], currentNumOfTables,myFlags[j][counter],partLength, d,myRef[j][counter], currentLength));
//                        minorLocalPaths[minorLocalPaths.size()-1].setTable(j, j);
//                        minorLocalPaths[minorLocalPaths.size()-1].addScore(myTables[2*j][counter]);
//                        minorLocalPaths[minorLocalPaths.size()-1].setIndexTable(j,myTables[2 * j + 1][counter + 1]);
//                        minorLocalPaths[minorLocalPaths.size()-1].addToCigars(mySeq[j][counter],j);
//                        long long pos_dam=minorLocalPaths[minorLocalPaths.size()-1].getBirth();

//                        if(reverse_Detection(minorLocalPaths[minorLocalPaths.size()-1].getFlag())){
//                            pos_dam=pos_dam-currentLength+(j*(partLength-d));
//                        }
//                        if (std::abs(truePath[j]-pos_dam)<=errorInIndex) {
//                            minorLocalPaths[minorLocalPaths.size()-1].setTruePath();
//                        }
//                    }
//                    counter++;
//                }
//                // copying the minor array into global array of Paths
//                if (!minorLocalPaths.empty()) {
//                    localPaths.insert(localPaths.end(), minorLocalPaths.begin(), minorLocalPaths.end());
//                    {
//                        vector<Path>().swap(minorLocalPaths);
//                    }
//                }
//            }
//            else if(j==currentNumOfTables-1 && localPaths.size()==0){
//                // do sth for failed indices;
//                faileds.push_back(myTables[1][counter]);
//                for(int i=0;i<currentNumOfTables;i++){
//                    vector<long long>().swap(myTables[2*i]);
//                    vector<long long>().swap(myTables[2*i+1]);
//                    vector<string>().swap(myRef[i]);
//                    vector<string>().swap(mySeq[i]);
//                    vector<int>().swap(myFlags[i]);
//                }
//                vector<Path>().swap(localPaths);
//                continue;
//            }
//        }
//    }
//    return;
//}

//void Assignment::Filters(vector<Path> &localPaths, vector<Path> &tempArray, vector<Path> &chosenArray) {
//    std::sort(localPaths.begin(), localPaths.end(), myComparison);
//    //Choosing the scores
//    for (int score = (int)localPaths.size()-1;score!=-1; score--) {
//        //testing << "before filtering: " << localPaths[score].getMain() << " score: " << localPaths[score].getScore() << endl;
//        // Filtering the scores by basicThreshold
//        if(localPaths[localPaths.size()-1].getScore()<0){
//            //if(localPaths[score].getScore() >= 2*(localPaths[localPaths.size() - 1].getScore())){
//            tempArray.push_back(localPaths[score]);
//            //}
//        }
//        else{
//            //if (localPaths[score].getScore() >= (localPaths[localPaths.size() - 1].getScore() / 2)) {
//            tempArray.push_back(localPaths[score]);
//            //}
//        }
//        if (localPaths[score].getScore() == 0) {
//            break;
//        }
//    }
//    // Filtering by computation of differences and consecutive threshold
//    for (int t = 0; t < tempArray.size(); t++) {
//        if ((t + 1) == tempArray.size()) {
//            chosenArray.push_back(tempArray[t]);
//            break;
//        }
//        chosenArray.push_back(tempArray[t]);
//        /*double proportion = abs((double)((double)(tempArray[t].getScore() - tempArray[t + 1].getScore()))
//         / (double) (tempArray[t + 1].getScore()));
//         if (proportion > consecutiveThreshold) {
//         break;
//         }*/
//    }
//    return;
//}

//void Assignment::merge_Method(vector<int> &mergeParts, int currentNumOfTables, int completing, vector<Path> &chosenArray,vector<bool> &quantArray, vector<string> &final_Cigars, int &M_cigar, int &I_cigar) {
//    //forward
//    {
//        vector<int>().swap(mergeParts);
//    }
//    long long current_Start=0;
//    long long current_End=0;
//    long long pointer_Start=0;
//    long long pointer_End=0;
//    long long startIndexB=0;
//    bool terminate=false;
//    bool save=false;
//    bool first=true;
//    int savePart = 0;
//    int tempFlag = chosenArray[completing].getFlag();
//    for(int part=0;part<currentNumOfTables;part++){
//        if(chosenArray[completing].getIndices()[part]!=0 && !terminate){
//            startIndexB=part*(partLength-d);
//            if(save){
//                mergeParts.push_back(savePart);
//                save=false;
//                savePart=0;
//            }
//            current_Start = chosenArray[completing].getIndices()[part];
//            int fake1, fake2;
//            if(forward_Detection(tempFlag)){
//                current_End = current_Start + get_mD_Length(chosenArray[completing].getCigars()[part], fake1, fake2)-1;
//            }
//            else if(reverse_Detection(tempFlag)){
//                current_End = current_Start - get_mD_Length(chosenArray[completing].getCigars()[part], fake1, fake2)-1;
//            }
//            if((current_Start<=pointer_End && current_Start>=pointer_Start && forward_Detection(tempFlag)) || first || (current_Start>=pointer_End && current_Start<=pointer_Start && reverse_Detection(tempFlag))){
//                if(first){
//                    first=false;
//                }
//                pointer_End=current_End;
//                pointer_Start=current_Start;
//                mergeParts.push_back(part);
//                if(part==currentNumOfTables-1){
//                    terminate=true;
//                }
//            }
//            else{
//                terminate=true;
//                pointer_Start=current_Start;
//                pointer_End=current_End;
//                savePart=part;
//                save=true;
//            }
//            for(int i=0;i<partLength;i++){
//                quantArray[startIndexB+i]=true;
//            }
//        }
//        else if(chosenArray[completing].getIndices()[part]==0 && !terminate){
//            if(part==currentNumOfTables-1){
//                terminate=true;
//            }
//            else{
//                if(save){
//                    mergeParts.push_back(savePart);
//                    save=false;
//                    savePart=0;
//                }
//                continue;
//            }
//        }
//        if(terminate){
//            std::string mergeStr,mergeMD;
//            if(mergeParts.size()==1){
//                get_mD_Length(chosenArray[completing].getCigars()[mergeParts[0]], M_cigar,I_cigar);
//                final_Cigars.push_back(chosenArray[completing].getCigars()[mergeParts[0]]);
//            }
//            if(mergeParts.size()>1){
//                std::string temp;
//                long long first, second;
//                for(int i=0;i<mergeParts.size()-1;i++){
//                    temp=mD2String(chosenArray[completing].getCigars()[mergeParts[i]]);
//                    first = chosenArray[completing].getIndices()[mergeParts[i]];
//                    second = chosenArray[completing].getIndices()[mergeParts[i+1]];
//                    mergeStr+=slice(0,(int)abs(second-first),temp);
//                }
//                int size = (int) mergeParts.size();
//                temp=mD2String(chosenArray[completing].getCigars()[mergeParts[size-1]]);
//                mergeMD = string2mD(mergeStr);
//                mergeMD = Glue_cigars(mergeMD,chosenArray[completing].getCigars()[mergeParts[size-1]]);
//                get_mD_Length(mergeMD,M_cigar,I_cigar);
//                final_Cigars.push_back(mergeMD);
//            }
//            if(part==currentNumOfTables-1 && save){
//                get_mD_Length(chosenArray[completing].getCigars()[savePart],M_cigar, I_cigar);
//                final_Cigars.push_back(chosenArray[completing].getCigars()[savePart]);
//            }
//            terminate=false;
//            {
//                vector<int>().swap(mergeParts);
//            }
//        }
//    }
//    return;
//}

//void Assignment::fill_The_Gaps(vector<Path> &chosenArray, int &completing, vector<long> &emptyIntervals, string &currentSeq, double &IndelShift, int &tempFlag, double &totalScore, int &seqSize, int &localCalling, bool &hasRemained, vector<string> &final_empty_Cigars, int &Total_M, int &Total_I) {
//    long long startIndexA , startIndexB = 0, startIndexA_end, startIndexB_end = 0, refPos;
//    string seQA,seQB;
//    string seqConstantB="";
//    LocalAligner *local = NULL;
//    for(int j=0;j<emptyIntervals.size()-1;){
//        refPos = multiFasta[chosenArray[completing].getRef()];
//        if(forward_Detection(tempFlag)){
//            startIndexA = chosenArray[completing].getMain()+emptyIntervals[j]-1;
//            startIndexA_end = chosenArray[completing].getMain()+emptyIntervals[j+1]-1;
//            startIndexB = emptyIntervals[j];
//            startIndexB_end = emptyIntervals[j+1];
//            seQB = currentSeq.substr(startIndexB,startIndexB_end - startIndexB);
//            seqConstantB.clear();
//            for(int i=0;i<(int)(IndelShift*(abs(startIndexB_end - startIndexB)));i++){
//                seqConstantB += "A";
//            }
//            seQB += seqConstantB;
//            seQA = accessGenome(startIndexA-(int)(IndelShift*(abs(startIndexB_end - startIndexB))),tempFlag,startIndexA_end-startIndexA+2*(int)(IndelShift*(abs(startIndexB_end - startIndexB))), refPos);
//        }
//        else if(reverse_Detection(tempFlag)){
//            startIndexA = chosenArray[completing].getMain()-emptyIntervals[j]-1;
//            startIndexA_end = chosenArray[completing].getMain()-emptyIntervals[j+1]-1;
//            startIndexB = emptyIntervals[j];
//            startIndexB_end = emptyIntervals[j+1];
//            seQB = currentSeq.substr(startIndexB,(abs(startIndexB_end - startIndexB)));
//            seqConstantB.clear();
//            for(int i=0;i<(int)(IndelShift*(abs(startIndexB_end - startIndexB)));i++){
//                seqConstantB += "A";
//            }
//            seQB += seqConstantB;
//            seQA = accessGenome(startIndexA_end-(int)(IndelShift*(abs(startIndexB_end - startIndexB))),tempFlag,abs(startIndexA_end-startIndexA)+2*(int)(IndelShift*(abs(startIndexB_end - startIndexB))), refPos);
//        }
//        local = new LocalAligner(seQB,seQA,GAPPEN,MISPEN,MATCH,(int)(IndelShift*(abs(startIndexB_end - startIndexB))),(int)(K_factor*(abs(startIndexB_end - startIndexB))));
//        local->process();
//        local->backtrack();
//        local->produceCigar();
//        totalScore = totalScore+local->mScore;
//        Total_M+=local->totalMismatch;
//        Total_I+=local->totalGap;
//        seqSize+=seQB.size()-(int)(IndelShift*(abs(startIndexB_end - startIndexB)));
//        final_empty_Cigars.push_back(local->cigar);
//        //delete local;
//        localCalling++;
//        j+=2;
//    }//**********************************************************avalo akhar nabas por she.
//    /*if(hasRemained){
//     refPos = multiFasta[chosenArray[completing].getRef()];
//     if(forward_Detection(tempFlag)){
//     startIndexA = chosenArray[completing].getMain()+emptyIntervals[emptyIntervals.size()-1]-1;
//     startIndexB = emptyIntervals[emptyIntervals.size()-1];
//     seQB = currentSeq.substr(startIndexB);
//     seqConstantB.clear();
//     for(int i=0;i<(int)(IndelShift*(seQB.length()));i++){
//     seqConstantB += "A";
//     }
//     seQB += seqConstantB;
//     seQA = accessGenome(startIndexA-(int)(IndelShift*(seQB.length())), tempFlag,seQB.size()+(int)(IndelShift*(seQB.length())), refPos);
//     }
//     else if(reverse_Detection(tempFlag)){
//     startIndexA = chosenArray[completing].getMain()-emptyIntervals[emptyIntervals.size()-1]-1;
//     startIndexB = emptyIntervals[emptyIntervals.size()-1];
//     seQB = currentSeq.substr(startIndexB);
//     seqConstantB.clear();
//     for(int i=0;i<(int)(IndelShift*(seQB.length()));i++){
//     seqConstantB += "A";
//     }
//     seQB += seqConstantB;
//     seQA = accessGenome(startIndexA-seQB.size(),tempFlag,seQB.size()+(int)(IndelShift*(seQB.length())),refPos);
//     }
//     local = new LocalAligner(seQB,seQA,GAPPEN,MISPEN,MATCH,(int)(IndelShift*(seQB.length())),(int)(K_factor*(seQB.length())));
//     local->process();
//     local->backtrack();
//     local->produceCigar();
//     totalScore = totalScore + local->mScore;
//     Total_M+=local->totalMismatch;
//     Total_I+=local->totalGap;
//     seqSize+=seQB.size()-(int)(IndelShift*currentSeq.substr(startIndexB).size());
//     final_empty_Cigars.push_back(local->cigar);
//     delete local;
//     localCalling++;
//     }*/
//    return;
//}

//bool Assignment::forward_Detection(int tempFlag) {
//    if(tempFlag==0 || tempFlag==256)
//        return true;
//    else{
//        return false;
//    }
//}

//bool Assignment::reverse_Detection(int tempFlag) {
//    if(tempFlag==16 || tempFlag==272){
//        return true;
//    }
//    else{
//        return false;
//    }
//}

//void Assignment::remain_Local(int V_L, Path &currentPath, int &remain_M, int &remain_I, int &remain_Score, string &localStr, string seq) {
//    LocalAligner *local = NULL;
//    string read_Seq, path_Seq, seqConstant;
//    long long startPos_Read, startPos_Path;
//    startPos_Read = (checkNumOfTables(partLength,V_L)-1)*(partLength-d)+partLength;
//    read_Seq = seq.substr(startPos_Read);
//    startPos_Path = startPos_Read+currentPath.getMain();
//    seqConstant.clear();
//    for(int i=0;i<(int)(IndelShift*(read_Seq.length()));i++){
//        seqConstant += "A";
//    }
//    read_Seq += seqConstant;
//    path_Seq = accessGenome(startPos_Path-(int)(IndelShift*(read_Seq.length())), currentPath.getFlag(),read_Seq.size()+(int)(IndelShift*(read_Seq.length())), multiFasta[currentPath.getRef()]);
//    local = new LocalAligner(read_Seq,path_Seq,GAPPEN,MISPEN,MATCH,(int)(IndelShift*(read_Seq.length())),(int)(K_factor*read_Seq.length()));
//    local->process();
//    local->backtrack();
//    local->produceCigar();
//    remain_M=local->totalMismatch;
//    remain_I=local->totalGap;
//    remain_Score=local->mScore;
//    localStr=local->cigar;
//    //delete local;
//    return;
//}

//long long Assignment::NextChar(std::string fileData, long long pos,string str) {
//    long long position;
//    position = pos;
//    while (fileData.substr(position,1) != str) {
//        position++;
//    }
//    return position;
//}

//vector<string> Assignment::split(const string & str, const string & delimiters) {
//    vector<string> v;
//    string::size_type start = 0;
//    string::size_type pos = str.find_first_of(delimiters, start);
//    while (pos != string::npos) {
//        if (pos != start) { // ignore empty tokens
//            const string temp = str.substr(start, pos - start);
//            v.push_back(temp);
//        }
//        start = pos + 1;
//        pos = str.find_first_of(delimiters, start);
//    }
//    if (start < str.length()) { // ignore trailing delimiter
//        const string temp = str.substr(start, str.length() - start);
//        v.push_back(temp);
//    } // add what's left of the string
//    return v;
//}

//std::string Assignment::NumToStr(long long number) {
//    std::string str;
//    std::ostringstream ss;
//    ss << number;
//    return ss.str();
//}

//string Assignment::mD2String(string mD) {
//    bool misMatch=false,indel=false,match=false;
//    string str;
//    string temp;
//    for(int i=0;i<mD.length();i++){
//        if(isalpha(mD[i])){
//            if(indel){
//                temp+="^";
//                temp+=mD[i];
//                continue;
//            }
//            else{
//                if(match){
//                    int matchNumber = atoi(temp.c_str());
//                    temp.clear();
//                    for(int i=0;i<matchNumber;i++){
//                        str+='X';
//                    }
//                    match=false;
//                }
//                misMatch=true;
//                temp+=mD[i];
//                continue;
//            }
//        }
//        else if(!isalnum(mD[i])){
//            indel=true;
//            if(match){
//                int matchNumber = atoi(temp.c_str());
//                for(int i=0;i<matchNumber;i++){
//                    str+='X';
//                }
//                match=false;
//            }
//            else if(misMatch){
//                str+=temp;
//                misMatch=false;
//            }
//            temp.clear();
//            continue;
//        }
//        else{
//            if(misMatch){
//                str+=temp;
//                misMatch=false;
//                temp.clear();
//            }
//            else if(indel){
//                str+=temp;
//                indel=false;
//                temp.clear();
//            }
//            match=true;
//            temp+=mD[i];
//            continue;
//        }
//    }
//    if(match){
//        int matchNumber = atoi(temp.c_str());
//        for(int i=0;i<matchNumber;i++){
//            str+='X';
//        }
//    }
//    else if(misMatch){
//        str+=temp;
//    }
//    else if(indel){
//        str+=temp;
//    }
//    return str;
//}

//string Assignment::slice(int begin, int length,string str){
//    int counter=0;
//    int i=0;
//    bool indel=false, misOrMatch=false;
//    string outString;
//    //find the right place to begin
//    while(counter<begin){
//        if(!isalnum(str[i])){
//            i=i+2;
//            counter++;
//            continue;
//        }
//        else if(isalpha(str[i])){
//            i++;
//            counter++;
//            continue;
//        }
//    }
//    // collect all the character including indels
//    counter=0;
//    while(counter<length){
//        if(!isalnum(str[i])){
//            if(!indel && !misOrMatch){
//                indel=true;
//                outString+="^";
//                outString+=str[i+1];
//                counter++;
//            }
//            else if(indel){
//                outString+=str[i+1];
//                counter++;
//            }
//            else if(misOrMatch){
//                indel=true;
//                outString+="^";
//                outString+=str[i+1];
//                counter++;
//            }
//            i=i+2;
//            continue;
//        }
//        else{
//            misOrMatch=true;
//            if(indel)
//                indel=false;
//            outString+=str[i];
//            i++;
//            counter++;
//            continue;
//        }
//    }
//    return outString;
//}

//string Assignment::string2mD(string str) {
//    string outString;
//    int i=0, counter=0;
//    bool match=false;
//    bool indel=false;
//    char sign='X';
//    while(i<str.length()){
//        if((int)str[i]==(int)(sign)){
//            if(indel){
//                indel=false;
//            }
//            if(!match)
//                match=true;
//            counter++;
//        }
//        else{
//            if(match){
//                outString+=NumToStr(counter);
//                match=false;
//                counter=0;
//            }
//            if(!indel && !isalnum(str[i])){
//                indel=true;
//                outString+=str[i];
//            }
//            if(indel && isalpha(str[i]))
//                outString+=str[i];
//            if(!match && !indel)
//                outString+=str[i];
//        }
//        i++;
//    }
//    if(match){
//        outString+=NumToStr(counter);
//    }
//    return outString;
//}

//string Assignment::mD2Cigar(string mD){
//    bool misMatch=false,indel=false,match=false;
//    int I_Counter=0, M_Counter=0;
//    string str;
//    string outStr="";
//    for(int i=0;i<mD.length();i++){
//        if(isalpha(mD[i])){
//            if(indel){
//                I_Counter++;
//                continue;
//            }
//            else{
//                if(match){
//                    int matchNumber = atoi(str.c_str());
//                    M_Counter+=matchNumber;
//                    match=false;
//                    str.clear();
//                }
//                misMatch=true;
//                continue;
//            }
//        }
//        else if(!isalnum(mD[i])){
//            indel=true;
//            if(match){
//                int matchNumber = atoi(str.c_str());
//                M_Counter+=matchNumber;
//                match=false;
//                str.clear();
//            }
//            else if(misMatch){
//                misMatch=false;
//            }
//            outStr+=NumToStr(M_Counter) + "M";
//            M_Counter=0;
//            continue;
//        }
//        else{
//            if(misMatch){
//                misMatch=false;
//            }
//            else if(indel){
//                indel=false;
//                outStr+=NumToStr(I_Counter)+"I";
//                I_Counter=0;
//            }
//            match=true;
//            str+=mD[i];
//            continue;
//        }
//    }
//    if(match){
//        int matchNumber = atoi(str.c_str());
//        M_Counter+=matchNumber;
//        outStr+=NumToStr(M_Counter)+"M";
//        str.clear();
//    }
//    return outStr;
//}

//void Assignment::readFasta(string filename,std::map<string ,long long> &multiFasta) {
//    ifstream infile(filename.c_str());
//    string line, header;
//    long long faSize = 0, genomeSize = 0;
//    map <string,long long> multiFA;
//    ofstream headerFile(outputDir+"headerFile.txt");
//    headerFile << "@HD\tVN:1.0\tSO:unsorted\n";
//    while(getline(infile, line))
//    {
//        faSize += line.size()+1;
//        if (line[0]=='>'){
//            multiFA[line.substr(1)] = faSize;
//            if (genomeSize){
//                headerFile << "@SQ\tSN:" << header << "\tLN:" << genomeSize << endl;
//                genomeSize = 0;
//                header = line.substr(1);
//            }
//            else
//                header = line.substr(1);
//        }
//        else
//            genomeSize += line.size();
//    }
//    headerFile << "@SQ\tSN:" << header << "\tLN:" << genomeSize << endl;
//    infile.close();
//    headerFile.close();
//    multiFasta = multiFA;
//    return;
//}

//void Assignment::Split_Header(vector<string> &splitted, string &Header,string &mainHeader, long long &index, long long &number, int &V_L, int &row, bool drop_Fragment, bool samFile, bool clear_splitted, bool pacbio, bool saveRow) {
//    std::size_t pos;
//    if(samFile){
//        pos = Header.find("r") + 1;
//    }
//    else{
//        pos = Header.find("@") + 2;
//    }
//    {
//        vector<string>().swap(splitted);
//    }
//    Header=Header.substr(pos);
//    splitted = split(Header,"_");
//    index = atoll(splitted[0].c_str());
//    number = atoll(splitted[1].c_str());
//    V_L = atoi(splitted[2].c_str());
//    if(pacbio){
//        main_Header(Header,mainHeader);
//    }
//    if(drop_Fragment && !pacbio){
//        Header = Header.substr(0,Header.length()-splitted[3].length()-1);
//    }
//    else if(drop_Fragment && pacbio){
//        ostringstream tempHeader;
//        tempHeader << index << "_" << number << "_" << V_L << "_" << mainHeader;
//        Header = tempHeader.str();
//    }
//    if(saveRow){
//        row = atoi(splitted[3].c_str());
//    }
//    if(clear_splitted){
//        vector<string>().swap(splitted);
//    }
//    if(!pacbio){
//        mainHeader = Header;
//    }
//    return;
//}

//void Assignment::main_Header(string &Header, string &mainHeader) {
//    mainHeader = Header.substr(Header.find("@"));
//    return;
//}

//void Assignment::reform_Header(vector<string> &splitted, string &Header, long long &index, long long &number, int new_V_L,  bool clear_splitted, bool pacbio, bool buildIndex, bool newLength, int i) {
//    string mainHeader;
//    ostringstream tempHeader;
//    {
//        vector<string>().swap(splitted);
//    }
//    splitted = split(Header,"_");
//    index = atoll(splitted[0].c_str());
//    number = atoll(splitted[1].c_str());
//    if(pacbio){
//        main_Header(Header, mainHeader);
//    }
//    if(newLength){
//        tempHeader << index << "_" << number << "_" << new_V_L;
//        if(pacbio){
//            tempHeader << "_" << mainHeader;
//        }
//        Header = tempHeader.str();
//    }
//    if(buildIndex){
//        tempHeader << index << "_" << number << "_" << splitted[2] << "_" << indexMaker(i);
//        if(pacbio){
//            tempHeader << "_" << mainHeader;
//        }
//        Header = tempHeader.str();
//    }
//    if(clear_splitted){
//        vector<string>().swap(splitted);
//    }
//    return;
//}

//void Assignment::align_Command(int mode,string indexAddress, string fastQAddress, int V, int depth, int numberOfThreads, string AlignedAddress, int option, int minScore, int N_bowtie, int L_bowtie) {
//    string command = "";
//    if(mode==BOWTIE1){
//        command.append(string(BOWTIE));
//        command.append(" -S ");
//        command.append(indexAddress.c_str());
//        command.append(" ");
//        command.append(fastQAddress.c_str());
//        command.append(" -v ");
//        command.append(NumToStr(V));
//        command.append(" -a -m ");
//        command.append(NumToStr(depth - 1));
//        command.append(" -t --mm -p ");
//        command.append(NumToStr(numberOfThreads));
//        command.append(" -so ");
//        command.append(" > ");
//        command.append(AlignedAddress.c_str());
//        command.append(" 2>> "+outputDir+"log-bowtie.txt");
//        system(command.c_str());
//    }
//    else if(mode==BOWTIE2){
//        command.append("bowtie2");
//        if(option!=0){
//            command.append(" -k ");
//            command.append(NumToStr(option));
//        }
//        else{
//            command.append(" -a ");
//        }
//        command.append(" -t --threads ");
//        command.append(NumToStr(numberOfThreads));
//        command.append(" --local --reorder ");
//        // command.append(NumToStr(minScore));
//        command.append(" -N " + NumToStr(N_bowtie) + " -L " +  NumToStr(L_bowtie) + " -x ");
//        command.append(indexAddress.c_str());
//        command.append(" -U ");
//        command.append(fastQAddress.c_str());
//        command.append(" -S ");
//        command.append(AlignedAddress.c_str());
//        system(command.c_str());
//    }
//    return;
//}

//void Assignment::Split_Token(vector<string> &splitted,string &token, string &Header, string &seq, string &quality) {
//    {
//        vector<string>().swap(splitted);
//    }
//    splitted = split(token,"\t");
//    Header = splitted[0];
//    // This is just for the "@r" sign in most headers
//    seq = splitted[1];
//    quality = splitted[2];
//    {
//        vector<string>().swap(splitted);
//    }
//    return;
//}

//void Assignment::system_sort(string semiSortedFinalReadsAddress, string sortedFinalReadsAddress) {
//    string command = "";
//    command.append("LC_COLLATE=C sort -k 1 ");
//    command.append(semiSortedFinalReadsAddress.c_str());
//    command.append(" > ");
//    command.append(sortedFinalReadsAddress.c_str());
//    system(command.c_str());
//    return;
//}

//void Assignment::system_remove(string fileName){
//    string command="rm " +  fileName;
//    system(command.c_str());
//    return;
//}

//string Assignment::currentPath(){
//    char path[FILENAME_MAX];
//    getcwd(path, sizeof(path));
//    string currentPath(path);
//    return currentPath;
//}

//tuple<bool, long long, long long, long long> Assignment::calculateIntronOrExonAndDistance(long long position) {
//    bool isInExon = false;
//    long long index = 0;
//    long long fromRight = 0;
//    long long fromLeft = 0;

//    long long i;
//    for (i = 1; i < this->iDepth->size(); i++) {
//        if (this->iDepth->at(i).start <= position && position <= this->iDepth->at(i).end) {
//            isInExon = true;
//            index = i;
//            fromLeft = abs(position - this->iDepth->at(i).start);
//            fromRight = abs(position - this->iDepth->at(i).end);
//            return make_tuple(isInExon, index, fromLeft, fromRight);
//        } else if (position < this->iDepth->at(i).start && this->iDepth->at(i-1).end < position) {
//            isInExon = false;
//            index = i;
//            fromLeft = abs(position - this->iDepth->at(i-1).end);
//            fromRight = abs(position - this->iDepth->at(i).start);
//            return make_tuple(isInExon, index, fromLeft, fromRight);
//        }
//    }

//    if (position < this->iDepth->at(0).start) {
//        isInExon = false;
//        index = 0;
//        fromLeft = 0;
//        fromRight = abs(position - this->iDepth->at(0).start);
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    } else if (this->iDepth->at(0).start <= position && position <= this->iDepth->at(0).end) {
//        isInExon = true;
//        index = 0;
//        fromLeft = abs(position - this->iDepth->at(0).start);;
//        fromRight = abs(position - this->iDepth->at(0).end);
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    } else if (this->iDepth->at(this->iDepth->size()-1).end < position) {
//        isInExon = false;
//        index = this->iDepth->size();
//        fromLeft = abs(position - this->iDepth->at(this->iDepth->size()-1).end);
//        fromRight = 0;
//        return make_tuple(isInExon, index, fromLeft, fromRight);
//    }
//    return make_tuple(isInExon, index, fromLeft, fromRight);
//}

//void Assignment::addIntervalToDepth(long long a, long long b, long long insert_index) {

//    this->iDepth->insert(this->iDepth->begin()+insert_index, interval());
//    this->iDepth->at(insert_index).start = a;
//    this->iDepth->at(insert_index).end = b;
//    return;

//}

//void Assignment::addIntervalToDepth(long long a, long long b) {

//    if (this->iDepth->size() == 0) {
//        this->iDepth->push_back(interval());
//        this->iDepth->at(0).start = a;
//        this->iDepth->at(0).end = b;
//        return;
//    }
//    for (long long i = 0; i < this->iDepth->size(); i++) {
//        if (this->iDepth->at(i).start <= a && this->iDepth->at(i).end >= b) {
//            return;
//        }
//        if (this->iDepth->at(i).start > b || this->iDepth->at(i).end < a) {
//            continue;
//        } else {
//            long long tempA = this->iDepth->at(i).start;
//            long long tempB = this->iDepth->at(i).end;
//            this->iDepth->erase(this->iDepth->begin()+i);
//            this->addIntervalToDepth(min(tempA, a), max(tempB, b));
//            return;
//        }
//    }
//    long long insert_index = 0;
//    while (insert_index < this->iDepth->size() && this->iDepth->at(insert_index).start < b) {
//        insert_index++;
//    }
//    this->iDepth->insert(this->iDepth->begin()+insert_index, interval());
//    this->iDepth->at(insert_index).start = a;
//    this->iDepth->at(insert_index).end = b;
//    return;

//}

//void Assignment::updateDepth(long long firstPosition, long long lastPosition, bool isForward, int partLen) {
//    bool isInExon;
//    long long indexNumber;
//    long long fromRight;
//    long long fromLeft;

//    tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(firstPosition);

//    if (isInExon) {
//        if (isForward) {
//            if (lastPosition + partLen > this->iDepth->at(indexNumber).end) {
//                this->iDepth->at(indexNumber).end = lastPosition + partLen;
//                this->concatExons(indexNumber);
//            }
//        } else {
//            if (this->iDepth->at(indexNumber).start > lastPosition) {
//                this->iDepth->at(indexNumber).start = lastPosition;
//                this->concatExons(indexNumber);
//            }
//        }
//    } else {
//        if (isForward) {
//            tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(lastPosition+partLen);
//            if (isInExon) {
//                this->iDepth->at(indexNumber).start = firstPosition;
//                this->concatExons(indexNumber);
//            } else {
//                this->addIntervalToDepth(firstPosition, lastPosition+partLen, indexNumber);
//            }
//        } else {
//            tie(isInExon,indexNumber,fromLeft,fromRight) = this->calculateIntronOrExonAndDistance(lastPosition);
//            if (isInExon) {
//                this->iDepth->at(indexNumber).end = firstPosition + partLen;
//                this->concatExons(indexNumber);
//            } else {
//                this->addIntervalToDepth(lastPosition, firstPosition+partLen, indexNumber);
//            }
//        }
//    }
//}

//void Assignment::concatExons(long long index) {

//    vector<interval>::iterator current = this->iDepth->begin()+index;
//    vector<interval>::iterator next = current + 1;
//    vector<interval>::iterator previous = current - 1;

//    if (current->end >= next->start) {
//        long long end = max(current->end, next->end);
//        current->end = end;
//        this->iDepth->erase(next);
//    }
//    if (current->start <= previous->end) {
//        long long start = min(current->start, previous->start);
//        current->start = start;
//        this->iDepth->erase(previous);
//    }

//    return;
//}

//string Assignment::convertNumToStr(long long number) {
//    ostringstream ostr;
//    ostr<<number;
//    return ostr.str();
//}

//void Assignment::runForOnGenes() {
//    for (long long i = 0; i < genes->size(); i++) {
//        bool goToNextGene = false;
//        for (int j = 0; j < genes->at(i).inExon.size() && !goToNextGene; j++) {
//            for (long long k = 0; k < iDepth->size() && !goToNextGene; k++) {
//                long long start = iDepth->at(k).start;
//                long long end = iDepth->at(k).end;

//                long long compare = genes->at(i).inExon.at(j).start;
//                if (start < compare && compare < end) {
//                    goToNextGene = true;
//                    for (int p = 0; p < genes->at(i).inExon.size(); p++) {
//                        addIntervalToDepth(genes->at(i).inExon.at(p).start, genes->at(i).inExon.at(p).end);
//                    }
//                    for (int p = 0; p < genes->at(i).inRead.size(); p++) {
//                        long long index = genes->at(i).inRead.at(p).index-1;
//                        long size = genes->at(i).inRead.at(p).fragments.size();
//                        this->reads[index].firstPosition = genes->at(i).inRead.at(p).fragments.at(0).firstPosition;
//                        this->reads[index].lastPosition = genes->at(i).inRead.at(p).fragments.at(size-1).lastPosition;
//                        this->reads[index].firstFragment = genes->at(i).inRead.at(p).fragments.at(0).firstFragment;
//                        this->reads[index].lastFragment = genes->at(i).inRead.at(p).fragments.at(size-1).lastFragment;
//                        this->reads[index].cigar = genes->at(i).inRead.at(p).cigar;
//                    }

//                } else {
//                    compare = genes->at(i).inExon.at(j).end;
//                    if (start < compare && compare < end) {
//                        goToNextGene = true;
//                        for (int p = 0; p < genes->at(i).inExon.size(); p++) {
//                            addIntervalToDepth(genes->at(i).inExon.at(p).start, genes->at(i).inExon.at(p).end);
//                        }
//                        for (int p = 0; p < genes->at(i).inRead.size(); p++) {
//                            long long index = genes->at(i).inRead.at(p).index-1;
//                            long size = genes->at(i).inRead.at(p).fragments.size();
//                            this->reads[index].firstPosition = genes->at(i).inRead.at(p).fragments.at(0).firstPosition;
//                            this->reads[index].lastPosition = genes->at(i).inRead.at(p).fragments.at(size-1).lastPosition;
//                            this->reads[index].firstFragment = genes->at(i).inRead.at(p).fragments.at(0).firstFragment;
//                            this->reads[index].lastFragment = genes->at(i).inRead.at(p).fragments.at(size-1).lastFragment;
//                            this->reads[index].cigar = genes->at(i).inRead.at(p).cigar;

//                        }
//                    }
//                }
//            }
//        }
//    }
//    return;
//}

//void Assignment::createNewRead(int d, string outputNeme) {
//    string name = outputDir+this->firstOutputName+"_reads_d="+this->NumToStr(this->shift)+".fq";
//    if (d == 0) {
//        name = outputDir+this->firstOutputName+"_reads_d=0.fq";
//    }
//    ifstream ifstr(name.c_str());
//    name = outputDir+outputNeme+"_reads_d="+this->convertNumToStr(d)+".fq";
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
//                line = line.substr(d);
//                ofstr<<line<<"\n";
//            }
//        } else if (lineCounter % 4 == 3) {
//            if (this->reads[stoll(segment.c_str())-1].flag == -1) {
//                line = line.substr(d);
//                ofstr<<"+\n"<<line<<"\n";
//            }
//        }
//        ++lineCounter;
//    }
//    ifstr.close();
//    ofstr.close();
//    return;
//}

//void Assignment::AddToJenes(gene &tempGene){
//    for (int geneIterator = 0; geneIterator < genes->size(); geneIterator++) {
//        for (int trueExonIterator = 0; trueExonIterator < genes->at(geneIterator).inExon.size(); trueExonIterator++) {
//            for (int tempExonIteraor = 0; tempExonIteraor < tempGene.inExon.size(); tempExonIteraor++) {
//                if ((genes->at(geneIterator).isForward == tempGene.isForward) && ((tempGene.start <= genes->at(geneIterator).end && tempGene.start >= genes->at(geneIterator).start) || (tempGene.end <= genes->at(geneIterator).end  && tempGene.end >= genes->at(geneIterator).start))) {
//                    if ((tempGene.inExon[tempExonIteraor].start <= genes->at(geneIterator).inExon[trueExonIterator].end &&tempGene.inExon[tempExonIteraor].start >= genes->at(geneIterator).inExon[trueExonIterator].start) || (tempGene.inExon[tempExonIteraor].end <= genes->at(geneIterator).inExon[trueExonIterator].end &&tempGene.inExon[tempExonIteraor].end >= genes->at(geneIterator).inExon[trueExonIterator].start) ) {
//                        genes->at(geneIterator).start = min(genes->at(geneIterator).start, tempGene.start);
//                        genes->at(geneIterator).end = max(genes->at(geneIterator).end, tempGene.end);
//                        for (int readIterator = 0; readIterator < tempGene.inRead.size(); readIterator++) {
//                            genes->at(geneIterator).inRead.push_back(tempGene.inRead[readIterator]);
//                        }
//                        //join to exonVector
//                        for (int mergeIterator = 0; mergeIterator < tempGene.inExon.size(); mergeIterator++) {
//                            this->mergeExonsInGenesAndTempGene(tempGene.inExon[mergeIterator].start, tempGene.inExon[mergeIterator].end, geneIterator);
//                        }
//                        return;
//                    }
//                }
//                else if (tempGene.end <= genes->at(geneIterator).start){
//                    //add temgene to this position of vector
//                    genes->insert(genes->begin()+geneIterator, tempGene);
//                    return;
//                }
//            }
//        }
//    }
//    this->genes->push_back(tempGene);
//    return;
//}

//void Assignment::mergeExonsInGenesAndTempGene(long long a, long long b, long long index) {

//    for (int i = 0; i < genes->at(index).inExon.size(); i++) {
//        if (genes->at(index).inExon.at(i).start <= a && genes->at(index).inExon.at(i).end >= b) {
//            return;
//        }
//        if (genes->at(index).inExon.at(i).start > b || genes->at(index).inExon.at(i).end < a) {
//            continue;
//        } else {
//            long long tempA = genes->at(index).inExon.at(i).start;
//            long long tempB = genes->at(index).inExon.at(i).end;
//            genes->at(index).inExon.erase(genes->at(index).inExon.begin()+i);
//            this->mergeExonsInGenesAndTempGene(min(tempA, a), max(tempB, b), index);
//            return;
//        }
//    }
//    long long insert_index = 0;
//    while (insert_index < genes->at(index).inExon.size() && genes->at(index).inExon.at(insert_index).start < b) {
//        insert_index++;
//    }
//    genes->at(index).inExon.insert(genes->at(index).inExon.begin()+insert_index, interval());
//    genes->at(index).inExon.at(insert_index).start = a;
//    genes->at(index).inExon.at(insert_index).end = b;
//    return;
//}
