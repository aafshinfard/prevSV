#include "Signalling.h"
#include <string>
long long informCounts = 0;
Signalling::Signalling(iRead* reads , string readName, string genomeName, int chunkSize, int numberOfThread, int numGap, string outputDir,long long numberOfReads1)
{
    this->reads = reads;
    this->readName = readName;
    this->genomeName = genomeName;
    this->chunkSize = chunkSize;
    this->numberOfThread = numberOfThread;
    this->numGap = numGap;
    this->outputDir = outputDir;
    this->numberOfReads1 = numberOfReads1;
    this->mainSignal = new int[this->numberOfReads1];
    this->genomeLength = getGenomeLength();
    eventBinGenome.resize(this->genomeLength);
    //bGLeft.resize(this->genomeLength);

}
void Signalling::automation(){
    cout<<"\n=============================================\n";
    cout<<  "============= Part 2 : Calling ==============";
    cout<<"\n=============================================\n";
    cout<<  "=============================================\n";

    cout<<  " building Events Binary Genome...";
    buildEventsBinGenome();
    cout<<"\n ...........................Done!\n";

    cout<<  " building Events base on BinGenome...";
    buildEvents();
    cout<<"\n ...............................Done!\n";

    cout<<  " finding Paired Events and conections ...";
    buildPairedEvents();
    cout<<"\n ...................................Done!\n";

    cout<<  " printing  Events Intervals from BinGenome ...";
    printEventIntervals();
    cout<<"\n ........................................Done!\n";

    cout<<  " printing Paired Events and conections ...";
    printPaiedEvents();
    cout<<"\n ....................................Done!\n";

}

void inline Signalling::changeBinGenome(long long start, long long end){

    informCounts++;
    if(start > end){
        long long temp = start;
        start = end;
        end = temp;
    }
    if(start == 33442)
        cout<<"test";
    //bGRChangeMutex.lock();
    for(long long i = start; i < end+1 ; i++)
        eventBinGenome[i] = true;
    //bGRChangeMutex.unlock();

}

void Signalling::buildEventsBinGenome(){
    for (long long i = 0; i < this->numberOfReads1; i++) {
        int lastFragNum = ((reads[i].length-1)/chunkSize)+1 ;
        if( reads[i].firstFragment == 0 && reads[i].lastFragment == 0 ){
            // unMapped reads:
            unMappedReads.push_back(reads[i]);
        }else{
            // Mapped reads:
            if( ! (reads[i].firstFragment == 1 && reads[i].lastFragment == lastFragNum ) ){
                // Informative Reads and ...
                informativeReads.push_back(reads[i]);
                // Informative l-mers:
                long long a;
                if( i == 584 ){
                    a = informativeReads.size();
                    cout<<a;
                }
                if( reads[i].firstFragment != 1  )
                    changeBinGenome(reads[i].firstPosition, reads[i].firstPosition+chunkSize-1);
                if( reads[i].lastFragment != lastFragNum )
                    changeBinGenome(reads[i].lastPosition, reads[i].lastPosition+chunkSize-1);
                for(int j=1 ; j<reads[i].firstFragClus2.size() ; j++ ){
                    if( reads[i].firstFragClus2[j] != 1 )
                        changeBinGenome(reads[i].firstPosClus2[j], reads[i].firstPosClus2[j]+chunkSize-1);
                    if( reads[i].lastFragClus2[j] != lastFragNum )
                        changeBinGenome(reads[i].lastPosClus2[j], reads[i].lastPosClus2[j]+chunkSize-1);

                }

            }
        }
    }
    cout<<"\n==> There is "<<informativeReads.size()<<" informative and "<<unMappedReads.size()<<" unmapped Reads.";
}
inline long long Signalling::correspondingEvent(long long location){
    for(int i = 0 ; i < eventIntervals.size() ; i++ ){
        if( location >=  eventIntervals[i].start && location <=  eventIntervals[i].end  )
            return i;
    }
    return -1;
}

void Signalling::buildEvents(){
    // Step1: finding Intervals of Events using Bingenome
    EventInterval tempEventInterval;
    bool flag = false;
    for( long long i = 0 ; i < eventBinGenome.size() ; i++){
        if(flag == false && eventBinGenome[i] == true ){
            flag = true;
            tempEventInterval.start = i;
        }
        if(flag == true && eventBinGenome[i] == false){
            flag = false;
            tempEventInterval.end = i-1;
            eventIntervals.push_back(tempEventInterval);
        }
    }
    cout<<"\n==> There is "<<eventIntervals.size()<<" Intervals of  true on BinGenome.";
    //eventConnectedEvents.resize( eventIntervals.size() );
    //eventInformativeLMers.resize( eventIntervals.size() );

    //vector<int> tempReadCE;
    //vector<int> connectedEvents;
    //vector<fragNum_Sign> readFragNumSign;

    //Cluster tempCluster;

    // Step2: Building Events Information:
    eventInformativeClusters.resize( eventIntervals.size() );
    eventPairedEvents.resize( eventIntervals.size() );
    eventSingleCountsMinus.resize( eventIntervals.size() , 0 );
    eventSingleCountsPlus.resize( eventIntervals.size() , 0 );

    vector<cCluster> tempFirst;
    vector<cCluster> tempLast;
    long long corespEvent;

    for (long long i = 0; i < informativeReads.size() ; i++) {
        if(i==52)
            cout<<"test";

        // informativeReads[i] clusters to Events :
        //cout<<"\n Read "<<informativeReads[i].index<<" with Length "<<informativeReads[i].length<<endl;
        int lastFragNum = ((informativeReads[i].length-1)/chunkSize)+1 ;
        //cout<<lastFragNum<<endl;
        tempFirst.clear();
        tempLast.clear();
        //bool tempoflag = true;
        iRead tempRead = informativeReads[i];
        if( informativeReads[i].firstFragment != 1 ){
            //tempoflag = false;
            corespEvent = correspondingEvent(informativeReads[i].firstPosition);
            cout<<"\n R "<<informativeReads[i].index<<" F-Informative in Event "<<corespEvent;
            if(corespEvent != -1 ){
                tempFirst.push_back(cCluster(Cluster(informativeReads[i].index
                                                     ,informativeReads[i].firstFragment,informativeReads[i].lastFragment
                                                     ,informativeReads[i].firstPosition,informativeReads[i].lastPosition
                                                     , informativeReads[i].length,informativeReads[i].flag
                                                     ,(informativeReads[i].flag == 0 ? +1:-1)
                                                     )
                                             ,informativeReads[i].firstFragment
                                             ,corespEvent
                                             )
                                    );
//                tempCluster.index = informativeReads[i].index;
//                tempCluster.firstFragment = informativeReads[i].firstFragment;
//                tempCluster.firstPosition = informativeReads[i].firstPosition;
//                tempCluster.lastFragment = informativeReads[i].lastFragment;
//                tempCluster.lastPosition = informativeReads[i].lastPosition;
//                tempCluster.length = informativeReads[i].length;
//                tempCluster.flag = informativeReads[i].flag;
//                tempCluster.sign = (tempCluster.flag == 1 ? +1:-1);
                //eventInformativeLMers[corespEvent].push_back(tempCluster);
            }else cout<<"\nErr:no Coressponding Event!";
        }
        if( informativeReads[i].lastFragment != lastFragNum ){
            //tempoflag = false;
            corespEvent = correspondingEvent(informativeReads[i].lastPosition  );
            cout<<"\n R "<<informativeReads[i].index<<" L-Informative in Event "<<corespEvent;
            if(corespEvent != -1 ){
                tempLast.push_back(cCluster(Cluster(informativeReads[i].index
                                                     ,informativeReads[i].firstFragment,informativeReads[i].lastFragment
                                                     ,informativeReads[i].firstPosition,informativeReads[i].lastPosition
                                                     , informativeReads[i].length,informativeReads[i].flag
                                                     ,(informativeReads[i].flag == 0 ? -1:+1)
                                                     )
                                             ,informativeReads[i].lastFragment,corespEvent
                                             )
                                    );
            }else cout<<"\nErr:no Coressponding Event!";
        }

        //bool tempoflag2 = false;
//        if(informativeReads[i].index == 6)
//            cout<<"6";
        for(int j=1 ; j < informativeReads[i].firstFragClus2.size() ; j++ ){

            if( informativeReads[i].firstFragClus2[j] != 1  ){
               // tempoflag2 = false;
                corespEvent = correspondingEvent( informativeReads[i].firstPosClus2[j] );
                cout<<"\n\t\t and "<<" F-Informative in Event "<<corespEvent;
                cout<<"\n\t\t sign is : "<<(informativeReads[i].flagClus2[j] == 0 ? +1:-1)<<" flag: "<<informativeReads[i].flagClus2[j];
                if(corespEvent != -1 ){
                    tempFirst.push_back(cCluster(Cluster(informativeReads[i].index
                                                         ,informativeReads[i].firstFragClus2[j],informativeReads[i].lastFragClus2[j]
                                                         ,informativeReads[i].firstPosClus2[j],informativeReads[i].lastPosClus2[j]
                                                         , informativeReads[i].length,informativeReads[i].flagClus2[j]
                                                         ,(informativeReads[i].flagClus2[j] == 0 ? +1:-1)
                                                         )
                                                 ,informativeReads[i].firstFragClus2[j],corespEvent
                                                 )
                                        );
                }else cout<<"\nErr:no Coressponding Event!";
            }
            if( informativeReads[i].lastFragClus2[j] != lastFragNum ){
                //tempoflag2 = false;
                corespEvent = correspondingEvent( informativeReads[i].lastPosClus2[j] );
                cout<<"\n\t\t and "<<" L-Informative in Event "<<corespEvent;
                cout<<"\n\t\t sign is : "<<(informativeReads[i].flagClus2[j] == 0 ? -1:+1)<<" flag: "<<informativeReads[i].flagClus2[j];
                if(corespEvent != -1 ){
                    tempLast.push_back(cCluster(Cluster(informativeReads[i].index
                                                         ,informativeReads[i].firstFragClus2[j],informativeReads[i].lastFragClus2[j]
                                                         ,informativeReads[i].firstPosClus2[j],informativeReads[i].lastPosClus2[j]
                                                         , informativeReads[i].length,informativeReads[i].flagClus2[j]
                                                         ,(informativeReads[i].flagClus2[j] == 0 ? -1:+1)
                                                         )
                                                 ,informativeReads[i].lastFragClus2[j],corespEvent
                                                 )
                                        );
                }else cout<<"\nErr:no Coressponding Event!";
            }
        }

        // Step3: adding Connectivity information using tempFirst and tempLast
        long long connectedWith = -1;
        for(int j = 0 ; j < tempFirst.size() ; j++){
            cout<<"\n----FirstFragment "<<tempFirst[j].fragNum<<" from Event "<<tempFirst[j].corespEvent;
        }
        for(int j = 0 ; j < tempLast.size() ; j++){
            cout<<"\n----LastFragment "<<tempLast[j].fragNum<<" from Event "<<tempLast[j].corespEvent;
        }
        for(int j = 0 ; j < tempFirst.size() ; j++){
            int maxLast = -1;
            for(int k = 0 ; k < tempLast.size() ; k++)
                if(tempLast[k].fragNum < tempFirst[j].fragNum
                        &&
                        tempLast[k].fragNum > (maxLast > -1 ? tempLast[maxLast].fragNum : -1)
                   )
                    maxLast = k;

            if(maxLast>-1)
                cout<<"\n\t\t -->firtFrag"<<tempFirst[j].fragNum<<" and lastFrag "<<tempLast[maxLast].fragNum;
            else
                cout<<"\n\t\t -->firtFrag"<<tempFirst[j].fragNum<<" having no match ";
//            cout<<"\n-----sign is :"<<tempFirst[j].cluster.sign;
            // tempFirst[j] is connected with tempLast[maxLast] so :
            if(maxLast > -1){
                connectedWith = tempLast[maxLast].corespEvent;
                eventInformativeClusters[tempFirst[j].corespEvent].push_back(
                            ClusterAndConnection(
                                tempFirst[j].cluster,               // this cluster (in event tempFirst[j].corespEvent) is connected with
                                tempLast[maxLast].corespEvent,      // this Event of sign
                                tempLast[maxLast].cluster.sign,
                                abs(tempLast[maxLast].fragNum - tempFirst[j].fragNum)
                                )
                            ); // this
                eventInformativeClusters[tempLast[maxLast].corespEvent].push_back(
                            ClusterAndConnection(
                                // cluster (tempLast[maxLast]) connected
                                // with Event ( tempFirst[j].corespEvent with sign tempFirst[j].cluster.sign)
                                tempLast[maxLast].cluster,     // this cluster (in event tempLast[maxLast].corespEvent) is connected with
                                tempFirst[j].corespEvent,      // this Event of sign
                                tempFirst[j].cluster.sign,
                                abs(tempLast[maxLast].fragNum - tempFirst[j].fragNum)
                                )
                            ); // this
                tempLast.erase( tempLast.begin()+maxLast );

            }else eventInformativeClusters[tempFirst[j].corespEvent].push_back(
                        ClusterAndConnection(
                                             tempFirst[j].cluster,  // this cluster (in event tempLast[maxLast].corespEvent) is connected with
                                             -1,                    // nobody
                                             0,
                                             -1
                                             )
                        );
        }
        for(int j = 0 ; j < tempLast.size() ; j++){
            if( j == 1){
                cout<<"\n ============ Why ??? ============\n";
                break;
            }
            eventInformativeClusters[tempLast[j].corespEvent].push_back(
                        ClusterAndConnection(
                            tempLast[j].cluster,     // this cluster (in event tempLast[maxLast].corespEvent) is connected with
                            -1,                      // nobody
                            0,                       //
                            -1
                            )
                        ); // this
        }
        if(false){
        // ///////////////////////////////////////////////
        // connectedEvents : events in wich this read shares a cluster !!!
        // there can be another data named: vector<vector<vector<int>>> readConnectedEvents;
        // add this data to eventConnectedEvents !!
        // ///////////////////////////////////////////////
        // Step3 - mode 1 :
//        for(int j = 0 ; j < connectedEvents.size() ; j++)
//            for(int k = 0 ; k < connectedEvents.size() ; k++)
//                if( connectedEvents[j] != connectedEvents[k] )
//                    eventConnectedEvents[connectedEvents[j]].push_back(connectedEvents[k]);
        // Step3 - mode 2 :
//        tempReadCE.clear();
//        for(int j = 0 ; j < connectedEvents.size() ; j++){
//            for(int k = 0 ; k < connectedEvents.size() ; k++)
//                if( connectedEvents[j] != connectedEvents[k] )
//                    tempReadCE.push_back(connectedEvents[k]);
//            readConnectedEvents[connectedEvents[j]].push_back(tempReadCE);
//        }
        // ///////////////////////////////////////////////
        }
    }
    cout<<"\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
//    for(int i = 0 ; i < eventInformativeClusters.size() ; i++ ){
//        cout<<"\n =============\nEvent "<<i;
//        for(int j=0 ; j<eventInformativeClusters[i].size(); j++ ){
//            cout<<"\n -- of sign "<<eventInformativeClusters[i][j].cluster.sign <<" connected with "
//               <<eventInformativeClusters[i][j].connectedEvent<<" of sign "<<eventInformativeClusters[i][j].connectedSign;
//        }
//    }
}

void Signalling::buildPairedEvents(){
    //eventInformativeClusters
    int lastElement;
    vector<long long> tempPairs;
    for( long long i = 0 ; i < eventInformativeClusters.size() ; i++ ){
        vector<ClusterAndConnection>  temp;
        if(i==5)
            temp = eventInformativeClusters[i];
        tempPairs.clear();
        cout<<"\n Event "<<i;
        for( int j = 0 ; j < eventInformativeClusters[i].size() ; j++ ){
            //eventSingleCounts
            cout<<"\nof sign"<<eventInformativeClusters[i][j].cluster.sign<<" with event "<<eventInformativeClusters[i][j].connectedEvent ;
            if( eventInformativeClusters[i][j].connectedEvent == -1 ){
                //eventSingleCountsPlus
                if(eventInformativeClusters[i][j].cluster.sign == 0 ){
                    eventSingleCountsMinus[i]++;
                    eventSingleCountsPlus[i]++;
                }
                else if(eventInformativeClusters[i][j].cluster.sign == -1 )
                    eventSingleCountsMinus[i]++;
                else if(eventInformativeClusters[i][j].cluster.sign == +1)
                    eventSingleCountsPlus[i]++;
            }

            if( eventInformativeClusters[i][j].connectedEvent >= i ){
                if( std::find( tempPairs.begin() , tempPairs.end(), eventInformativeClusters[i][j].connectedEvent ) == tempPairs.end() ){
                    eventPairedEvents[i].push_back(PairedEvents(eventInformativeClusters[i][j].connectedEvent));
                    //pairedEvents.push_back( pairedEventGraph(i, eventInformativeClusters[i][j].connectedEvent) );
                    tempPairs.push_back( eventInformativeClusters[i][j].connectedEvent );
                    short e1Sign = eventInformativeClusters[i][j].cluster.sign,
                          e2Sign = eventInformativeClusters[i][j].connectedSign;
                    lastElement = eventPairedEvents[i].size()-1;
                    if( e1Sign == -1 && e2Sign == -1 )
                        eventPairedEvents[i][ lastElement ].countMinusMinus++;
                    if( e1Sign == -1 && e2Sign == +1 )
                        eventPairedEvents[i][ lastElement ].countMinusPlus++;
                    if( e1Sign == +1 && e2Sign == -1 )
                        eventPairedEvents[i][ lastElement ].countPlusMinus++;
                    if( e1Sign == +1 && e2Sign == +1 )
                        eventPairedEvents[i][ lastElement ].countPlusPlus++;
                    for( int k = j+1 ; k < eventInformativeClusters[i].size() ; k++ ){
                        if(eventInformativeClusters[i][j].connectedEvent == eventInformativeClusters[i][k].connectedEvent){
                            short e1Sign = eventInformativeClusters[i][k].cluster.sign,
                                  e2Sign = eventInformativeClusters[i][k].connectedSign;
                            if( e1Sign == -1 && e2Sign == -1 )
                                eventPairedEvents[i][ lastElement ].countMinusMinus++;
                            if( e1Sign == -1 && e2Sign == +1 )
                                eventPairedEvents[i][ lastElement ].countMinusPlus++;
                            if( e1Sign == +1 && e2Sign == -1 )
                                eventPairedEvents[i][ lastElement ].countPlusMinus++;
                            if( e1Sign == +1 && e2Sign == +1 )
                                eventPairedEvents[i][ lastElement ].countPlusPlus++;
                        }
                    }
                }else{
                    //std::find( pairedEvents.begin() , pairedEvents.end(), eventInformativeClusters[i][j].connectedEvent )
                }

            }
        }
    }
    cout<<"\n=========================================\n";
    for(int i = 0 ; i < eventPairedEvents.size() ; i++ ){
        cout<<"\n --- Event "<<i<<" pairs: ";
        for(int j = 0 ; j < eventPairedEvents[i].size() ; j++ ){
            cout<<"\n\t\t=> "<<eventPairedEvents[i][j].destination;
            cout<<"\n\t\t  "<<eventPairedEvents[i][j].countMinusMinus<<" - "<<eventPairedEvents[i][j].countMinusPlus<<" - "<<eventPairedEvents[i][j].countPlusMinus<<" - "<<eventPairedEvents[i][j].countPlusPlus;
        }
    }

}

void Signalling::printEventIntervals(){

    string name = outputDir+"Z_EventIntervals.txt";

    ofstream ofstr(name.c_str());
    ofstr << " Event\tStart\tEnd\tLength\tClusterCount"<<endl;
    for( int i = 0 ; i < eventIntervals.size() ; i++ ){
        ofstr<<i<<"\t"<<eventIntervals[i].start<<"\t"<<eventIntervals[i].end<<"\t"<<eventIntervals[i].end - eventIntervals[i].start<<"\t"<< eventInformativeClusters[i].size() <<endl;
        for( int j = 0 ; j < eventInformativeClusters[i].size() ; j++ ){

        }
    }
}

void Signalling::printPaiedEvents(){

    string name = outputDir+"Z_PairedEvents.txt";

    ofstream ofstr(name.c_str());
    ofstr << " Event1\t[From\tto]\t\tEvent2\t[- -]\t[- +]\t[+ -]\t[+ +]\t[- 0]\t[+ 0]"<<endl;
    //cout << " Event1\tEvent2\t[- -]\t[- +]\t[+ -]\t[+ +]\t[- 0]\t[+ 0]"<<endl;
    string delim = " - ";
    for (long long i = 0; i < eventPairedEvents.size() ; i++) {
        //ofstr <<"\n";
        //cout <<"\n";
        if(eventPairedEvents[i].size() ==0 ){
            ofstr <<" "<< i <<"\t"<<eventIntervals[i].start<<"\t"<<eventIntervals[i].end<<"\t\t"<< delim <<"\t"
                  <<delim<<"\t"<<delim<<"\t"
                  <<delim<<"\t"<<delim<<"\t"
                  << (eventSingleCountsMinus[i] == 0 ? delim : to_string(eventSingleCountsMinus[i]))<<"\t"
                  << (eventSingleCountsPlus[i]  == 0 ? delim : to_string(eventSingleCountsPlus[i])) <<endl;
        }
        for (long long j = 0; j < eventPairedEvents[i].size() ; j++) {
            ofstr <<" "<< i <<"\t"<<eventIntervals[i].start<<"\t"<<eventIntervals[i].end<<"\t\t"
                  << eventPairedEvents[i][j].destination <<"\t"
                  << (eventPairedEvents[i][j].countMinusMinus == 0 ? delim :to_string(eventPairedEvents[i][j].countMinusMinus) )  <<"\t"
                  << (eventPairedEvents[i][j].countMinusPlus  == 0 ? delim :to_string(eventPairedEvents[i][j].countMinusPlus) )<<"\t"
                  << (eventPairedEvents[i][j].countPlusMinus  == 0 ? delim :to_string(eventPairedEvents[i][j].countPlusMinus) )<<"\t"
                  << (eventPairedEvents[i][j].countPlusPlus   == 0 ? delim :to_string(eventPairedEvents[i][j].countPlusPlus) )<<"\t"
                  << (eventSingleCountsMinus[i] == 0 ? delim : to_string(eventSingleCountsMinus[i]))<<"\t"
                  << (eventSingleCountsPlus[i]  == 0 ? delim : to_string(eventSingleCountsPlus[i])) <<endl;
        }

    }
    ofstr.close();
    cout<<"\n==============\n InformCounts: "<<informCounts<<endl;
}

void Signalling::typeMatching(){
    for( long long i = 0 ; i < eventPairedEvents.size() ; i++){
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        // //================== WE ARE HERE ===================//
        ;
    }
}

void Signalling::buildSignal(){



    cout<<"HEHE";
//    for (long long i = 0; i < this->numberOfReads1; i++) {
//        //-/ Ameerosein edits
//        //ofstr << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" << this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition << "\t" << this->reads[i].flag << "\t" << this->reads[i].cigar << "\t";
//        ofstr << this->reads[i].index << "\t" << this->reads[i].length << "\t" << this->reads[i].d << "\t" <<
//                 this->reads[i].firstFragment << "\t" << this->reads[i].lastFragment << "\t" << this->reads[i].firstPosition << "\t" << this->reads[i].lastPosition <<"\t"<<reads[i].flag<<
//                 "\t"<< (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].firstFragClus2.at(1):0 ) << "\t" << (this->reads[i].firstFragClus2.size()>1 ? this->reads[i].lastFragClus2.at(1):0)<<
//                 "\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].firstPosClus2.at(1):0)<<"\t"<<(this->reads[i].firstFragClus2.size()>1?this->reads[i].lastPosClus2.at(1):0)<<
//                 "\t"<< (this->reads[i].firstFragClus2.size()>1?this->reads[i].flagClus2.at(1):0); //<< "\t" << this-   >reads[i].cigar << "\t";
//    }
}

long long Signalling::getGenomeLength(){

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

Signalling::~Signalling() {
    delete [] this->reads;
}
