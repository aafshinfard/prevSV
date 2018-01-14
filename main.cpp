//
//  main.cpp
//  SV-Detection
//
//  Created by Ameerhosein Afshinfard - BRLab - Sharif Uiversity of echnology
//  Copyright (c) 2016 Ameerhosein Afshinfard. All rights reserved.
//

#include "SVAnchor.h"
void help();
int main(int argc, const char * argv[]) {
    /*
     There is some assumption for running the code:
        1. The read length must be at least 3 times greater than chunk size (the other for shift)
        2. There must be no _ in header of reads
     */

    // ================================|
    // ================================|
    //> Settings:                      |
    // ================================|
    string readName = "readst";    //string readName = "reads_k12.fq";
    string genomeName = "Ref_SV.fa"; //string genomeName = "E_coli_K12_DH10B.fa";    //string genomeName = "chr19.fa";
    string indexName = "Ref_SV"; //string indexName = "E_coli_K12_DH10B";        //string indexName = "chr19";
    // ====
    long long numberOfReads = 20000/4;
    int readLength = 148;
    int chunkSize = 25;
    // ====
    int numberOfThread = 8;
    // ====
    int dashV = 1;
    int numGap = 5;//numGap in Anchoring | error in Assignment
    int anchoringShift = 5;
    int assignmentShift = 10; //d = chunkSize - assignmentShift
    bool isSlideNotShift = true;
    string outputDir = "/home/ameer/SVAnchoring/SV_out2/";
    //string outputDir = "/home/ubuntu/SVAnchoring/SV_out2/";
    bool runAnchoring = true;

    // ================================|
    //> End of Settings:               |
    // ================================|


    cout<<"\n=========|================================|=========\n";
    cout<<  "_________|           WELCOME TO           |_________\n";
    cout<<  "_________|     SV-Detection Algorithm     |_________\n";
    cout<<  "_________|      by BRLab @ SharifUni.     |_________\n";
    cout<<  "=========|================================|=========\n";
    SVAnchor *mySVAnchor = new SVAnchor(readName, genomeName, indexName, chunkSize, numberOfThread, dashV, numGap, anchoringShift, assignmentShift, isSlideNotShift, readLength, outputDir, numberOfReads);
    mySVAnchor-> automation(runAnchoring);
    delete mySVAnchor;


    for (int i = 1; i < argc; ) {
        if (string(argv[i]) == "-o") {
            outputDir = atoi(argv[i + 1]);
            i += 2;
            if (i==argc) break;
        }
        if (string(argv[i]) == "-v") {
            dashV = atoi(argv[i + 1]);
            i += 2;
            if (i==argc) break;
        }
        if (string(argv[i]) == "-h") {
            help();
            break;
        }
    }


    return 0;

}


void help() {
    cout << "RNASeq maps short sequences from spliced transcripts to whole genomes.\n\n";


    cout << "Usage:\n";
    cout << "\t" << "./rnaseq.out [options]";

    cout << "\n\nOptions (Common):\n";
    cout << "\t" << "-r/--read\t"         << "<filename>\t"    << "[ default: 30           ]\n";
    cout << "\t" << "-fa/--genome\t"      << "<filename>\t"    << "[ default: 30           ]\n";
    cout << "\t" << "-x/--index\t"        << "<filename>\t"    << "[ default: 30           ]\n";

    cout << "\t" << "-l1/--chunk-size\t"  << "<int>\t"         << "[ default: 30           ]\n";
    cout << "\t" << "-v/--dash-v\t"       << "<int>\t"         << "[ default: 3            ]\n";
    cout << "\t" << "-ng/--num-gap\t"     << "<int>\t"         << "[ default: 5            ]\n";
    cout << "\t" << "-p/--thread\t"       << "<int>\t"         << "[ default: 1            ]\n";
    cout << "\t" << "-o/--output-dir\t"   << "<string>\t"      << "[ default: ./rnaseq_out ]\n";
    cout << "\t" << "-h/--help\n";

    cout << "\n\nOptions (Anchoring):\n";
    cout << "\t" << "-ansh/--an-shift\t"  << "<int>\t"         << "[ default: 5            ]\n";

    cout << "\n\nOptions (Assignment):\n";
    cout << "\t" << "-s/--is-slide\t"     << "<0-1>\t"         << "[ default: 1            ]\n";
    cout << "\t" << "-assh/--as-shift\t"  << "<int>\t"         << "[ default: 10            ]\n";

/*Options:
    -a/--min-anchor                <int>       [ default: 8                ]
    -m/--splice-mismatches         <0-2>       [ default: 0                ]
    --insertions                   <filename>*/
}
