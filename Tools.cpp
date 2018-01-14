//
//  Tools.cpp
//  RNA-Seq
//
//  Created by Farid Rashidi on 8/6/16.
//  Copyright © 2016 Farid Rashidi. All rights reserved.
//

#include "Tools.h"

string convertNumToStr(long long number) {
    ostringstream ostr;
    ostr<<number;
    return ostr.str();
}

void pF(string s) {
    FILE *myOutput;
    myOutput = fopen("log.txt", "a");
    fprintf(myOutput, "%s", s.c_str());
    fclose(myOutput);
}

void pFS(string s) {
    FILE *myOutput;
    myOutput = fopen("log.txt", "a");
    fprintf(myOutput, "%s", s.c_str());
    fclose(myOutput);
    fprintf(stdout, "%s", s.c_str());
}

string get_current_time() {
    time_t rawtime;
    struct tm * timeinfo;

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    char output[30];
    strftime(output, 30, "[%Y-%m-%d %H:%M:%S]", timeinfo);

    return string(output);
}
