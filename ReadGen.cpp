#include <iostream>
#include <cstdlib>
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <math.h>
#include <ctime>
#include <random>
#include <chrono>

using namespace std;

int NUM_THREADS;
int READS_PER_THREAD;
int ReverseShare;
string inFasta;

string complement(string read){

     string temp = "";
     for (int i=0;i<read.size();i++){
	switch	(read[i]){
	case 'A':
		temp+='T';
		break;
	case 'T':
		temp+='A';
		break;
	case 'C':
		temp+='G';
		break;
	case 'G':
		temp+='C';
		break;
	case 'a':
		temp+='t';
		break;
	case 't':
		temp+='a';
		break;
	case 'c':
		temp+='g';
		break;
	case 'g':
		temp+='c';
		break;
	default:
		temp += 'N';
	}
    }
	return temp;
}
	
string NumToStr(int number)
{
    ostringstream ostr;
    ostr<<number;
    return ostr.str();
}

string LongToStr(long long number)
{
	ostringstream ostr;
	ostr<<number;
	return ostr.str();
}

string whole_genome;

void readFASTA(string filename)
{
	ifstream infile(filename.c_str());
	string line;
	while(getline(infile, line))
	{
		if (line[0]=='>')
			continue;
		whole_genome += line;		
	}
	cout << "Genome read." << endl;
}

struct thread_data
{
   string ref;
   float PM;
   float PI;
   int frag_num;
   int READ_LENGTH;
};

void *noise(void *threadarg)
{
   struct thread_data *my_data;
   my_data = (struct thread_data *) threadarg;
   string ref = my_data->ref;
   float PI = my_data->PI;
   float PM = my_data->PM;
   int frag_num = my_data->frag_num;
   int READ_LENGTH = my_data->READ_LENGTH;
   //long long count = floor(ref.size() * log(ref.size())/(READ_LENGTH*log(2)));
   long long count = READS_PER_THREAD;
   vector<long long> begin;
   begin.resize(count);
   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
   default_random_engine generator (seed);
   uniform_real_distribution<double> begindist(0,1);  
   string filename = "reads" + NumToStr(frag_num) + ".fq";
   ofstream outfile(filename.c_str()); 
   string quality = "\n+\n";
   for (int i=0;i<READ_LENGTH;i++)
	quality += "x";
   long long frag_index = frag_num*floor(whole_genome.size()/(NUM_THREADS));
   long long avg_len = floor(whole_genome.size()/(NUM_THREADS));
   //long long previous_count = frag_num*floor(avg_len*log(avg_len)/(READ_LENGTH*log(2)));
   long long previous_count = frag_num*READS_PER_THREAD;
   for (long long i=0;i<count;i++)
   {
	   double a = begindist(generator);	
	   int Ncount = 0;
	   long long begin = floor(a*(ref.size() - 2*READ_LENGTH)+READ_LENGTH);
	   string noisy_frag;
	   string base = "ATCG";
	   long long indel = 0;
	   long long insertion = 0;
	   long long deletion = 0;
	   long long mismatch = 0;
	   long long b_ref = 0;
	   long long b_tar = 0;
	   long long ref_length = ref.size();
	   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	   default_random_engine generator (seed);
	   uniform_real_distribution<double> unidist(0,100);
	   uniform_real_distribution<double> unidist2(0,100);
	   uniform_int_distribution<int> intdist3(0,3);
	   uniform_int_distribution<int> intdist2(0,2);
	   uniform_real_distribution<double> ReverseOrNot (0,100);
	   double RorF = ReverseOrNot(generator);
	   
		while (b_tar < READ_LENGTH)
		{
			int flag = 1;
			double U = unidist(generator);
			if (U < PI)
			{
				double V = unidist2(generator);
				indel ++;
				if (V < 49.8)
				{
					deletion ++;
					b_ref ++;
					flag = 0;
				}
				else
				{	
					insertion ++;
					noisy_frag.push_back(base[intdist3(generator)]);
				}
			}		
			else
			{
				double W = unidist(generator);
				if (W < PM)
				{
					mismatch ++;
					string temp;
					for (int j=0;j<base.size();j++)
					{
						if (base[j] != toupper(ref[b_ref+begin]))
							temp.push_back(base[j]);
					}
					char newbase = temp[intdist2(generator)];
					noisy_frag.push_back(newbase);
	    				b_ref ++;
				}
				else
				{
					noisy_frag.push_back(ref[b_ref+begin]);
					if (ref[b_ref+begin] == 'N')
						Ncount++;
					b_ref ++;
				}

			}
			if (flag)
				b_tar ++;
		}
	    if (RorF < ReverseShare){
		noisy_frag = complement(noisy_frag);
		reverse(noisy_frag.begin(),noisy_frag.end());
	    }
	    if ((double) Ncount/READ_LENGTH > 0.15)
			i--;
	    else
			outfile<<"@r"<<LongToStr(begin+frag_index)<<"_"<<LongToStr(i+previous_count)<<"^"<<"\n"<<noisy_frag<<quality<<endl;
    }
    outfile.close();
    pthread_exit(NULL);
}

int main (int argc, char* argv[])
{
   
   double mismatch, indel;
   int READ_LENGTH;
   int total_reads;
   
   
   for (int i = 1; i < argc; )
	{
		if (string(argv[i]) == "-readlen")
		{READ_LENGTH = atoi(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-mismatch")
		{mismatch = stod(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-indel")
		{indel = stod(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-p")
		{NUM_THREADS = atoi(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-num")
		{total_reads = atoi(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-rev")
		{ReverseShare = atoi(argv[i + 1]); i += 2;
			if (i==argc)
				break;}
		if (string(argv[i]) == "-fa")
		{inFasta = (argv[i + 1]); i += 2;
			if (i==argc)
				break;}
	}
	
	READS_PER_THREAD = total_reads / NUM_THREADS;
   
   cout << "start : reading noisy genome..."<<endl;
   readFASTA(inFasta.c_str());
   
   pthread_t threads[NUM_THREADS];
   pthread_attr_t attr;
   struct thread_data td[NUM_THREADS];
   int rc;
   int i;
   void *status;

   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   long long frag = floor(whole_genome.size()/(NUM_THREADS));
   for( i=0; i < NUM_THREADS; i++ ){
      td[i].ref = whole_genome.substr(i*frag,frag);
      if (i == NUM_THREADS-1)
	td[i].ref = whole_genome.substr(i*frag,whole_genome.size());
      td[i].PM = mismatch;
      td[i].PI = indel;
      td[i].frag_num = i;
      td[i].READ_LENGTH = READ_LENGTH;
      rc = pthread_create(&threads[i], NULL, noise, (void *)&td[i]);
      if (rc)
      {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
      }
   }
   
   pthread_attr_destroy(&attr);
   for( i=0; i < NUM_THREADS; i++ )
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
   string command = "cat ";
   for (int i=0;i<NUM_THREADS;i++)
	command += "reads" + NumToStr(i) + ".fq ";
   command += "> reads.fq";
   cout << "Now concatenating read files..." << endl;
   system(command.c_str());
   cout << "Main: program exiting." << endl;
   pthread_exit(NULL); 
}
