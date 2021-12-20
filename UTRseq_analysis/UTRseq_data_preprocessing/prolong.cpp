//
//  main.cpp
//  gene_extractor
//
//  Created by QuRihao on 7/29/16.
//  Copyright © 2016 QuRihao. All rights reserved.
//

#include <iostream>
#include <cstdlib> //atoi, atof
#include <string>
#include <cstring>
#include <fstream>
#include <unistd.h>
#include <getopt.h>

using namespace std;

string KERNEL = "";//>=7
int LEN = 0; //KERNEL size
const int MAX_ARRAY_LENGTH = 2000;
int ITERATE = 1;
int ITERATE_TIMER = 0;
long COUNTER = 0;
long** RECORD = new long*[MAX_ARRAY_LENGTH]; //100~2000
//int SIZE = 100;
int THRESHOLD[5] = {15, int('0'), 75, int('5'), 500};
int start_KERNEL = 1000;
int stop_KERNEL = 1000;
string left_end_KERNEL = "";
string right_end_KERNEL = "";
string S = "";
string T = "";
int LEFT_MOVE = 0;
int START_CHECK_POSITION;
int STOP_CHECK_POSITION;
int LENGTH = 0;
int LEVEL = 1;
int *LABEL = new int[MAX_ARRAY_LENGTH];

inline bool counter(const int i, const int &position)
{
    if (T[i] > THRESHOLD[3])
    {
        switch (S[i])
        {
            case 'A':
            {
                RECORD[position][0]+=LEVEL;
                break;
            }
            case 'T':
            {
                RECORD[position][1]+=LEVEL;
                break;
            }
            case 'C':
            {
                RECORD[position][2]+=LEVEL;
                break;
            }
            case 'G':
            {
                RECORD[position][3]+=LEVEL;
                break;
            }
        }
    }
    return true;
}

inline bool double_check(int move_s)
{
    for (int i=START_CHECK_POSITION; i<=STOP_CHECK_POSITION; i++)
    {
        if (!(S[i+move_s]==KERNEL[i] || T[i+move_s]<THRESHOLD[1])) return false;
    }
    return true;
}

bool left_match(int left_search)
{
    LEFT_MOVE = start_KERNEL-left_search;
    START_CHECK_POSITION = THRESHOLD[0];
    STOP_CHECK_POSITION = ((LENGTH-left_search)>LEN)?(LEN-1):(LENGTH-1-left_search);
    if (START_CHECK_POSITION <= STOP_CHECK_POSITION && double_check(left_search))
    {
        for (int i=0; i<left_search; i++)
        {
            if (counter(i, LEFT_MOVE+i)) continue;
        }
        return true;
    }
    return false;
}

bool right_match(int right_search)
{
    LEFT_MOVE = stop_KERNEL+1-THRESHOLD[0]-right_search;
    START_CHECK_POSITION = ((right_search + THRESHOLD[0]) > LEN)? 0:(LEN-THRESHOLD[0]-right_search);
    STOP_CHECK_POSITION = LEN-1-THRESHOLD[0];
    if (STOP_CHECK_POSITION >= START_CHECK_POSITION && double_check(right_search-(LEN-THRESHOLD[0])))
    {
        for (int i = right_search+THRESHOLD[0]; i<LENGTH; i++)
        {
            if (counter(i, LEFT_MOVE+i)) continue;
        }
        return true;
    }
    return false;
}

bool prolong()
{
    LENGTH = int(S.length());
    size_t left_search = S.find(left_end_KERNEL, 0);
    size_t right_search = S.find(right_end_KERNEL, 0);
    bool x = (left_search < LENGTH);
    bool y = (right_search < LENGTH);
    if (!(x || y))
    {
        return false;
    }
    else if (x && !y)
    {
        LEVEL = 1;
        if (left_match(int(left_search))) return true;
    }
    else if (!x && y)
    {
        LEVEL = 1;
        if (right_match(int(right_search))) return true;
    }
    else
    {
        LEVEL = 4;
        if (left_match(int(left_search)) && right_match(int(right_search))) return true;
    }
    return false;
}

inline char exchange(int a)
{
    switch (a)
    {
        case 0:
        {
            return 'A';
            break;
        }
        case 1:
        {
            return 'T';
            break;
        }
        case 2:
        {
            return 'C';
            break;
        }
        case 3:
        {
            return 'G';
            break;
        }
    }
    return 'N';
}

inline bool string_convert(string & a, char A[], int len_array=100)
{
    if (len_array<a.length())
    {
        cout << "Invalid format change between string and char array!" << endl;
        exit(0);
    }
    else
    {
        strcpy(A, a.c_str());
        return true;
    }
}

inline bool system_process(string &a)
{
    if (a.length()>0)
    {
        char * CMD = new char[200];
        if (string_convert(a, CMD, 200)) system(CMD);
        delete [] CMD;
        return true;
    }
    else return false;
}


int main(int argc, char * argv[]) {
    
    string in_file = "";
    string outdir = "";
    string pre_dir = "gene-extract";
    string filt_file = "";
    string out_name = "extract.fa";
    string log_name = "extract.log";
    
    //int number = 100;
    int existence_dir = -1;
    unsigned long lines = 0;
    
    //receive parameters from command line
    int ch;
    opterr = 0; //do not declare any error
    while ((ch = getopt(argc, argv, "d:i:e:s:l:t:n:p:g:")) != -1)
    {
        switch (ch)
        {
            case 'i': in_file = string(optarg); cout << "Input fastq file: " << in_file << endl; break;
            //case 'f': filt_file = string(optarg); cout << "Sequence file used as matched reference: " << filt_file << endl; break;
            case 's': KERNEL = string(optarg); cout << "Kernel sequence: " << KERNEL << endl; break;
            case 'd': outdir = string(optarg); cout << "Output directory: " << outdir << endl; break;
            case 'e': existence_dir = atoi(optarg); cout << "Whether output directory need to be built: " << existence_dir << endl; break;
            case 't': ITERATE = atoi(optarg); cout << "Iteration times: " << ITERATE << endl; break;
            case 'l': THRESHOLD[2] = atoi(optarg); cout << "Confidence level: " << THRESHOLD[2] << "%" << endl; break;
            case 'n': THRESHOLD[4] = atoi(optarg); cout << "Minimal count requirement: " << THRESHOLD[4] << endl; break;
            case 'p': out_name = string(optarg); cout << "File name for output sequence: " << out_name <<endl; break;
            case 'g': log_name = string(optarg); cout << "File name for log: " << log_name <<endl; break;
        }
    }
    cout << endl;
    
    
    if (in_file.length() == 0 || KERNEL.length() == 0)
    {
        cout << "At least TWO parameters should be provided." << endl;
        cout << "Format: ./prolong -i filename.fastq -s KERNEL_SEQUENCE [-t Iterate_times] [-l Confidence_level (default 75%)] [-n Minimal_count_requirement (default >500)] [-d outdir [-e dir_existence]] [-p file name] [-g log file name]" << endl;
        exit(0);
    }
    
    if (outdir.length()>0)
    {
        if (existence_dir == -1 || existence_dir > 0)
        {
            string a = "mkdir "+outdir;
            system_process(a);
            pre_dir = outdir+"/"+pre_dir;
        }
        else pre_dir = outdir+"/"+pre_dir;
    }
    
    else
    {
        if (existence_dir > 0)
        {
            cout << "New directory name must be provided." << endl;
            cout << "Format: ./prolong -i filename.fastq -s KERNEL_SEQUENCE [-t Iterate_times] [-l Confidence_level (default 75%)] [-n Minimal_count_requirement (default >500)] [-d outdir [-e dir_existence]] [-p file name] [-g log file name]" << endl;
            exit(0);
        }
    }
    
    string mkdir = "mkdir "+pre_dir;
    system_process(mkdir);
    
    
    //process.....
    
    char filename[100];
    string_convert(in_file, filename);
    
    char searchfile_name[100];
    string_convert(filt_file, searchfile_name);
    ifstream searchfile;
    searchfile.open(searchfile_name);
    
    char main_output_name[100];
    string main_output = pre_dir+"/"+out_name;
    string_convert(main_output, main_output_name);
    ofstream out;
    out.open(main_output_name);
    
    char log_output_name[100];
    string log_output = pre_dir+"/"+log_name;
    string_convert(log_output, log_output_name);
    ofstream log;
    log.open(log_output_name);
    
    LEN = int(KERNEL.length());
    string line;
    
    int left_half = MAX_ARRAY_LENGTH/2;
    
    
    for (ITERATE_TIMER=0; ITERATE_TIMER<ITERATE; ITERATE_TIMER++)
    {
        ifstream in;
        in.open(filename);
        
        COUNTER = 0;
        
        for (int i=0; i<MAX_ARRAY_LENGTH; i++)
        {
            LABEL[i] = -1;
            RECORD[i] = new long[4];
            for (int j=0; j<4; j++)
            {
                RECORD[i][j] = 0;
            }
        }
        
        start_KERNEL = left_half-LEN/2;
        stop_KERNEL = start_KERNEL+LEN-1;
        
        
        left_end_KERNEL=KERNEL.substr(0, THRESHOLD[0]);
        right_end_KERNEL = KERNEL.substr(LEN-THRESHOLD[0], THRESHOLD[0]);
        
        
        int count = 0;
        
        string tmp1;
        string tmp2;
        string tmp3;
        string tmp4;
        
        if (in) // file exists
        {
            string b = pre_dir+"/tmp1";
            string c = "wc -l "+in_file+" >"+ b;
            system_process(c);
            
            char tmp_name[100];
            string_convert(b, tmp_name);
            ifstream tmp_read;
            tmp_read.open(tmp_name);
            
            if (getline(tmp_read, line))
            {
                char number_line[100];
                string_convert(line, number_line);
                sscanf(number_line,"%lu",&lines);
                //cout << "Successfully submitted, processing..." << endl;
                
                string d = "rm "+b;
                system_process(d);
            }
            else
            {
                cout << "No write permission." << endl;
                exit(0);
            }
            
            unsigned long process_timer = 0;
            int report_timer = 0;
            unsigned long compare[10];
            for (int i=0; i<10; i++) compare[i]=(lines/10)*(i+1);
            
            
            while (getline (in, line)) // '\n' is not included
            {
                switch (count)
                {
                    case 0:
                        tmp1 = line; break;
                    case 1:
                        S = line; break;
                    case 2:
                        tmp3 = line; break;
                    case 3:
                        T = line; break;
                }
                
                if (count == 3 && prolong()) COUNTER++;
                
                count = (count+1)%4;
                
                
                if (process_timer >= compare[report_timer])
                {
                cout << "ITERATION-" << ITERATE_TIMER+1 << " | Processed "<< 10*(report_timer+1) << "%" <<endl;
                report_timer++;
                }
                
                process_timer++;
            }
            
            cout << "ITERATION-" << ITERATE_TIMER+1 << " | Terminate and processing..." << endl;
            cout << "Counter:" << COUNTER << endl;
            log << "Counter:" << COUNTER <<endl;
            
            int left_edge = start_KERNEL;
            int right_edge = stop_KERNEL;
            
            for (int i=start_KERNEL-1; i>=0; i--)
            {
                cout << "Position "<<i<<": ";
                long scale = 0;
                long sum = 0;
                for (int j=0; j<4; j++)
                {
                    cout << exchange(j) << "-count: "<<RECORD[i][j]<<" |";
                    if (RECORD[i][j]>scale && RECORD[i][j]>THRESHOLD[4])
                    {
                        LABEL[i] = j;
                        scale = RECORD[i][j];
                    }
                    sum += RECORD[i][j];
                }
                cout << endl;
                if (scale > THRESHOLD[2]*sum/100 && scale > 0)
                {
                    left_edge--;
                    //cout << "Position " << i-start_KERNEL << ": "<< scale << "/" << sum << endl;
                }
                else break;
            }
            
            for (int i=stop_KERNEL+1; i<MAX_ARRAY_LENGTH; i++)
            {
                cout << "Position "<<i<<": ";
                long scale = 0;
                long sum = 0;
                for (int j=0; j<4; j++)
                {
                    cout << exchange(j) << "-count: "<<RECORD[i][j]<<" |";
                    if (RECORD[i][j]>scale && RECORD[i][j]>THRESHOLD[4])
                    {
                        LABEL[i] = j;
                        scale = RECORD[i][j];
                    }
                    sum += RECORD[i][j];
                }
                cout << endl;
                if (scale > THRESHOLD[2]*sum/100 && scale > 0)
                {
                    right_edge++;
                    //cout << "Position " << i-start_KERNEL << ": "<< scale << "/" << sum << endl;
                }
                else break;
            }
            
            for (int i=start_KERNEL-1; i>=left_edge; i--)
            {
                switch (LABEL[i])
                {
                    case 0:
                    {
                        KERNEL = "A"+KERNEL;
                        break;
                    }
                    case 1:
                    {
                        KERNEL = "T"+KERNEL;
                        break;
                    }
                    case 2:
                    {
                        KERNEL = "C"+KERNEL;
                        break;
                    }
                    case 3:
                    {
                        KERNEL = "G"+KERNEL;
                        break;
                    }
                }
            }
        
            for (int i=stop_KERNEL+1; i<=right_edge; i++)
            {
                switch (LABEL[i])
                {
                    case 0:
                    {
                        KERNEL = KERNEL+"A";
                        break;
                    }
                    case 1:
                    {
                        KERNEL = KERNEL+"T";
                        break;
                    }
                    case 2:
                    {
                        KERNEL = KERNEL+"C";
                        break;
                    }
                    case 3:
                    {
                        KERNEL = KERNEL+"G";
                        break;
                    }
                }
            }
            int original_LEN = LEN ;
            LEN = int(KERNEL.length());
            cout << "LENGTH change from " << original_LEN << "nt to " << LEN << "nt." << endl;
            
            out << "ITERATION-" << ITERATE_TIMER+1 << ":" << endl;
            out << KERNEL << endl;
            cout << endl;
        }
        
        else // file doesn't exist
        {
            cout <<"No such file!" << endl;
        }
        
        in.close();
        
        for (int i=0; i<MAX_ARRAY_LENGTH; i++)
        {
            delete RECORD[i];
        }
    }
    
    out.close();
    log.close();
    searchfile.close();
    
    return 0;
}
