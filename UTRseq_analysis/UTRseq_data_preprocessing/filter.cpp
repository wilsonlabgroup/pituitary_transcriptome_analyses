//
//  main.cpp
//  test
//
//  Created by QuRihao on 7/14/16.
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


int judge(string &s, string match_sequence[], unsigned long len[], unsigned long l, int number)
{
    
    int start_pos_bait = int(l-10)/2;
    string bait = s.substr(start_pos_bait, 10);
    int flag_bait = 0;
    int potential = 0;
    for (int k = 0; k < number; k++)
    {
        if (match_sequence[k].find(bait, 0) < len[k])
        {
            flag_bait = 1;
            potential = k;
            break;
        }
        else continue;
    }
    
    if (flag_bait == 1)
    {
        int start_pos = int(l-24)/2;
        string t = s.substr(start_pos-4, 24);
        int flag = 0;
        int mark = 0;
        for (int i = potential; i < number; i++)
        {
            if (match_sequence[i].find(t, 0) < len[i])
            {
                flag = 1;
                mark = i+1;
                break;
            }
            else continue;
        }
        
        if (flag == 1)
        {
            return mark;
        }
        else return 0;
    }
    
    else return 0;
}

inline bool string_convert(string & a, char A[], int len_array=100)
{
    if (len_array<a.length())
    {
        //return false;
        throw ("Invalid format change between string and char array!");
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

int main(int argc, char* argv[])
{
    
    ////////////////////////////////////////////////////////////////////////////////
    //////////////Receive and initialize parameters/////////////////////////////////
    
    
    string in_file = "";
    string outdir = "";
    string pre_dir = "2-filter";
    string filt_file = "";
    
    int number = 100;
    int MINIMAL_LENGTH = 36;
    int existence_dir = -1;
    unsigned long lines = 0;
    
    //receive parameters from command line
    int ch;
    opterr = 0; //do not declare any error
    while ((ch = getopt(argc, argv, "d:i:e:f:l:")) != -1)
    {
        switch (ch)
        {
            case 'i': in_file = string(optarg); cout << "Input fastq file: " << in_file << endl; break;
            case 'f': filt_file = string(optarg); cout << "Sequence file used as matched reference: " << filt_file << endl; break;
            case 'd': outdir = string(optarg); cout << "Output directory: " << outdir << endl; break;
            case 'e': existence_dir = atoi(optarg); cout << "Whether output directory need to be built: " << existence_dir << endl; break;
            case 'l': MINIMAL_LENGTH = atoi(optarg); cout << MINIMAL_LENGTH << endl; break;
        }
    }
    cout << endl;
    
    if (in_file.length() == 0 || filt_file.length() == 0)
    {
        cout << "At least TWO parameters should be provided." << endl;
        cout << "Example: ./filter -i filename.fastq -f match_sequence.fa [-l minimal_length] [-d outdir [-e dir_existence]]" << endl;
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
            cout << "Example: ./filter -i filename.fastq -f match_sequence.fa [-l minimal_length] [-d outdir [-e dir_existence]]" << endl;
            exit(0);
        }
    }
    
    string mkdir = "mkdir "+pre_dir;
    system_process(mkdir);
    
    
    //process.....
    char filename[100];
    string_convert(in_file, filename);
    ifstream in;
    in.open(filename);
    
    char matchfile_name[100];
    string_convert(filt_file, matchfile_name);
    ifstream matchfile;
    matchfile.open(matchfile_name);
    
    char main_output_name[100];
    string main_output = pre_dir+"/filtered.fastq";
    string_convert(main_output, main_output_name);
    ofstream out;
    out.open(main_output_name);
    
    char removed_output_name[100];
    string removed_output = pre_dir+"/removed.fastq";
    string_convert(removed_output, removed_output_name);
    ofstream rm;
    rm.open(removed_output_name);
    
    char log_output_name[100];
    string log_output = pre_dir+"/filter.log";
    string_convert(log_output, log_output_name);
    ofstream log;
    log.open(log_output_name);
    
    
    
    ////////////////////////////////////////////////////////////////////////////////
    ////Receive matched fasta file as the filtered reference////////////////////////
    
    string * match_sequence = new string[number];
    for (int i = 0; i < number; i++) match_sequence[i] = "";
    string * sequence_name = new string[number];
    for (int i = 0; i < number; i++) sequence_name[i] = "";
    unsigned long * len = new unsigned long[number];
    for (int i = 0; i < number; i++) len[i] = 0;
    
    string line;
    
    int timer = 0;
    
    while (getline (matchfile, line))
    {
        unsigned long lg = line.length();
        char marker = line[0];
        if (marker == '>')
        {
            if (lg == 1)
            {
                cout << "Please provide name for all the sequences." << endl;
                exit(0);
            }
            else
            {
                if (line[1]==' ') sequence_name[timer] = line.substr(2, lg-2);
                else sequence_name[timer] = line.substr(1, lg-1);
                cout << "Receive sequence: " << sequence_name[timer] << endl;
            }
        }
        
        else if (marker=='A' || marker=='G' || marker=='T' || marker=='C')
        {
            match_sequence[timer] = line;
            len[timer] = line.length();
            
            timer++;
        }
        else continue;
    }
    cout << endl;
    
    
    ////////////////////////////////////////////////////////////////////////////////
    ///////////Create corresponding-named filtered output files/////////////////////
    
    if (timer > 0) number = timer;
    
    long * counter = new long[number];
    for (int i = 0; i < number; i++) counter[i] = 0;
    
    ofstream output[100];
    string suffix = ".fastq";
    for (int i = 0; i < number; i++)
    {
        char * name = new char[100];
        
        string prename = pre_dir+"/"+sequence_name[i]+suffix;
        string_convert(prename, name);
        for (int j = 0; j<100; j++)
        {
            if (name[j] == ' ') name[j] = '_';
        }
        output[i].open(name);
        delete [] name;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////
    //////////////////////////PROCESS FASTQ/////////////////////////////////////////
    
    long counter_del = 0;
    int count = 0;
    int bit = 0;
    unsigned long l = 0;
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
            cout << "Fastq file line number: " << lines << endl;
            cout << endl;
            cout << "Successfully submitted, processing..." << endl;
            
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
        unsigned long compare[100];
        for (int i=0; i<100; i++)
        {
            compare[i]=(lines/100)*(i+1);
        }
        
        while (getline (in, line)) // '\n' is not included
        {
            switch (count)
            {
                case 0:
                    tmp1 = line; break;
                case 1:
                    tmp2 = line; break;
                case 2:
                    tmp3 = line; break;
                case 3:
                    tmp4 = line; break;
            }
            
            if (count == 3)
            {
                l = tmp2.length();
                
                if (l >= MINIMAL_LENGTH)
                {
                    bit = judge(tmp2, match_sequence, len, l, number);
                    log << bit <<endl;
                    
                    if (bit == 0)
                    {
                        out << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                    }
                    else
                    {
                        output[bit-1] << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                        counter[bit-1] += 1;
                    }
                }
                else
                {
                    rm << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                    counter_del += 1;
                }
            }
            
            count = (count+1)%4;
            
            
            if (process_timer >= compare[report_timer])
            {
                cout << "Processed "<< (report_timer+1) << "%" <<endl;
                report_timer++;
            }
            
            process_timer++;
        }
        
        cout << "Filter result:" <<endl;
        log << "Filter result:" <<endl;
        for (int i = 0; i < number; i ++)
        {
            cout << sequence_name[i] << " : " << counter[i] << endl;
            log << sequence_name[i] << " : " << counter[i] << endl;
        }
        
        cout << "Deleted sequence: " << counter_del << endl;
        log << "Deleted sequence: " << counter_del << endl;
    }
    
    else // file doesn't exist
    {
        cout <<"No such file!" << endl;
    }
    
    for (int i = 0; i < number; i++) output[i].close();
    in.close();
    out.close();
    rm.close();
    log.close();
    matchfile.close();
    
    return 0;
}

