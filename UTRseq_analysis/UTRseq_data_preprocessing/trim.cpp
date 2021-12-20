//
//  trim.cpp
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

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////FUNCTION DEFINITION/////////////////////////////////

int judge_position(const string &s, const string &t, const int start_backsearch_site, const int threshold[])
{
    int recorder = start_backsearch_site; //record the positon of poly'A' START site
    int subrecorder = 0;
    int flag = 0;
    
    for (int i = start_backsearch_site-1; i >= 0; i--)
    {
        if (s[i]=='A')
        {
            recorder--;
            if (flag > 0) subrecorder++;
        }
        else
        {
            if (flag < 2)
            {
                if (t[i]<threshold[2]) flag+=1;
                else flag+=2;
                recorder--;
                subrecorder++;
            }
            else break;
        }
    }
    
    int totrim = start_backsearch_site-recorder;
    if (flag > 0)
    {
        if (subrecorder > max(3, totrim/2)) return recorder;
        else
        {
            recorder += subrecorder;
            return recorder;
        }
    }
    else return recorder;
}

int quality_control(const string &t, const int &start_backsearch_site, const int threshold[], const int &stop_backsearch_site=0)
{
    if (start_backsearch_site <= stop_backsearch_site) return -1;
    
    int recorder = start_backsearch_site; //record the positon of cut START site
    int subrecorder = 0;
    int subsubrecorder = 0;
    int flag = 0;
    
    for (int i = start_backsearch_site-1; i >= stop_backsearch_site; i--)
    {
        if (t[i]<threshold[7])
        {
            recorder--;
            if (flag > 0) subrecorder++;
            if (flag == 2) subsubrecorder++;
        }
        else
        {
            if (flag < 2)
            {
                flag+=1;
                recorder--;
                subrecorder++;
                if (flag == 1) subsubrecorder++;
            }
            else break;
        }
    }
    
    if (flag > 0)
    {
        if (subrecorder > 4)
        {
            if (subsubrecorder > 2) return recorder;
            else
            {
                recorder += subsubrecorder;
                return recorder;
            }
        }
        else
        {
            recorder += subrecorder;
            return recorder;
        }
    }
    else return recorder;
}


bool start_adapter_match(const string &s, const string &t, const string adapter[], const int threshold[], const int adapter_number)
{
    int score = 0;
    for (int i=0; i<adapter_number; i++)
    {
        int len = int(adapter[i].length());
        if (len <= threshold[9])
        {
            cout << "Invalid adapter sequence, too short (<=" << threshold[9] << ")!" << endl;
            exit(0);
        }
        else
        {
            int max_length_match = min(len, int(s.length()));
            for (int j=threshold[9]+1; j< max_length_match; j++)
            {
                score = 0;
                for (int k=0; k<j; k++)
                {
                    if (s[k]==adapter[i][len-j+k] || t[k]<threshold[2]) score++;
                }
                if (100*score/j > threshold[8]) return true;
            }
        }
    }
    return false;
}
                         

bool maximal_A_length_search(const string &s, const string &t, int record[], const int threshold[], int &marker_point)
{
    //record[0]: length of s
    //record[2N-1]: length of continuous 'A'
    //record[2N]: the positon of poly'A' START site
    
    int len = record[0];
    record[1] = 0; // to reset the value
    int counter = 0;
    int N = 1;
    int flag = 0;
    int tmp = -1;
    
    for (int i=len-1; i>=tmp; i--)
    {
        if (i!=-1 && (s[i]=='A' || t[i]<threshold[2]))
        {
            counter++;
        }
        else
        {
            if (counter>=max(1,record[1]) || counter>threshold[6])
            {
                if (counter>threshold[6] && record[1]>threshold[6]) N++;
                record[2*N-1] = counter;
                record[2*N] = i+1;
                flag = min(counter,threshold[6]);
            }
            counter = 0;
        }
        tmp = ((flag-counter-1)>-1)?(flag-counter-1):-1;
    }
    
    marker_point = N;
    
    return true;
}


bool match_accuracy(const string &s, const string &t, const string adapter[], const int threshold[], const int start_site, const int scale, const int adapter_number, const int search_limit[])
{
    //threshold[0]: the minimal pairing length of overlapping area for terminal site adapter matching
    //threshold[1]: initial parameter of required matching accuracy
    //threshold[2]: minimal sequencing quality of each base while pairing and searching for poly-A
    //threshold[3]: maximal searching diameter
    //threshold[4]: maximal matching length
    
    int counter = 0;
    int uplimit = 0;
    if (scale >= threshold[0])
    {
        for (int i=0; i<adapter_number; i++)
        {
            uplimit = (scale > search_limit[i])?search_limit[i]:scale;
            
            for (int j=0; j<uplimit; j++)
            {
                if (t[start_site+j]<threshold[2] || s[start_site+j]==adapter[i][j]) counter++;
            }
            if ((100*counter)/uplimit>threshold[1]) return true;
            
            counter = 0;
        }
    }
    else
    {
        int adjusted_accuracy_threshold=((threshold[0]-scale)*(100-threshold[1])/threshold[0])+threshold[1];
        for (int i=0; i<adapter_number; i++)
        {
            for (int j=0; j<scale; j++)
            {
                if (t[start_site+j]<threshold[2] || s[start_site+j]==adapter[i][j]) counter++;
            }
            if ((100*counter)/scale>adjusted_accuracy_threshold) return true;
            
            counter = 0;
        }
    }
    return false;
}


int judge_adapter(const string &s, const string &t, const int record[], const string adapter[], const int threshold[], const int adapter_number, const int search_limit[], const int marker_point, const int offset)
{
    //threshold[0]: the minimal pairing length of overlapping area for terminal site adapter matching
    //threshold[1]: initial parameter of required matching accuracy
    //threshold[2]: minimal sequencing quality of each base while pairing and searching for poly-A
    //threshold[3]: maximal searching diameter
    //threshold[4]: maximal matching length
    //threshold[5]: if contain continuous 'A' chain longer than threshold[5], then set as a cutting point
    //record[0]: length of s
    //record[2N-1]: length of continuous 'A'
    //record[2N]: the positon of poly'A' START site
    
    int sign = 1;
    int diameter = 0;
    int term_pos = 0;
    int scale = 0;
    int maximal_searching_times = 2*threshold[3]+1;
    int start_search_site = 0;

    for (int N=marker_point; N>0; N--)
    {
        if (record[2*N-1]>threshold[5]) return N;
        else
        {
            term_pos = record[2*N]+record[2*N-1]-1;
            start_search_site = term_pos+1-offset;
            diameter = 0;
            
            for (int i=0; i<maximal_searching_times; i++)
            {
                start_search_site += (sign*diameter);
                scale = record[0]-start_search_site;
                if (scale<=record[0])
                {
                    if (scale>0)
                    {
                        if (match_accuracy(s, t, adapter, threshold, start_search_site, scale, adapter_number, search_limit)) return N;
                        else
                        {
                            diameter++;
                            sign*=(-1);
                        }
                    }
                    else if (diameter>0)
                    {
                        if (record[2*N-1]>max(2,(diameter+diameter%2))) return N; //minimal length of embedded poly-A tail is 3
                    }
                    else return N;
                }
            }
        }
    }
    return 0;
}


inline bool string_convert(string & a, char A[], int len_array=100)
{
    if (len_array<a.length())
    {
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
    
    string in_file;
    string adapter_file_name;
    string outdir;
    string pre_dir = "1-trim";
    int existence_dir = -1;
    int HEADCUTOFF = 12;
    int lines = 0;
    int threshold[10] = {5, 50, int('0'), 4, 20, 25, 3, int('6'), 80, 8};
    //threshold[0]: the minimal pairing length of overlapping area for terminal site adapter matching
    //threshold[1]: initial parameter of required matching accuracy
    //threshold[2]: minimal sequencing quality of each base while pairing and searching for poly-A
    //threshold[3]: maximal searching diameter
    //threshold[4]: maximal matching length
    //threshold[5]: if contain continuous 'A' chain longer than threshold[5], then set as a cutting point
    //threshold[6]: if contain continuous 'A' chain longer than threshold[6], then push into waiting stack
    //threshold[7]: minimal sequencing quality of each base under quality control
    //threshold[8]: confidence level for start site adapter matching
    //threshold[9]: the minimal pairing length of overlapping area for start site adapter matching
    
    
    int ch;
    opterr = 0; //do not declare any error
    while ((ch = getopt(argc, argv, "d:i:e:a:c:")) != -1)
    {
        switch (ch)
        {
            case 'i': in_file = string(optarg); cout << "Input FASTQ file: " << in_file << endl; break;
            case 'a': adapter_file_name = string(optarg); cout << "Input adapter FASTA file: " << adapter_file_name << endl; break;
            case 'd': outdir = string(optarg); cout << "Output directory: " << outdir << endl; break;
            case 'e': existence_dir = atoi(optarg); cout << "Whether output directory need to be built: " << existence_dir << endl; break;
            case 'c': HEADCUTOFF = atoi(optarg); cout << "Head-cutoff: " << HEADCUTOFF << endl; break;
        }
    }
    cout << "Waiting to be processed..." << endl;
    
    
    if (in_file.length() == 0 || adapter_file_name.length() == 0)
    {
        cout << "At least TWO parameters should be provided." << endl;
        cout << "Example: ./trim -i filename.fastq -a adapter.fa [-c HEADCUTOFF] [-d outdir [-e dir_existence]]" << endl;
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
            cout << "Example: ./trim -i filename.fastq -a adapter.fa [-c HEADCUTOFF] [-d outdir [-e dir_existence]]" << endl;
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
    
    char adapter_filename[100];
    string_convert(adapter_file_name, adapter_filename);
    ifstream adapter_file;
    adapter_file.open(adapter_filename);
    
    char main_output_name[100];
    string main_output = pre_dir+"/trimmed.fastq";
    string_convert(main_output, main_output_name);
    ofstream out;
    out.open(main_output_name);
    
    char removed_output_name[100];
    string removed_output = pre_dir+"/removed.fastq";
    string_convert(removed_output, removed_output_name);
    ofstream rm;
    rm.open(removed_output_name);
    
    char length_output_name[100];
    string length_output = pre_dir+"/polyA_length.txt";
    string_convert(length_output, length_output_name);
    ofstream length;
    length.open(length_output_name);
    
    char log_output_name[100];
    string log_output = pre_dir+"/trim.log";
    string_convert(log_output, log_output_name);
    ofstream log;
    log.open(log_output_name);
    
    
    ////////////////////////////////////////////////////////////////////////////////
    ////Receive matched fasta file as the filtered reference////////////////////////
    
    const int MAX_ADAPTER_NUMBER = 10;
    string adapter[MAX_ADAPTER_NUMBER];
    string adapter_name[MAX_ADAPTER_NUMBER];
    int adapter_length[MAX_ADAPTER_NUMBER];
    for (int i = 0; i < MAX_ADAPTER_NUMBER; i++) adapter[i] = "";
    for (int i = 0; i < MAX_ADAPTER_NUMBER; i++) adapter_name[i] = "";
    for (int i = 0; i < MAX_ADAPTER_NUMBER; i++) adapter_length[i] = 0;
    
    string line;
    
    int timer = 0;
    
    while (getline (adapter_file, line))
    {
        int lg = int(line.length());
        char marker = line[0];
        if (marker == '>')
        {
            if (lg == 1)
            {
                //cout << "Please provide name for all the adapter sequences." << endl;
                //exit(0);
            }
            else
            {
                if (line[1]==' ') adapter_name[timer] = line.substr(2, lg-2);
                else adapter_name[timer] = line.substr(1, lg-1);
                cout << "Receive adapter sequence: " << adapter_name[timer] << endl;
            }
        }
        
        else if (marker=='A' || marker=='G' || marker=='T' || marker=='C')
        {
            adapter[timer] = line;
            adapter_length[timer] = lg;
            timer++;
        }
        else continue;
    }
    cout << endl;
    
    int adapter_number = timer;
    int search_limit[timer];
    for (int i=0; i<timer; i++)
    {
        search_limit[i]=min(threshold[4],adapter_length[i]);
    }
    
    ////////////////////////////////////////////////////////////////////////////////
    //////////////////////////PROCESS FASTQ/////////////////////////////////////////
    
    string tmp1;
    string tmp2;
    string tmp3;
    string tmp4;
    string t;
    string q;
    string bait = "AGATCGGAAGAGCAC";
    int offset = 1;
    int bait_length = int(bait.length());
    bool whether_skip = false;
    bool whether_to_cutoff_terminal_base = false;
    long counter_trimmed = 0;
    //long counter_filt = 0;
    long counter_del = 0;
    int count = 0;
    unsigned long len = 0;
    int marker_point = 0;
    int pointer = 0;
    int polyA_len = 0;
    int final_len = 0;
    int lll = 150-HEADCUTOFF;
    long length_recorder[lll];
    for (int k = 0; k < lll; k ++) length_recorder[k] = 0;
    
    int cutoff = 0;
    int cutsite = 0;
    int cutsite_qualitycontrol = 0;
    
    int size = 2*(lll/threshold[6]+1)+2;
    int *record = new int[size];
    for (int i=0; i<size; i++) record[i]=0;
    //record[0]: length of s
    //record[2N-1]: length of continuous 'A'
    //record[2N]: the positon of poly'A' START site
    
    
    if (in) // file exists
    {
        ///////////////////////////////REPORT TIMER/////////////////////////////////
        
        string b = pre_dir+"/tmp2";
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
            sscanf(number_line,"%d",&lines);
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
        unsigned long compare[20];
        for (int i=0; i<20; i++) compare[i]=(lines/20)*(i+1);
        
        
        ///////////////////////////////PROCESSING///////////////////////////////////
        
        log << "ALL Sequences would be trimmed off " << HEADCUTOFF << " HEAD-base(s) in first quality control step and then processed." << endl;
        
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
                len = tmp2.length();
                cutoff = quality_control(tmp4, int(len), threshold, HEADCUTOFF);
                record[0] = cutoff - HEADCUTOFF;
                
                if (record[0] <= threshold[9])
                {
                    rm << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                    counter_del += 1;
                    log << "Sequence " << (process_timer+1)/4 << " is removed in first quality control step (too short for next-step processing)." << endl;
                }
                else
                {
                    log << "Sequence " << (process_timer+1)/4 << " is trimmed off " << len-cutoff << " TAIL-base(s) in first quality control step and ";
                    if (len-cutoff == 0) whether_to_cutoff_terminal_base = true;
                    else whether_to_cutoff_terminal_base = false;
                    
                    t = tmp2.substr(HEADCUTOFF, record[0]);
                    q = tmp4.substr(HEADCUTOFF, record[0]);
                    whether_skip = false;
                    
                    if (bait_length)
                    {
                        size_t bait_test = t.find(bait, 0);
                        if (bait_test < record[0])
                        {
                            cutsite = judge_position(t, q, int(bait_test), threshold);
                            cutsite_qualitycontrol = quality_control(q, cutsite-1, threshold);//trim off another base
                            
                            if (cutsite_qualitycontrol > 0)
                            {
                                t = t.substr(0, cutsite_qualitycontrol);
                                q = q.substr(0, cutsite_qualitycontrol);
                                out << tmp1 << "\n" << t << "\n" << tmp3 << "\n" << q << endl;
                                counter_trimmed++;
                                polyA_len = record[0]-cutsite;
                                log << "another " << polyA_len+1 << " TAIL-base(s) in adapter-directed poly-A trimming step." << endl;
                                length_recorder[polyA_len-1] += 1;
                            }
                            else
                            {
                                rm << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                                counter_del += 1;
                                polyA_len = record[0]-cutsite;
                                log << "is then removed in adapter-directed poly-A trimming step." << endl;
                                length_recorder[polyA_len-1] += 1;
                            }
                            whether_skip = true;
                        }
                    }
                    if (!whether_skip)
                    {
                        if (start_adapter_match(t, q, adapter, threshold, adapter_number))
                        {
                            rm << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                            counter_del += 1;
                            log << "is then removed in start site adapter matching step." << endl;
                        }
                        else
                        {
                            if (maximal_A_length_search(t, q, record, threshold, marker_point))
                            {
                                pointer = judge_adapter(t, q, record, adapter, threshold, adapter_number, adapter_length, marker_point, offset);
                                
                                if (pointer)
                                {
                                    cutsite_qualitycontrol = quality_control(q, record[2*pointer]-1, threshold);//trim off another base
                                    
                                    if (cutsite_qualitycontrol > 0)
                                    {
                                        t = t.substr(0, cutsite_qualitycontrol);
                                        q = q.substr(0, cutsite_qualitycontrol);
                                        out << tmp1 << "\n" << t << "\n" << tmp3 << "\n" << q << endl;
                                        counter_trimmed++;
                                        polyA_len = record[0]-record[2*pointer];
                                        log << "another " << polyA_len+1 << " TAIL-base(s) in array-traversal poly-A trimming step." << endl;
                                        length_recorder[polyA_len-1] += 1;
                                    }
                                    else
                                    {
                                        rm << tmp1 << "\n" << tmp2 << "\n" << tmp3 << "\n" << tmp4 << endl;
                                        counter_del += 1;
                                        polyA_len = record[0]-record[2*pointer];
                                        log << "is then removed in array-traversal poly-A trimming step." << endl;
                                        length_recorder[polyA_len-1] += 1;
                                    }
                                }
                                else
                                {
                                    if (whether_to_cutoff_terminal_base)
                                    {
                                        final_len = int(t.length());
                                        t = t.substr(0, final_len-1);
                                        q = q.substr(0, final_len-1);
                                        out << tmp1 << "\n" << t << "\n" << tmp3 << "\n" << q << endl;
                                        log << "trimmed off 1 TAIL-base as offset in second trimming step." << endl;
                                    }
                                    else
                                    {
                                        out << tmp1 << "\n" << t << "\n" << tmp3 << "\n" << q << endl;
                                        log << "intact in second trimming step." << endl;
                                    }
                                }
                            }
                            
                            else
                            {
                                cout << "Cannot process the FASTQ file at line " << 4*process_timer+1 << "." << endl;
                                exit(0);
                            }
                        }
                    }
                }
            }
            
            count = (count+1)%4;
            
            //print job status
            if (process_timer >= compare[report_timer])
            {
                cout << "Processed "<< 5*(report_timer+1) << "%" <<endl;
                report_timer++;
            }
            
            process_timer++;
        }
        
        for (int j = 0; j < lll; j++)
        {
            int l = j + 1;
            length << l << " : " << length_recorder[j] << endl;
        }
        
        cout << endl;
        cout << "Filter result:" << endl;
        cout << "Trimmed: " << counter_trimmed << endl;
        //cout << "Filtered: " << counter_filt << endl;
        cout << "Del: " << counter_del << endl;
        
        log << "Filter result:" << endl;
        log << "Trimmed: " << counter_trimmed << endl;
        //log << "Filtered: " << counter_filt << endl;
        log << "Del: " << counter_del << endl;
        
    }
    
    else // file doesn't exist
    {
        cout <<"No such file!" << endl;
    }
    
    ///////////////////////////////CLOSE FILES//////////////////////////////////

    in.close();
    out.close();
    log.close();
    rm.close();
    adapter_file.close();
    length.close();
    
    return 0;
}

