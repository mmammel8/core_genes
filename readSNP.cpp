/*
readSNP.cpp
read alignment files and produce SNP fasta file
v1.1 MKM 4/5/2018
*/
//#include "stdafx.h"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <cstring>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include <dirent.h>
//#include <windows.h>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
using namespace std;

#define WIN 0

const int MAXREP4SNP = 804;
string dir = "/mnt/dmb/Mark_backup/SLSNP/aln1/";
string file1 = "AE006468.2.txt";
string file2 = "labeltableSL2.txt";
string fileOut1 = "SL2.fasta";
string fileOut2 = "tallySL2.csv";
string fileOut3 = "SL2_annot.txt";

const int NTALLYCHAR = 18;

#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
	#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
	#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif
typedef unsigned long long ui64;

typedef struct
{
	string accession;
	string label;
	string group;
	string seq;
	int tally[NTALLYCHAR];

} Accession;

typedef struct
{
	string locus;
	int start, end, ostart;
	string strand, gene, function;

} Locus;

typedef struct
{
	int pos, chromosome_pos;
	string snps;

} SNP;

struct by_pos
{
	bool operator()(SNP const &a, SNP const &b)
	{
		return a.chromosome_pos < b.chromosome_pos;
	}
};

//globals
vector <Accession> access_list;
vector <Locus> locus_list;
unordered_map<string, int> map_acc;
unordered_map<string, int> map_locus;
unordered_set<string> names;

#if WIN == 1
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	ui64 tmpres = 0;
	static int tzflag;

	if (NULL != tv)
	{
		GetSystemTimeAsFileTime(&ft);

		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;

		/*converting file time to unix epoch*/
		tmpres -= DELTA_EPOCH_IN_MICROSECS;
		tmpres /= 10;  /*convert into microseconds*/
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (long)(tmpres % 1000000UL);
	}

	return 0;
}
#endif

bool get_directory(string path)
{
	string locus;
#if WIN == 1
	path += "*.*";
	WIN32_FIND_DATAA fd;
	HANDLE hFind = ::FindFirstFileA(LPCSTR(path.c_str()), &fd);
	if (hFind != INVALID_HANDLE_VALUE) 
	{
		do 
		{
			//if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			locus = fd.cFileName;
			size_t found = locus.find(".txt");
			if (found != string::npos)
			{
				//cout << fd.cFileName << endl;
				names.insert(locus.substr(0,found));
			}
		} while (::FindNextFileA(hFind, &fd));
		::FindClose(hFind);
	}

#else
	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(path.c_str())) != NULL) 
	{
		/* print all the files and directories within directory */
		while ((ent = readdir(dir)) != NULL) 
		{
			printf("%s\n", ent->d_name);
			locus = ent->d_name;
			size_t found = locus.find(".txt");
			if (found != string::npos)
			{
				names.insert(locus.substr(0,found));
			}
		}
		closedir(dir);
	}
	else 
	{
		/* could not open directory */
		perror("");
		return false;
	}
#endif
	return true;
}

bool load_data2()
{ //load labeltable
	ifstream fin;
	string acc, lab, dbs, group;
	int perc, i;
	//not currently using database, percentage, group
	Accession A;

	int num_acc = 0;
	A.seq = "";
	for (i = 0; i < NTALLYCHAR; i++)
		A.tally[i] = 0;
	fin.open(file2);
	while (fin >> acc >> lab >> dbs >> perc >> group)
	{
		A.accession = acc;
		A.label = lab;
		A.group = group;
		access_list.push_back(A);
		map_acc[acc] = num_acc++;
		//cout << acc << " " << num_acc << endl;
	}
	fin.close();
	cout << "accessions: " << num_acc << endl;

	return true;
}

bool load_data1()
{ //load annotation table
	ifstream fin;
	string locus, strand, gene, function;
	int start, end, length;
	int prev_end = 0;
	Locus L;

	int num_loc = 0;
	fin.open(file1);
	while (fin >> locus >> start >> end >> strand >> length >> gene)
	{
		getline(fin, function);
		function.erase(0, 1); //remove tab at start
		L.ostart = start;
		if (start <= prev_end) //no overlap
			start = prev_end + 1;
		prev_end = end;
		L.locus = locus;
		L.start = start;
		L.end = end;
		L.strand = strand;
		L.gene = gene;
		L.function = function;
		if (names.find(locus) != names.end())
		{
			locus_list.push_back(L);
			map_locus[locus] = num_loc++;
			//cout << locus << " " << num_loc << endl;
		}
	}
	fin.close();
	cout << "genes: " << num_loc << endl;

	return true;
}

bool load_data3(vector<string> &seq_rep, string fname)
{ //load alignment of accessions for this gene fname
	ifstream fin;
	string strain, sequence, seq2;
	string::iterator it1;
	unordered_map<string, int>::iterator it2;
	int idx;

	fin.open(fname);
	while (getline(fin, strain))
	{
		strain.erase(remove_if(strain.begin(), strain.end(), ::isspace), strain.end());
		strain.erase(0, 1); //remove  '>'
		getline(fin, sequence);
		it2 = map_acc.find(strain);
		if (it2 != map_acc.end())
		{
			idx = it2->second; //value
			sequence.erase(remove_if(sequence.begin(), sequence.end(), ::isspace), sequence.end());
			transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
			seq_rep[idx] = sequence;
			//cout << "add" << idx << endl;
		}
	}
	fin.close();

	return true;
}

unsigned int code(char b1)
{ //return binary code A=1, C=2, G=4, T=8
	unsigned int v1, vals[26] = { 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 8, 7, 9, 0, 10, 0 };

	if (b1 >= 'A' && b1 <= 'Z')
		v1 = vals[b1 - 'A'];
	else if (b1 == '?' || b1 == '-')
		v1 = 15; //not a snp
	else
		v1 = 0;

	return v1;
}

unsigned int tally_char(char b1)
{
	//return index for tally of chars in each accession
	//A	C	G	T	N	X	R	Y	S	W	K	M	B	D	H	V -
	unsigned int v1, vals[26] = { 1, 13, 2, 14, 0, 0, 3, 15, 0, 0, 11, 0, 12, 5, 0, 0, 0, 7, 9, 4, 4, 16, 10, 6, 8, 0 };

	if (b1 >= 'A' && b1 <= 'Z')
		v1 = vals[b1 - 'A'];
	else if (b1 == '?')
		v1 = 5; //N
	else if (b1 == '-')
		v1 = 17; //N
	else
		v1 = 0;

	return v1;
}

void find_snp(vector<SNP> &snp_list, vector<string> &seq_rep, Locus *locus)
{
	//find snps in alignment seq_rep for locus. 
	//bbonepos is cumulative gene length
	int seq_len = seq_rep[0].length();
	int acc_len = access_list.size();
	int pos, acc, mm, mmacc, pos2, chromosome_pos;
	unsigned int v1, v2;
	char b1, b2;
	string snps;
	SNP S;

	for (pos = 0; pos < seq_len; pos++)
	{
		S.snps = "";
		b1 = seq_rep[0].at(pos);
		mm = 0;
		v1 = code(b1);
		for (acc = 0; acc < acc_len; acc++)
		{
			if (pos < seq_rep[acc].length())
			{
				b2 = seq_rep[acc].at(pos);
				v2 = code(b2);
				if ((v1 & v2) == 0 && acc <= MAXREP4SNP)
				{
					mm++;
					mmacc = acc;
				}
				S.snps += b2;
			}
			else
				S.snps += '-';
		} //acc
		if (mm == 1)
		{
			//look for mismatch due to error from repeat / gap
			int st = pos - 2;
			int en = pos + 2;
			if (st < 0) { st = 0; }
			if (en >= seq_len) { en = seq_len - 1; }
			for (pos2 = st; pos2 <= en; pos2++)
			{
				b1 = seq_rep[0].at(pos2);
				if (b1 == '-') { mm = 0; }
				b2 = seq_rep[mmacc].at(pos2);
				if (b2 == '-') { mm = 0; }
			}
		}
		if (mm > 0)
		{
			if (locus->strand == "+")
			{
				chromosome_pos = pos + locus->start;
			}
			else
			{
				chromosome_pos = locus->end - pos;
			}
			S.pos = pos;
			S.chromosome_pos = chromosome_pos;
			snp_list.push_back(S);
		}

	} //pos
}

int main()
{
	double time0, time1;
	struct timeval tv;
	int idx1, bbonepos, acc;
	int concat_pos, locus_pos, chromosome_pos;
	vector<Locus>::iterator itl;
	vector<Accession>::iterator ita;
	vector<SNP>::iterator its;
	vector<SNP> snp_list;
	string mst1, mst2;
	ofstream snpfile, tallyfile, annotfile;
	char b1;
	unsigned int v1;

	gettimeofday(&tv, NULL);
	time0 = tv.tv_sec + 1e-6 * tv.tv_usec;

	get_directory(dir);
	load_data1();
	load_data2();
	vector<string> seq_rep(access_list.size());
	gettimeofday(&tv, NULL);
	time1 = tv.tv_sec + 1e-6 * tv.tv_usec;
	printf("load data time: %f\n", time1-time0);
	time0 = time1;

	bbonepos = 0;
	annotfile.open(fileOut3);
	for (itl = locus_list.begin(); itl != locus_list.end(); ++itl)
	{
		cout << itl->locus << endl;
		string fname = dir + itl->locus + ".txt";
		snp_list.clear();
		for (idx1 = 0; idx1 < access_list.size(); idx1++)
			seq_rep[idx1] = "";
		load_data3(seq_rep, fname);
		find_snp(snp_list, seq_rep, &*itl);
		sort(snp_list.begin(), snp_list.end(), by_pos());
		for (its = snp_list.begin(); its != snp_list.end(); ++its)
		{
			for (acc = 0; acc < access_list.size(); acc++)
			{
				b1 = its->snps.at(acc);
				v1 = tally_char(b1);
				access_list[acc].tally[v1]++;
				access_list[acc].seq += b1;
			}
			concat_pos = its->pos + bbonepos + 1;
			if (itl->strand == "+")
			{
				chromosome_pos = its->pos + itl->start;
				locus_pos = its->pos + 1 + itl->start - itl->ostart; //adjust for trim
			}
			else
			{
				chromosome_pos = itl->end - its->pos;
				locus_pos = its->pos + 1; //no adjustment needed
			}
			annotfile << itl->locus << '\t' << chromosome_pos << '\t' << concat_pos << '\t' << locus_pos << '\t' << itl->gene << '\t' << itl->function << endl;
		}
		bbonepos += itl->end - itl->start + 1;
	}
	gettimeofday(&tv, NULL);
	time1 = tv.tv_sec + 1e-6 * tv.tv_usec;
	printf("calculation time: %f\n", time1 - time0);
	annotfile.close();

	snpfile.open(fileOut1);
	tallyfile.open(fileOut2);
	tallyfile << "strain,0,A,C,G,T,N,X,R,Y,S,W,K,M,B,D,H,V,-" << endl;
	for (ita = access_list.begin(); ita != access_list.end(); ++ita)
	{
		snpfile << ">" << ita->label << endl << ita->seq << endl;
		tallyfile << ita->label;
		for (idx1 = 0; idx1 < NTALLYCHAR; idx1++)
			tallyfile << ',' << ita->tally[idx1];
		tallyfile << endl;
	}
	snpfile.close();
	tallyfile.close();

	return 0;
}

