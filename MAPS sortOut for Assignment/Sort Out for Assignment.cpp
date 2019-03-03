/***********************************************************************
MAPS Assignment 1 - Starting Point for OMP/SIMD Assignment

SortOut - integer sorting and file output program.

The data is held in 2000 rows of 1000 numbers each. Each row is sorted independently,
currently using a simple bubble sort.

Your version should also include sorting all of the data at once, i.e. all 2,000,000 numbers,
and this can be done after the rows have been sorted.

Outputting the strings should employ SIMD to exploit hardware parallelism.

This version includes basic timing information and uses strings to create the file output.

S.Andrews / A.Oram
Revised: A.Oram Feb 2018
************************************************************************
PLEASE ADD YOUR NAME AND STUDENT NUMBER HERE:
Joshua Windsor 
24008182
************************************************************************/

#include <fstream>			//for file output
#include <iostream>			//for console output
#include <conio.h>			//for kbhit
#include "hr_time.h"		//for stopwatches
#include <stdio.h>			//for fputs
#include <omp.h>
#include <immintrin.h>

using namespace std;

//2000 * 1000
#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6			// numbers are in the range 1- 32,767, so 6 digits is enough.
#define ASCII 0x30	// add this to a single digit 0-9 to make an ASCII character

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"		// for future use

int _data [MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
int _temp [MAX_COLS];

const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2;

void getData(void);
void sortEachRow(void);
void displayCheckData(void);
void outputTimes(void);
void outputDataAsString(void);
void merge(int iArray[], int iBeginning, int iMiddle, int iEnd);
void mergesort(int iArray[], int iBeginning, int iEnd);
void ItoA_SIMD(int *num, char * numStr);

int main(void)
{
	getData();

	cout << "Sorting data...";
	s1.startTimer();
	sortEachRow();
	s1.stopTimer();

	displayCheckData();

	cout << "\n\nOutputting data to " << SortedRows << "...";
	s2.startTimer();
	outputDataAsString();
	s2.stopTimer();

	outputTimes();

	while (!_kbhit());  //to hold console
}

//*********************************************************************************
void getData()		// Generate the same sequence of 'random' numbers.
{
	srand(123); //random number seed PLEASE DON'T CHANGE!
	for (int i = 0; i<MAX_ROWS; i++)
	{
		for (int j = 0; j < MAX_COLS; j++)
		{
			_data[i][j] = rand(); //RAND_MAX = 32767
		}
	}
}


//*********************************************************************************
void sortEachRow()
{
	for(int i=0; i<MAX_ROWS; i++)
	{
		mergesort(_data[i], 0, MAX_COLS - 1);
	}
}


//*********************************************************************************
void displayCheckData()
{
	cout << "\n\ndata[0][0]                   = " << _data[0][0] << "\t" << (_data[0][0] == checkBeg ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS/2][MAX_COLS/2] = " << _data[MAX_ROWS / 2][MAX_COLS / 2] << "\t" << (_data[MAX_ROWS / 2][MAX_COLS / 2] == checkMid ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS-1][MAX_COLS-1] = " << _data[MAX_ROWS - 1][MAX_COLS - 1] << "\t" << (_data[MAX_ROWS - 1][MAX_COLS - 1] == checkEnd ? "OK" : "BAD");
}




//*********************************************************************************
void outputTimes()
{
	ofstream os;
	os.open("StoredTimes.txt", ios::app);
	os << "\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	os << "\nTime for outputting to file (s) : " << s2.getElapsedTime();
	os << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() << "\n\n-----------------------------";
	os.close();


	cout << "\n\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	cout << "\nTime for outputting to file (s) : " << s2.getElapsedTime();
	cout << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() << "\n\n\nPress a key to terminate.";
}


//*********************************************************************************
//Builds a sorted number list as a long string then outputs the whole thing in one big fputs!

void outputDataAsString()
{
	char numString[MAX_CHARS * 8];
	string odata;

	for (int i = 0; i<MAX_ROWS; i++){
		for (int j = 0; j<MAX_COLS; j+=8){
			ItoA_SIMD(&_data[i][j], numString);
			odata += numString;
			odata += "\t";
		}
		odata += "\n";
	}

	FILE * sodata;
	fopen_s(&sodata, SortedRows, "w");
	fputs(odata.c_str(), sodata);
	fclose(sodata);
}


//Ref: https://www.comrevo.com/2016/02/openmp-program-for-merge-sort.html
//TODO: SIMD

inline void merge(int iArray[], int iBeginning, int iMiddle, int iEnd)
{
	int a,
		b = iBeginning,
		c,
		d = iMiddle + 1;

	for (a = iBeginning; b <= iMiddle && d <= iEnd; a++)
	{
		if (iArray[b] <= iArray[d])
		{
			_temp[a] = iArray[b];
			b++;
		}
		else
		{
			_temp[a] = iArray[d];
			d++;
		}
	}
	if (b > iMiddle)
	{
		for (c = d; c <= iEnd; c++)
		{
			_temp[a] = iArray[c];
			a++;
		}
	}
	else
	{
		for (c = b; c <= iMiddle; c++)
		{
			_temp[a] = iArray[c];
			a++;
		}
	}
	for (c = iBeginning; c <= iEnd; c++)
	{
		iArray[c] = _temp[c];
	}
}


inline void mergesort(int iArray[], int iBeginning, int iEnd)
{
	if (iBeginning < iEnd)
	{
		int middle = (iBeginning + iEnd) / 2;

		#pragma omp parallel sections num_threads(2)
		{
			#pragma omp section
			{
				mergesort(iArray, iBeginning, middle);
			}

			#pragma omp section
			{
				mergesort(iArray, middle + 1, iEnd);
			}
		}
		merge(iArray, iBeginning, middle, iEnd);
	}

}

//Ref: https://github.com/SeanShortreed/ParallelProgramming/blob/22d4a97ccc6494ea8234f6f7a5005885b4e0df81/shortreed_sean_SIMD.cpp
//TODO: Make unique, rename stuff and fix up


inline void ItoA_SIMD(int *num, char * numStr)
{
	_declspec(align(16)) short DataToSort[8] = { num[0], num[1],  num[2],  num[3], num[4], num[5], num[6], num[7] };

	//Simd vars
	__m128i* OrigData = (__m128i*)DataToSort;
	__m128i SimdData;
	__m128i opNumber;
	__m128i DataToUse = *OrigData;
	__m128i SavedData;

	//start for loop
	for (int i = 0; i < 5; i++)
	{
		SimdData = DataToUse;
		SavedData = SimdData;
		//setup opNumber with magic number
		opNumber = _mm_set1_epi16(0x6667);
		//times by magic number
		SimdData = _mm_mulhi_epi16(opNumber, SimdData);
		//shift right
		SimdData = _mm_srai_epi16(SimdData, 2);
		//Save these numbers for next round
		DataToUse = SimdData;
		//setup opNumber with 10
		opNumber = _mm_set1_epi16(10);
		//times by 10
		SimdData = _mm_mullo_epi16(opNumber, SimdData);
		//take new num from orig num
		SimdData = _mm_sub_epi16(SavedData, SimdData);
		//SimdData now contains single digits
		//Convert SimdData to char numbers
		//Setup opNumber with char conversion
		opNumber = _mm_set1_epi16(ASCII);
		//do the addition
		SimdData = _mm_add_epi16(opNumber, SimdData);
		//SimdData now contains correct chars
		//sort letters into correct arrays
		numStr[4 - i] = SimdData.m128i_i16[0];
		numStr[6 + 4 - i] = SimdData.m128i_i16[1];
		numStr[12 + 4 - i] = SimdData.m128i_i16[2];
		numStr[18 + 4 - i] = SimdData.m128i_i16[3];
		numStr[24 + 4 - i] = SimdData.m128i_i16[4];
		numStr[30 + 4 - i] = SimdData.m128i_i16[5];
		numStr[36 + 4 - i] = SimdData.m128i_i16[6];
		numStr[42 + 4 - i] = SimdData.m128i_i16[7];
	}

	numStr[5] = '\t';
	numStr[6 + 5] = '\t';
	numStr[12 + 5] = '\t';
	numStr[18 + 5] = '\t';
	numStr[24 + 5] = '\t';
	numStr[30 + 5] = '\t';
	numStr[36 + 5] = '\t';
	numStr[42 + 5] = '\0';
}
//***********************************************************