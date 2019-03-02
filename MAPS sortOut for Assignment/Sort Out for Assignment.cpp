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


************************************************************************/

#include <fstream>			//for file output
#include <iostream>			//for console output
#include <conio.h>			//for kbhit
#include "hr_time.h"		//for stopwatches
#include <stdio.h>			//for fputs

using namespace std;

#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6			// numbers are in the range 1- 32,767, so 6 digits is enough.

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"		// for future use


int _data [MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2;

void getData(void);
void sortEachRow(void);
void displayCheckData(void);
void outputDataAsString(void);
void outputTimes(void);
void outputDataAsString(void);


int main(void)
{
	getData();

	s1.startTimer();
	sortEachRow();
	s1.stopTimer();

	displayCheckData();

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
		for (int j = 0; j<MAX_COLS; j++)
			_data[i][j] = rand(); //RAND_MAX = 32767
}


//*********************************************************************************
void sortEachRow()
{
	cout << "Sorting data...";
	for(int i=0; i<MAX_ROWS; i++)
	{	
		//Use a bubble sort on a row 
		for(int n=MAX_COLS-1; n>=0; n--)
		{   for(int j=0; j<n; j++)
			{
                if(_data[i][j] > _data[i][j+1])
				{
                    int temp = _data[i][j];
                    _data[i][j] = _data[i][j+1];
                    _data[i][j+1] = temp;
                }
            }
		}
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
	cout << "\n\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	cout << "\nTime for outputting to file (s) : " << s2.getElapsedTime();
	cout << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() << "\n\n\nPress a key to terminate.";
}


//*********************************************************************************
//Builds a sorted number list as a long string then outputs the whole thing in one big fputs!

void outputDataAsString()
{
	char numString[MAX_CHARS];
	string odata;
	cout << "\n\nOutputting data to " << SortedRows << "...";

	for (int i = 0; i<MAX_ROWS; i++){
		for (int j = 0; j<MAX_COLS; j++){
			_itoa_s<6>(_data[i][j], numString, 10);
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
//*********************************************************************************
