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

#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6			// numbers are in the range 1- 32,767, so 6 digits is enough.

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"

const __m256i ascii_magic_1_16 = _mm256_set1_epi16(0x6667);
const __m256i ascii_magic_2_16 = _mm256_set1_epi16(0x30);
const __m256i ten_16 = _mm256_set1_epi16(10);
const __m128i ascii_magic_1_8 = _mm_set1_epi16(0x6667);
const __m128i ascii_magic_2_8 = _mm_set1_epi16(0x30);
const __m128i ten_8 = _mm_set1_epi16(10);


int _data [MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
int _overflow = MAX_COLS - 8;		// overflow value for 1000/16
int _allData [MAX_ROWS][MAX_COLS];	// flat version of _data for outputting

const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2, s3, s4;

__forceinline void getData(void);
__forceinline void sortEachRow(void);
__forceinline void displayCheckData(void);
__forceinline void sortAll(void);
__forceinline void outputTimes(void);
__forceinline void outputDataAsString(int iOutputData[][MAX_COLS], string iFilename);
__forceinline void mergeSort(int iArray[], int iBeginning, int iEnd);
__forceinline void SIMDitoa16(int * iArray, char * oNumString);
__forceinline void SIMDitoa8 (int * iArray, char * oNumString);
__forceinline void radixSort(int iArray[]);
__forceinline void radixSortAll(int iArray[][MAX_COLS]);


//*********************************************************************************
int main(void)
{
	//remove previous files due to lab issues
	remove(SortedRows);
	remove(SortedAll);

	getData();

	//speeds up processing 
	omp_set_dynamic(1);
	//seems to slow down when turned on (from 0.3 to 0.45)
	omp_set_nested(0);

	cout << "Threads: " << omp_get_thread_num();
	cout << "\nDynamic: " << omp_get_dynamic();
	cout << "\nNested: " << omp_get_nested();
	cout << "\n\nSorting data...";

	s1.startTimer();
	sortEachRow();
	s1.stopTimer();

	displayCheckData();

	cout << "\n\nOutputting data to " << SortedRows << "...";
	s2.startTimer();
	outputDataAsString(_data, SortedRows);
	s2.stopTimer();

	cout << "\n\nSorting all...";
	s3.startTimer();
	sortAll();
	s3.stopTimer();

	cout << "\n\nOutputting data to " << SortedAll << "...";
	s4.startTimer();
	outputDataAsString(_allData, SortedAll);
	s4.stopTimer();

	outputTimes();

	while (!_kbhit());  //to hold console
}

//*********************************************************************************
__forceinline void getData()		// Generate the same sequence of 'random' numbers.
{
	srand(123); //random number seed PLEASE DON'T CHANGE!
	for (int i = 0; i<MAX_ROWS; ++i)
	{
		for (int j = 0; j < MAX_COLS; ++j)
		{
			_data[i][j] = rand(); //RAND_MAX = 32767
		}
	}
}

//*********************************************************************************
__forceinline void sortEachRow()
{
	//basic parallel for to speed up individual rows
	#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; ++i)
	{
		radixSort(_data[i]);
	}

}

//*********************************************************************************
__forceinline void sortAll()
{
	//runs a full radix without omp as it applies to the entire array
	radixSortAll(_data);
}

//*********************************************************************************
__forceinline void displayCheckData()
{
	cout << "\n\ndata[0][0]                   = " << _data[0][0] << "\t" << (_data[0][0] == checkBeg ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS/2][MAX_COLS/2] = " << _data[MAX_ROWS / 2][MAX_COLS / 2] << "\t" << (_data[MAX_ROWS / 2][MAX_COLS / 2] == checkMid ? "OK" : "BAD");
	cout << "\ndata[MAX_ROWS-1][MAX_COLS-1] = " << _data[MAX_ROWS - 1][MAX_COLS - 1] << "\t" << (_data[MAX_ROWS - 1][MAX_COLS - 1] == checkEnd ? "OK" : "BAD");
}

//*********************************************************************************
__forceinline void outputTimes()
{
	//output times to a file
	ofstream os;
	os.open("StoredTimes.txt", ios::app);
	os << "\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	os << "\nTime for outputting to file (s) : " << s2.getElapsedTime();
	os << "\nTime for sorting everything (s) : " << s3.getElapsedTime();
	os << "\nTime for outputting to file (s) : " << s4.getElapsedTime();
	os << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() + s3.getElapsedTime() + s4.getElapsedTime() << "\n\n-----------------------------";
	os.close();

	//output times to console
	cout << "\n\nTime for sorting all rows   (s) : " << s1.getElapsedTime();
	cout << "\nTime for outputting to file (s) : " << s2.getElapsedTime();
	cout << "\nTime for sorting everything (s) : " << s3.getElapsedTime();
	cout << "\nTime for outputting to file (s) : " << s4.getElapsedTime();
	cout << "\nCombined time               (s) : " << s1.getElapsedTime() + s2.getElapsedTime() + s3.getElapsedTime() + s4.getElapsedTime() << "\n\n\nPress a key to terminate.";
}

//*********************************************************************************
//Builds a sorted number list as a long string then outputs the whole thing in one big fputs!
__forceinline void outputDataAsString(int iOutputData[][MAX_COLS], string iFilename)
{
	char numString[MAX_CHARS * 16];
	string odata;

	//A parallel for could be implemeneted such as this to speed up 
	//#pragma omp parallel for private(iOutputData, numString)
	//Or even split into 8 sections similar to the radixSortAll 
	//and combined at the end
	for (int i = 0; i<MAX_ROWS; ++i){
		for (int j = 0; j<MAX_COLS; j+=16){
			//determines wether to do a 8 or 16 simd operation depending on how many left
			if (j != _overflow)
			{
				SIMDitoa16(&iOutputData[i][j], numString);
			}
			else
			{
				SIMDitoa8(&iOutputData[i][j], numString);
			}
			//loops through all 16 nums
			for (int k = 0; k < 16; ++k)
			{
				//gets the start index of the number
				int simdIdx = k * MAX_CHARS;
				//loops through each digit (only 4 to show 0 values)
				for (int m = 0; m < 4; ++m)
				{
					//if it's a 0
					if (numString[simdIdx + m] == '0')
					{
						//replace with a space
						numString[simdIdx + m] = ' ';
					}
					//or break as no more 0s to remove
					else
					{
						break;
					}
				}
			}
			//Store output data in output string
			odata += numString;
			odata += "\t";
		}
		odata += "\n";
	}

	FILE * sodata;
	fopen_s(&sodata, iFilename.c_str(), "w");
	fputs(odata.c_str(), sodata);
	fclose(sodata);
}

#pragma region RemovedMergeSort
//*********************************************************************************
//Ref: https://www.comrevo.com/2016/02/openmp-program-for-merge-sort.html
//Removed as it seems slower than a radix sort

__forceinline void merge(int iArray[], int iTemp[], int iBeginning, int iMiddle, int iEnd)
{
	int a,
		b = iBeginning,
		c,
		d = iMiddle + 1;

	for (a = iBeginning; b <= iMiddle && d <= iEnd; a++)
	{
		if (iArray[b] <= iArray[d])
		{
			iTemp[a] = iArray[b];
			b++;
		}
		else
		{
			iTemp[a] = iArray[d];
			d++;
		}
	}
	if (b > iMiddle)
	{
		for (c = d; c <= iEnd; c++)
		{
			iTemp[a] = iArray[c];
			a++;
		}
	}
	else
	{
		for (c = b; c <= iMiddle; c++)
		{
			iTemp[a] = iArray[c];
			a++;
		}
	}
	for (c = iBeginning; c <= iEnd; c++)
	{
		iArray[c] = iTemp[c];
	}
}

//*********************************************************************************
__forceinline void mergesort(int iArray[], int iBeginning, int iEnd)
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

		int temp[MAX_COLS];
		merge(iArray, temp, iBeginning, middle, iEnd);
	}
}

#pragma endregion RemovedMergeSort

#pragma region SIMDitoa16&8

//*********************************************************************************
//Uses 256 bit SIMD registers to convert 16 numbers to ascii at a time
__forceinline void SIMDitoa16(int *iArray, char * oNumString)
{
	//setup a memory aligned array with the 16 input ints
	_declspec(align(32)) short intArray[16] =
	{
		iArray[0], iArray[1],  iArray[2],   iArray[3],  iArray[4],  iArray[5],  iArray[6],  iArray[7],
		iArray[8], iArray[9],  iArray[10],  iArray[11], iArray[12], iArray[13], iArray[14], iArray[15] 
	};

	//declare some simd variables to use in the loop
	__m256i startingArray = _mm256_load_si256((__m256i*)intArray);
	__m256i outputData;
	__m256i inputData = startingArray;
	__m256i storedData;

	for (int i = 0; i < 5; ++i)
	{
		//starts the 2 vars which will be iterated upon
		outputData = inputData;
		storedData = outputData;

		//multiply by the first magic number
		outputData = _mm256_mulhi_epu16(ascii_magic_1_16, outputData);

		//shift the data right
		outputData = _mm256_srai_epi16(outputData, 2);

		//store these numbers for future iterations
		inputData = outputData;

		//multiply by 10
		outputData = _mm256_mullo_epi16(ten_16, outputData);

		//subtract the generated value from the original value
		outputData = _mm256_sub_epi16(storedData, outputData);

		//add the second magic number to complete single
		//digit to char conversion
		outputData = _mm256_add_epi16(ascii_magic_2_16, outputData);

		//sort letters into correct 16 output arrays backwards
		oNumString[4 - i] = outputData.m256i_i16[0];
		oNumString[10 - i] = outputData.m256i_i16[1];
		oNumString[16 - i] = outputData.m256i_i16[2];
		oNumString[22 - i] = outputData.m256i_i16[3];
		oNumString[28 - i] = outputData.m256i_i16[4];
		oNumString[34 - i] = outputData.m256i_i16[5];
		oNumString[40 - i] = outputData.m256i_i16[6];
		oNumString[46 - i] = outputData.m256i_i16[7];
		oNumString[52 - i] = outputData.m256i_i16[8];
		oNumString[58 - i] = outputData.m256i_i16[9];
		oNumString[64 - i] = outputData.m256i_i16[10];
		oNumString[70 - i] = outputData.m256i_i16[11];
		oNumString[76 - i] = outputData.m256i_i16[12];
		oNumString[82 - i] = outputData.m256i_i16[13];
		oNumString[88 - i] = outputData.m256i_i16[14];
		oNumString[94 - i] = outputData.m256i_i16[15];
	}

	//append the tabs for neat output
	oNumString[5] = '\t';
	oNumString[11] = '\t';
	oNumString[17] = '\t';
	oNumString[23] = '\t';
	oNumString[29] = '\t';
	oNumString[35] = '\t';
	oNumString[41] = '\t';
	oNumString[47] = '\t';
	oNumString[53] = '\t';
	oNumString[59] = '\t';
	oNumString[65] = '\t';
	oNumString[71] = '\t';
	oNumString[77] = '\t';
	oNumString[83] = '\t';
	oNumString[89] = '\t';
	//end with null terminator
	oNumString[95] = '\0';
}

//*********************************************************************************
//Uses 128 bit SIMD registers to convert 8 numbers to ascii at a time (Used to do overflow of 8 digits when 1000/16)
__forceinline void SIMDitoa8(int *iArray, char * oNumString)
{
	//setup a memory aligned array with the 8 input ints
	_declspec(align(16)) short intArray[8] = {
		iArray[0], iArray[1],  iArray[2],  iArray[3], iArray[4], iArray[5], iArray[6], iArray[7] 
	};
	//declare some simd variables to use in the loop
	__m128i startingArray = _mm_load_si128((__m128i*)intArray);
	__m128i outputData;
	__m128i inputData = startingArray;
	__m128i storedData;

	for (int i = 0; i < 5; ++i)
	{
		//starts the 2 vars which will be iterated upon
		outputData = inputData;
		storedData = outputData;

		//multiply by the first magic number
		outputData = _mm_mulhi_epi16(ascii_magic_1_8, outputData);

		//shift the data right
		outputData = _mm_srai_epi16(outputData, 2);

		//store these numbers for future iterations
		inputData = outputData;

		//multiply by 10
		outputData = _mm_mullo_epi16(ten_8, outputData);

		//subtract the generated value from the original value
		outputData = _mm_sub_epi16(storedData, outputData);

		//add the second magic number to complete single
		//digit to char conversion
		outputData = _mm_add_epi16(ascii_magic_2_8, outputData);

		//sort letters into correct 8 output arrays
		oNumString[4 - i] = outputData.m128i_i16[0];
		oNumString[10 - i] = outputData.m128i_i16[1];
		oNumString[16 - i] = outputData.m128i_i16[2];
		oNumString[22 - i] = outputData.m128i_i16[3];
		oNumString[28 - i] = outputData.m128i_i16[4];
		oNumString[34 - i] = outputData.m128i_i16[5];
		oNumString[40 - i] = outputData.m128i_i16[6];
		oNumString[46 - i] = outputData.m128i_i16[7];
	}

	//append the tabs for neat output
	oNumString[5] = '\t';
	oNumString[11] = '\t';
	oNumString[17] = '\t';
	oNumString[23] = '\t';
	oNumString[29] = '\t';
	oNumString[35] = '\t';
	oNumString[41] = '\t';
	//end with null terminator
	oNumString[47] = '\0';
}

#pragma endregion SIMDitoas

#pragma region RowSorting
//*********************************************************************************
//Ref: https://en.wikipedia.org/wiki/Radix_sort
//based on Java implementation on Wikipedia
__forceinline int getMax(int iArray[])
{
	// determines the highest number in the array
	int max = iArray[0];

	for (int i = 1; i < MAX_ROWS; ++i)
	{
		if (iArray[i] > max)
		{
			max = iArray[i];
		}
	}
	return max;
}

//*********************************************************************************
// A function to do counting sort of iArray[] according to 
// the digit represented by exp. (eg. 300 is represented by 100)
__forceinline void countSort(int iArray[], int iExp)
{
	int output[MAX_COLS], 
		i, 
		count[10] = { 0 };

	// Store count of occurrences in count[] 
	for (i = 0; i < MAX_COLS; ++i)
	{
		count[(iArray[i] / iExp) % 10]++;
	}

	// Change count[i] so that count[i] now contains 
	// actual position of this digit in output[] 
	for (i = 1; i < 10; ++i)
		count[i] += count[i - 1];

	// Build the output array 
	for (i = MAX_COLS - 1; i >= 0; --i)
	{
		output[count[(iArray[i] / iExp) % 10] - 1] = iArray[i];
		count[(iArray[i] / iExp) % 10]--;
	}

	// Copy the output array to iArray[], so that iArray[] now 
	// contains sorted numbers according to curent digit 
	// Apply omp to speed up array transfer
	#pragma omp parallel for
	for (i = 0; i < MAX_COLS; ++i)
	{
		iArray[i] = output[i];
	}
}

//*********************************************************************************
// Radix Sort
__forceinline void radixSort(int iArray[])
{
	// Find the maximum number to know number of digits
	int m = getMax(iArray);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
	{
		countSort(iArray, exp);
	}
}

#pragma endregion RowSorting

#pragma region AllSorting
//*********************************************************************************
// determines the highest number in the 2D array
__forceinline int getMaxAll(int iArray[][MAX_COLS])
{
	int max = iArray[0][0];

	for (int i = 0; i < MAX_ROWS; ++i)
	{
		for (int j = 0; j < MAX_COLS; ++j)
		{
			if (iArray[i][j] > max)
			{
				max = iArray[i][j];
			}
		}
	}
	return max;
}

//*********************************************************************************
// Does the main sort of the radix sort
__forceinline void countSortAll(int iArray[][MAX_COLS], int iExp)
{
	int output[10] = { 0 };

	//create 8 sets of buckets
	int bucket1[10] = { 0 },
		bucket2[10] = { 0 },
		bucket3[10] = { 0 },
		bucket4[10] = { 0 },
		bucket5[10] = { 0 },
		bucket6[10] = { 0 },
		bucket7[10] = { 0 },
		bucket8[10] = { 0 };
	
	
	int sectionSize = MAX_ROWS / 8;
	// Create 8 sections, and fill the 8 sets of buckets for 8threads
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i = 0; i < sectionSize; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					//store the occureces in the buckets
					bucket1[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize; i < sectionSize * 2; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket2[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 2; i < sectionSize * 3; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket3[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 3; i < sectionSize * 4; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket4[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 4; i < sectionSize * 5; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket5[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 5; i < sectionSize * 6; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket6[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 6; i < sectionSize * 7; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket7[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 7; i < MAX_ROWS; ++i)
			{
				for (int j = 0; j < MAX_COLS; ++j)
				{
					bucket8[(iArray[i][j] / iExp) % 10]++;
				}
			}
		}
	}

	//add the buckets together
	//OpenMP showed slowdowns when applied here
	for (int i = 0; i < 10; ++i)
	{
		output[i] = bucket1[i] + bucket2[i] + bucket3[i] + bucket4[i] + bucket5[i] + bucket6[i] + bucket7[i] + bucket8[i];
	}

	// Change output[i] so that output[i] now contains actual
	// position of this digit in output[]
	// OpenMP cannot be applied as array is being modified
	for (int i = 1; i < 10; ++i)
	{
		output[i] += output[i - 1];
	}

	// Build the output array
	// OpenMP cannot be applied as array is being modified
	for (int i = MAX_ROWS - 1; i >= 0; --i)
	{
		for (int j = MAX_COLS - 1; j >= 0; --j)
		{
			int index = output[(iArray[i][j] / iExp) % 10] - 1;

			_allData[index / MAX_COLS][index % MAX_COLS] = iArray[i][j];
			output[(iArray[i][j] / iExp) % 10]--;
		}
	}

	// Copy the output array to iArray[], so that iArray[] now 
	// contains sorted numbers according to curent digit 
	#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; ++i)
	{
		for (int j = 0; j < MAX_COLS; ++j)
		{
			iArray[i][j] = _allData[i][j];
		}
	}
}

//*********************************************************************************
// Radix Sort for all numbers
__forceinline void radixSortAll(int iArray[][MAX_COLS])
{
	// Find the maximum number to know number of digits
	int m = getMaxAll(iArray);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
	{
		countSortAll(iArray, exp);
	}
}
#pragma endregion AllSorting
