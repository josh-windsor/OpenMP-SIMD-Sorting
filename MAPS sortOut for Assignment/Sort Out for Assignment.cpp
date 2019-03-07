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

#define ASCII_MAGIC_1 0x6667	//ASCII magic numbers
#define ASCII_MAGIC_2 0x30

#define SortedRows "SortedRows.txt"
#define SortedAll  "SortedAll.txt"

int _data [MAX_ROWS][MAX_COLS];		// 2000 rows of 1000 numbers to sort!
int _overflow = MAX_COLS - 8;		// overflow value for 1000/16
int _allData [MAX_ROWS][MAX_COLS];	// flat version of _data for outputting

const int rseed = 123;				// arbitrary seed for random number generator - PLEASE DON'T ALTER
									// After sorting data generated with seed 123 these results should be true:
const int	checkBeg = 87,			// at [0][0]
			checkMid = 16440,		// at [MAX_ROWS/2][MAX_COLS/2]
			checkEnd = 32760;		// at [MAX_ROWS-1][MAX_COLS-1]

CStopWatch s1, s2, s3, s4;

void getData(void);
void sortEachRow(void);
void displayCheckData(void);
void sortAll(void);
void outputTimes(void);
void outputDataAsString(int outputData[][MAX_COLS], string file);
inline void mergeSort(int iArray[], int iBeginning, int iEnd);
inline void SIMDitoa16(int * iArray, char * oNumString);
inline void SIMDitoa8 (int * iArray, char * oNumString);
inline void radixSort(int arr[]);
inline void radixSortAll(int arr[][MAX_COLS]);


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
	outputDataAsString(_data, SortedRows);
	s2.stopTimer();

	s3.startTimer();
	sortAll();
	s3.stopTimer();

	s4.startTimer();
	outputDataAsString(_allData, SortedAll);
	s4.stopTimer();

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
	//basic parallel for to speed up individual rows
	#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; i++)
	{
		radixSort(_data[i]);
	}

}

void sortAll()
{
	//runs a full radix without omp as it applies to the entire array
	radixSortAll(_data);
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

void outputDataAsString(int outputData[][MAX_COLS], string file)
{
	char numString[MAX_CHARS * 16];
	string odata;

	for (int i = 0; i<MAX_ROWS; ++i){
		for (int j = 0; j<MAX_COLS; j+=16){
			//determines wether to do a 8 or 16 simd operation depending on how many left
			if (j != _overflow)
			{
				SIMDitoa16(&outputData[i][j], numString);
			}
			else
			{
				SIMDitoa8(&outputData[i][j], numString);
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
	fopen_s(&sodata, file.c_str(), "w");
	fputs(odata.c_str(), sodata);
	fclose(sodata);
}

#pragma region RemovedMergeSort
//Ref: https://www.comrevo.com/2016/02/openmp-program-for-merge-sort.html
//Removed as it seems slower than a radix sort

inline void merge(int iArray[], int iTemp[], int iBeginning, int iMiddle, int iEnd)
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

		int temp[MAX_COLS];
		merge(iArray, temp, iBeginning, middle, iEnd);
	}

}

#pragma endregion RemovedMergeSort

#pragma region SIMDitoa16&8

//*********************************************************************************
//Uses 256 bit SIMD registers to convert 16 numbers to ascii at a time
inline void SIMDitoa16(int *iArray, char * oNumString)
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
	__m256i magicNumber;
	__m256i inputData = startingArray;
	__m256i storedData;

	for (int i = 0; i < 5; i++)
	{
		//starts the 2 vars which will be iterated upon
		outputData = inputData;
		storedData = outputData;

		//setup and multiply by the first magic number
		magicNumber = _mm256_set1_epi16(ASCII_MAGIC_1);
		outputData = _mm256_mulhi_epu16(magicNumber, outputData);

		//shift the data right
		outputData = _mm256_srai_epi16(outputData, 2);

		//store these numbers for future iterations
		inputData = outputData;

		//setup and multiply by 10
		magicNumber = _mm256_set1_epi16(10);
		outputData = _mm256_mullo_epi16(magicNumber, outputData);

		//subtract the generated value from the original value
		outputData = _mm256_sub_epi16(storedData, outputData);

		//setup and add the second magic number to complete single
		//digit to char conversion
		magicNumber = _mm256_set1_epi16(ASCII_MAGIC_2);
		outputData = _mm256_add_epi16(magicNumber, outputData);

		//sort letters into correct 16 output arrays
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
inline void SIMDitoa8(int *iArray, char * oNumString)
{
	//setup a memory aligned array with the 8 input ints
	_declspec(align(16)) short intArray[8] = {
		iArray[0], iArray[1],  iArray[2],  iArray[3], iArray[4], iArray[5], iArray[6], iArray[7] 
	};
	//declare some simd variables to use in the loop
	__m128i startingArray = _mm_load_si128((__m128i*)intArray);
	__m128i outputData;
	__m128i magicNumber;
	__m128i inputData = startingArray;
	__m128i storedData;

	for (int i = 0; i < 5; i++)
	{
		//starts the 2 vars which will be iterated upon
		outputData = inputData;
		storedData = outputData;

		//setup and multiply by the first magic number
		magicNumber = _mm_set1_epi16(ASCII_MAGIC_1);
		outputData = _mm_mulhi_epi16(magicNumber, outputData);

		//shift the data right
		outputData = _mm_srai_epi16(outputData, 2);

		//store these numbers for future iterations
		inputData = outputData;

		//setup and multiply by 10
		magicNumber = _mm_set1_epi16(10);
		outputData = _mm_mullo_epi16(magicNumber, outputData);

		//subtract the generated value from the original value
		outputData = _mm_sub_epi16(storedData, outputData);

		//setup and add the second magic number to complete single
		//digit to char conversion
		magicNumber = _mm_set1_epi16(ASCII_MAGIC_2);
		outputData = _mm_add_epi16(magicNumber, outputData);

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
// determines the highest number in the array
int getMax(int iArray[])
{
	int max = iArray[0];
	#pragma omp parallel for
	for (int i = 1; i < MAX_ROWS; i++)
	{
		if (iArray[i] > max)
		{
			max = iArray[i];
		}
	}
	return max;
}

// Does the main sort of the radix sort
void countSort(int arr[], int exp)
{
	int output[MAX_COLS], // output array
		i, 
		count[10] = { 0 };

	//fills up the buckets
	for (i = 0; i < MAX_COLS; i++)
	{
		count[(arr[i] / exp) % 10]++;
	}

	// Change count[i] so that count[i] now contains actual
	//  position of this digit in output[]
	for (i = 1; i < 10; i++)
		count[i] += count[i - 1];

	// Build the output array
	for (i = MAX_COLS - 1; i >= 0; i--)
	{
		output[count[(arr[i] / exp) % 10] - 1] = arr[i];
		count[(arr[i] / exp) % 10]--;
	}

	//Copy the output array to the initial array, to have sorted numbers
	#pragma omp parallel for
	for (i = 0; i < MAX_COLS; i++)
	{
		arr[i] = output[i];
	}
}

// Radix Sort
void radixSort(int arr[])
{
	// Find the maximum number to know number of digits
	int m = getMax(arr);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
		countSort(arr, exp);
}

#pragma endregion RowSorting

#pragma region AllSorting
// determines the highest number in the 2D array
int getMaxAll(int iArray[][MAX_COLS])
{
	int max = iArray[0][0];
	#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; i++)
	{
		for (int j = 0; j < MAX_COLS; j++)
		{
			if (iArray[i][j] > max)
			{
				max = iArray[i][j];
			}
		}
	}
	return max;
}

// Does the main sort of the radix sort
void countSortAll(int iArray[][MAX_COLS], int exp)
{
	int count[10] = { 0 };

	//create 4 sets of buckets
	int bucket1[10] = { 0 },
		bucket2[10] = { 0 },
		bucket3[10] = { 0 },
		bucket4[10] = { 0 };

	int sectionSize = MAX_ROWS / 4;

	// Create 4 sections, and fill the 4 sets of buckets
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i = 0; i < sectionSize; i++)
			{
				for (int j = 0; j < MAX_COLS; j++)
				{
					bucket1[(iArray[i][j] / exp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize; i < sectionSize * 2; i++)
			{
				for (int j = 0; j < MAX_COLS; j++)
				{
					bucket2[(iArray[i][j] / exp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 2; i < sectionSize * 3; i++)
			{
				for (int j = 0; j < MAX_COLS; j++)
				{
					bucket3[(iArray[i][j] / exp) % 10]++;
				}
			}
		}

		#pragma omp section
		{
			for (int i = sectionSize * 3; i < MAX_ROWS; i++)
			{
				for (int j = 0; j < MAX_COLS; j++)
				{
					bucket4[(iArray[i][j] / exp) % 10]++;
				}
			}
		}
	}

	//add the buckets together
	for (int i = 0; i < 10; i++)
	{
		count[i] = bucket1[i] + bucket2[i] + bucket3[i] + bucket4[i];
	}

	// Change count[i] so that count[i] now contains actual
	// position of this digit in output[]
	for (int i = 1; i < 10; i++)
	{
		count[i] += count[i - 1];
	}

	//***********************************
	// Section 2 - No possible optimization due to array bring manpiulated
	//***********************************
	// Build the output array

	for (int i = MAX_ROWS - 1; i >= 0; i--)
	{
		for (int j = MAX_COLS - 1; j >= 0; j--)
		{
			int pos = count[(iArray[i][j] / exp) % 10] - 1;

			_allData[pos / MAX_COLS][pos % MAX_COLS] = iArray[i][j];
			count[(iArray[i][j] / exp) % 10]--;
		}
	}

	//***********************************
	// Section 3 - simple OMP FOR to speed up copying of arrays
	//***********************************
	// Copy output array to initial array
	#pragma omp parallel for
	for (int i = 0; i < MAX_ROWS; i++)
	{
		for (int j = 0; j < MAX_COLS; j++)
		{
			iArray[i][j] = _allData[i][j];
		}
	}
}

// Radix Sort Full
void radixSortAll(int arr[][MAX_COLS])
{
	// Find the maximum number to know number of digits
	int m = getMaxAll(arr);

	// Do counting sort for every digit. Note that instead
	// of passing digit number, exp is passed. exp is 10^i
	// where i is current digit number
	for (int exp = 1; m / exp > 0; exp *= 10)
	{
		countSortAll(arr, exp);
	}
}
#pragma endregion AllSorting