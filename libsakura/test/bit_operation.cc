#include <iostream>

#include <libsakura/sakura.h>
#include "gtest/gtest.h"

#define NUM_IN 8

using namespace std;

/*
 * Tests vaious bit operation of an uint8_t value and array
 * INPUTS:
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
class BitOperation8 : public ::testing::Test
{
protected:
	virtual void SetUp()
	{
		//cout << "SetUp method()" << endl;
		size_t ntype = 4;
		for (size_t i = 0; i < NUM_IN; i++){
			in_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
			bit_mask_ = 2; /* bit pattern of 0...010 */
			edit_mask_[i] = (i/ntype > 0); /*{F, F, F, F, T, T, T, T};*/
		}
	}
	virtual void TearDown()
	{
		//std::cout << "TearDown method()" << std::endl;
	}

	/* Converts an input uint8_t value to a bit pattern.*/
	char* B8ToS(uint8_t in_value) {
		int const bit_size = 8;
		static char buff[bit_size+1];
		buff[8] = '\0';
		for (size_t i = 0; i < bit_size ; i++){
			if((in_value>>i) % 2 == 0)
				buff[bit_size-1-i] = '0';
			else
				buff[bit_size-1-i] = '1';
		}
		//cout << "input = " << (int) in << " => " << buff << endl;
		return buff;
	}

	/* Converts an bit pattern (char) to a uint8_t value.*/
	uint8_t SToB8(char* in_bit) {
		uint8_t result(0);
		size_t i = 0;
		while (in_bit[i] != '\0') {
			//resutl = result * 2 + ((uint8_t) in_string[i]);
			result <<= 1;
			if (in_bit[i] == '1')
				result += 1;
			i++;
		}
		//cout << "INPUT (char) = " << in_string << "; OUTPUT (unsigned int) = " << (unsigned int) result << endl;
		return result;
	}

	void getResultUInt8(bool print, size_t num_in, char* const in_bit, uint8_t* result) {
		for (size_t i = 0 ; i < num_in ; i++){
			result[i] = SToB8(&in_bit[i]);
		}
		// Summary
		if (print){
			cout << "Conversion summary" << endl;
			cout << "Input (char) = [ ";
			for (size_t i = 0 ; i < num_in-1 ; i++)
				cout << in_bit[i] << ", ";
			cout << in_bit[num_in-1] << " ]" << endl;
			cout << "Result (unsigned int) = [ ";
			for (size_t i = 0 ; i < num_in-1 ; i++)
				cout << (unsigned int) result[i] << ", ";
			cout << (unsigned int) result[num_in-1] << " ]" << endl;
			cout << "Result (uint8_t) = [ '";
			for (size_t i = 0 ; i < num_in-1 ; i++)
				cout << (uint8_t) result[i] << "', '";
			cout << (uint8_t) result[num_in-1] << "' ]" << endl;
		}
	}

	void PrintInputs(){
		cout << "bit_mask = " << B8ToS(bit_mask_) << endl;
		PrintUInt8Array("in", NUM_IN, in_);
	}

	void PrintUInt8Array(char *name, size_t num_in, uint8_t *in){
		cout << name << " = [";
		for (size_t i = 0; i < num_in-1; i++)
			cout << B8ToS(in[i]) << ", " ;
		cout << B8ToS(in[num_in-1]) <<" ]" << endl;
	}

	uint8_t in_[NUM_IN];
	uint8_t bit_mask_;
	bool edit_mask_[NUM_IN];

};

/*
 * Test bit operation AND by sakura_OperateBitsUint8And
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation8, And) {
	//std::cout << "BitOperation8.And" << std::endl;
	uint8_t out[NUM_IN];
	uint8_t result[NUM_IN] = {0, 1, 2, 3, 0, 0, 2, 2};
	//uint8_t result[NUM_IN];
	size_t const num_in(NUM_IN);
	bool verbose(true);
	// conversion between bit pattern to uint8_t
	//char result_bit[NUM_IN][9] = {"00000000", "00000001", "00000010", "00000011", "00000000", "00000000", "00000010", "00000010"};
	//getResultUInt8(true, NUM_IN, result_bit, result);

	if (verbose) PrintInputs();

	sakura_OperateBitsUint8And(bit_mask_, num_in, in_, edit_mask_, out);

	if (verbose) PrintUInt8Array("out", num_in, out);

	// Verification
	for (size_t i = 0 ; i < num_in ; i++){
		ASSERT_EQ(out[i], result[i]);
	}
}
