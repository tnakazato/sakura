#include <iostream>
#include <string>

#include <libsakura/sakura.h>
#include "aligned_memory.h"
#include "gtest/gtest.h"

/* the number of elements in input/output array to test */
#define NUM_IN 8

using namespace std;

/*
 * A super class to test various bit operation of an value and array
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
template<typename DataType>
class BitOperation : public ::testing::Test
{
protected:

	BitOperation()
	: verbose(false)
	{
		bit_size = sizeof(DataType)*8;
	}

	virtual void SetUp()
	{
		size_t const ntype(4);
		for (size_t i = 0; i < NUM_IN; i++){
			in_[i] = i % ntype; /* repeat bit pattern of *00, *01, *10, *11,... */
			bit_mask_ = 2; /* bit pattern of 0...010 */
			edit_mask_[i] = (i/ntype > 0); /*{F, F, F, F, T, T, T, T};*/
		}

		// Initialize sakura
		LIBSAKURA_SYMBOL(Status) status = LIBSAKURA_SYMBOL(Initialize)();
		EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	}

	virtual void TearDown()
	{
		// Clean-up sakura
		LIBSAKURA_SYMBOL(CleanUp)();
	}

	/* Converts an input value to a bit pattern.*/
//	char* BToS(DataType in_value) {
	string BToS(DataType in_value) {
		char buff[bit_size+1];
		buff[bit_size] = '\0';
		for (size_t i = 0; i < bit_size ; i++){
			if((in_value>>i) % 2 == 0)
				buff[bit_size-1-i] = '0';
			else
				buff[bit_size-1-i] = '1';
		}
		return string(buff);
	}

	/* Converts an bit pattern (char) to a value of DataType.*/
	DataType SToB(char* in_bit) {
		DataType result(0);
		size_t i = 0;
		while (in_bit[i] != '\0') {
			//result = result * 2 + ((uint8_t) in_string[i]);
			result <<= 1;
			if (in_bit[i] == '1')
				result += 1;
			i++;
		}
		return result;
	}

	void PrintInputs(){
		cout << "bit_mask = " << BToS(bit_mask_) ;
		cout << endl;
		PrintArray("in", NUM_IN, in_);
	}

	void PrintArray(char const *name, size_t num_in, DataType *in){
		cout << name << " = [";
		for (size_t i = 0; i < num_in-1; i++)
			cout << BToS(in[i]) << ", " ;
		cout << BToS(in[num_in-1]) ;
		cout << " ]" << endl;
	}

	SIMD_ALIGN DataType in_[NUM_IN];
	DataType bit_mask_;
	SIMD_ALIGN bool edit_mask_[NUM_IN];
	bool verbose;
	size_t bit_size;
	//size_t const bit_size = sizeof(DataType)*8;

};

/*
 * Tests various bit operation of an uint32_t value and array
 * INPUTS:
 * - bit_mask = 0...010
 * - in = [ 0...000, 0...001, 0...010, 0...011, 0...000, 0...001, 0...010, 0...011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
class BitOperation32 : public BitOperation<uint32_t>
{};

/*
 * Tests various bit operation of an uint8_t value and array
 * INPUTS:
 * - bit_mask = 00000010
 * - in = [ 00000000, 00000001, 00000010, 00000011, 00000000, 00000001, 00000010, 00000011 ]
 * - edit_mask = [F, F, F, F, T, T, T, T]
 */
class BitOperation8 : public BitOperation<uint8_t>
{};


/*
 * Test bit operation AND by sakura_OperateBitsUint8And
 * RESULT:
 * out = [00000000, 00000001, 00000010, 00000011, 00000000, 00000000, 00000010, 00000010 ]
 */
TEST_F(BitOperation8, And) {
	SIMD_ALIGN uint8_t out[NUM_IN];
	uint8_t result[NUM_IN] = {0, 1, 2, 3, 0, 0, 2, 2};
	//uint8_t result[NUM_IN];
	size_t const num_in(NUM_IN);
	// conversion between bit pattern to uint8_t
	//char result_bit[NUM_IN][9] = {"00000000", "00000001", "00000010", "00000011", "00000000", "00000000", "00000010", "00000010"};
	//getResultUInt8(true, NUM_IN, result_bit, result);

	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint8And(bit_mask_, num_in, in_, edit_mask_, out);

	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_in ; i++){
		ASSERT_EQ(out[i], result[i]);
	}
}

/*
 * Test bit operation AND by sakura_OperateBitsUint32And
 * RESULT:
 * out = [0...000, 0...001, 0...010, 0...011, 0...000, 0...000, 0...010, 0...010 ]
 */
TEST_F(BitOperation32, And) {
	SIMD_ALIGN uint32_t out[NUM_IN];
	uint32_t result[NUM_IN] = {0, 1, 2, 3, 0, 0, 2, 2};
	size_t const num_in(NUM_IN);
	//uint32_t result[NUM_IN];
	// conversion between bit pattern to uint8_t
	//char result_bit[NUM_IN][9] = {"00000000", "00000001", "00000010", "00000011", "00000000", "00000000", "00000010", "00000010"};
	//getResultUInt(true, NUM_IN, result_bit, result);

	if (verbose) PrintInputs();

	LIBSAKURA_SYMBOL(Status) status = sakura_OperateBitsUint32And(bit_mask_, num_in, in_, edit_mask_, out);

	if (verbose) PrintArray("out", num_in, out);

	// Verification
	EXPECT_EQ(LIBSAKURA_SYMBOL(Status_kOK), status);
	for (size_t i = 0 ; i < num_in ; i++){
		ASSERT_EQ(out[i], result[i]);
	}
}
