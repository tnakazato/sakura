#ifndef _SAKURA_E2E_UTILS_H_
#define _SAKURA_E2E_UTILS_H_

#include <memory>

#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableParse.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/tables/Tables/TableCopy.h>
#include <casacore/tables/Tables/TableRow.h>

#include <libsakura/sakura.h>

#define ELEMENTSOF(x) (sizeof(x) / sizeof((x)[0]))
#define STATIC_ASSERT(x) static_assert((x), # x)

namespace {
#ifdef __GNUG__
/**
 * @~japanese
 * @brief @a ptr がsakuraのアライメント要件を満たしていると見なし、そのアドレスを返す。
 *
 * コンパイラがサポートしていれば、
 * コンパイラは、戻り値がsakuraのアライメント要件を満たしているものとして最適化を行う。
 *
 * @param ptr sakuraのアライメント要件を満たしているアドレス
 * @return @a ptr (sakuraのアライメント要件を満たしているというコンパイラ依存の属性付き)
 */
template<typename T>
inline T AssumeAligned(T ptr) {
	return reinterpret_cast<T>(__builtin_assume_aligned(ptr,
	32));
}
#else /* __GNUG__ */
/**
 * @~japanese
 * @brief @a ptr がsakuraのアライメント要件を満たしていると見なし、そのアドレスを返す。
 *
 * コンパイラがサポートしていれば、
 * コンパイラは、戻り値がsakuraのアライメント要件を満たしているものとして最適化を行う。
 *
 * @param ptr sakuraのアライメント要件を満たしているアドレス
 * @return @a ptr (sakuraのアライメント要件を満たしているというコンパイラ依存の属性付き)
 */
template<typename T>
inline /*alignas(LIBSAKURA_ALIGNMENT)*/T *AssumeAligned(T *ptr) {
	return ptr;
}
#endif /* __GNUG__ */

class ScopeGuard {
	typedef std::function<void(void) noexcept> Func;
public:
	ScopeGuard() = delete;
	explicit ScopeGuard(Func clean_up, bool enabled = true) :
	clean_up_(clean_up), engaged_(enabled), called_(false) {
	}
	ScopeGuard(ScopeGuard const &other) = delete;
	ScopeGuard &operator =(ScopeGuard const &other) = delete;
	ScopeGuard const &operator =(ScopeGuard const &other) const = delete;
	void *operator new(std::size_t) = delete;

	/**
	 * @~
	 * @brief Calls @ clean_up parameter provided to the constructor
	 * if the ScopeGuard is enabled.
	 */
	~ScopeGuard() {
		if (engaged_) {
			assert(! called_);
			clean_up_();
			// called_ = true;
		}
	}
	/**
	 * @~
	 * @brief Disables the ScopeGuard
	 *
	 * Destructor won't call @a clean_up parameter provided to the constructor
	 * if the ScopeGuard is disabled.
	 */
	void Disable() {
		engaged_ = false;
	}
	/**
	 * @~
	 * @brief Enables the ScopeGuard if it has not cleaned up yet
	 *
	 * Don't call @a Enable() after explicit call of @ref CleanUpNow().
	 */
	void Enable() {
		assert(! called_);
		engaged_ = true;
	}

	/**
	 * @~
	 * @brief Calls @a clean_up parameter provided to the constructor
	 * if the ScopeGuard is enabled
	 *
	 * Calling this method more than once is not allowed.
	 */
	void CleanUpNow() {
		if (engaged_) {
			assert(! called_);
			clean_up_();
			called_ = true;
			engaged_ = false;
		}
	}
private:
	Func clean_up_;
	bool engaged_;
	bool called_;
};

/**
 * Aligned array utility
 */
struct Deleter {
	inline void operator()(void *ptr) const {
		free(ptr);
	}
};

class AlignedArrayGenerator {
public:
	AlignedArrayGenerator() :
			alignment_(sakura_GetAlignment()), index_(0), pointer_holder_(128) {
	}
	/**
	 * MT-unsafe
	 */
	template<class T> inline T *GetAlignedArray(size_t num_elements)
			throw (std::bad_alloc) {
		size_t num_arena = num_elements * sizeof(T) + alignment_ - 1;
		if (index_ >= pointer_holder_.size()) {
			pointer_holder_.resize(pointer_holder_.size() + 128);
		}
		void *ptr = malloc(num_arena);
		if (ptr == nullptr) {
			throw std::bad_alloc();
		}
		pointer_holder_[index_++].reset(ptr);
		return AssumeAligned(
				reinterpret_cast<T *>(sakura_AlignAny(num_arena, ptr,
						num_elements * sizeof(T))));
	}
	inline size_t index() {
		return index_;
	}
	inline void *Release(size_t index) {
		return pointer_holder_[index].release();
	}
	inline void Transfer(size_t index, AlignedArrayGenerator *other) {
		if (index_ >= pointer_holder_.size()) {
			pointer_holder_.resize(pointer_holder_.size() + 128);
		}
		pointer_holder_[index_++].reset(other->Release(index));
	}
private:
	size_t alignment_;
	size_t index_;
	std::vector<std::unique_ptr<void, Deleter> > pointer_holder_;
};

/**
 * CASA related utility functions
 */
void CreateOutputTable(std::string const input_table,
		std::string const output_table) {
	auto func_start = sakura_GetCurrentTime();
	LOG4CXX_INFO(logger, "CreateOutputTable: Enter");
	casa::Table in(casa::String(input_table), casa::Table::Old);

	{
		auto start = sakura_GetCurrentTime();
		// create empty output table
		in.deepCopy(casa::String(output_table), casa::Table::New, casa::True,
				in.endianFormat(), casa::False);
		auto end = sakura_GetCurrentTime();
		LOG4CXX_INFO(logger, "CreateOutputTable: deep copied: " << end - start);
	}

	// copy subtables
	casa::Table out(casa::String(output_table), casa::Table::Update);
	{
		auto start = sakura_GetCurrentTime();
		casa::TableRecord const keyword_set = in.keywordSet();
		casa::TableRecord rw_keyword_set = out.rwKeywordSet();
		for (casa::uInt i = 0; i < keyword_set.nfields(); ++i) {
			if (keyword_set.type(i) == casa::DataType::TpTable) {
				casa::Table in_subtable = keyword_set.asTable(i);
				casa::Table out_subtable = rw_keyword_set.asTable(keyword_set.name(i));
				casa::TableCopy::copyRows(out_subtable, in_subtable);
			}
		}
		auto end = sakura_GetCurrentTime();
		LOG4CXX_INFO(logger, "CreateOutputTable: sub-tables copied: " << end - start);
	}

	// copy main table except SPECTRA, FLAGTRA, and TSYS
	{
		auto start = sakura_GetCurrentTime();
		out.addRow(in.nrow());
		casa::Vector<casa::String> excludes(3);
		excludes[0] = "SPECTRA";
		excludes[1] = "FLAGTRA";
		excludes[2] = "TSYS";
		casa::ROTableRow input_row(in, excludes, casa::True);
		casa::TableRow output_row(out, excludes, casa::True);
		for (casa::uInt i = 0; i < in.nrow(); ++i) {
			input_row.get(i);
			output_row.putMatchingFields(i, input_row.record());
		}
		auto end = sakura_GetCurrentTime();
		LOG4CXX_INFO(logger, "CreateOutputTable: main table copied: " << end - start);
	}
	auto func_end = sakura_GetCurrentTime();
	LOG4CXX_INFO(logger, "CreateOutputTable: Leave: " << func_end - func_start);
}

std::string GetTaqlString(std::string table_name, unsigned int ifno,
		unsigned int polno, bool on_source_only) {
	std::ostringstream oss;
	oss << "SELECT FROM \"" << table_name << "\" WHERE IFNO == " << ifno
			<< " && POLNO == " << polno;
	if (on_source_only) {
		oss << " && SRCTYPE == 0";
	}
	oss << " ORDER BY TIME";
	return oss.str();
}

casa::Table GetSelectedTable(std::string table_name, unsigned int ifno,
		unsigned int polno, bool on_source_only) {
	casa::String taql(GetTaqlString(table_name, ifno, polno, on_source_only));
	return tableCommand(taql);
}

void GetFrequencyLabelFromScantable(casa::Table const table, unsigned int ifno,
		size_t num_channel, double *frequency_label) {
	// Get FREQ_ID
	casa::Table const selected_table = table(table.col("IFNO") == ifno);
	if (selected_table.nrow() == 0) {
		return;
	}
	casa::TableColumn column(selected_table, "FREQ_ID");
	unsigned int freq_id = column.asuInt(0);

	// The table must have FREQUENCIES subtable
	casa::String const subtable_name = "FREQUENCIES";
	casa::TableRecord const header = selected_table.keywordSet();
	if (!header.isDefined(subtable_name)) {
		return;
	}
	casa::Table const frequencies_table = header.asTable(subtable_name);
	casa::Table const target_row = frequencies_table(
			frequencies_table.col("ID") == freq_id);
	if (target_row.nrow() == 0) {
		return;
	}
	auto GetValue = [&](casa::String const column_name) {
		column.attach(target_row, column_name);
		return column.asdouble(0);
	};
	double reference_pixel = GetValue("REFPIX");
	double reference_value = GetValue("REFVAL");
	double increment = GetValue("INCREMENT");
	for (size_t i = 0; i < num_channel; ++i) {
		frequency_label[i] = (static_cast<double>(i) - reference_pixel)
				* increment + reference_value;
	}
}

void GetFrequencyLabelFromScantable(std::string const table_name,
		unsigned int ifno, size_t num_channel, double *frequency_label) {
	GetFrequencyLabelFromScantable(casa::Table(table_name, casa::Table::Old),
			ifno, num_channel, frequency_label);
}

template<class T>
void GetScalarCells(T *array, size_t start_row, size_t nrow,
		casa::ROScalarColumn<T> const column) {
	casa::Vector < T > casa_array(casa::IPosition(1, nrow), array, casa::SHARE);
	casa::Slicer slice(casa::IPosition(1, start_row), casa::IPosition(1, nrow),
			casa::Slicer::endIsLength);
	column.getColumnRange(slice, casa_array);
}

template<class T>
void GetArrayCell(T *array, size_t row_index,
		casa::ROArrayColumn<T> const column, casa::IPosition const cell_shape) {
	casa::Array < T > casa_array(cell_shape, array, casa::SHARE);
	column.get(row_index, casa_array);
}

template<class T>
void GetArrayColumn(T *array, casa::ROArrayColumn<T> const column,
		casa::IPosition const column_shape) {
	casa::Array < T > casa_array(column_shape, array, casa::SHARE);
	column.getColumn(casa_array);
}

template<class T>
void GetScalarColumn(T *array, casa::ROScalarColumn<T> const column,
		size_t num_rows) {
	casa::Vector < T
			> casa_array(casa::IPosition(1, num_rows), array, casa::SHARE);
	column.getColumn(casa_array);
}

void GetFromCalTable(std::string const table_name, unsigned int ifno,
		unsigned int polno, std::string const data_name,
		AlignedArrayGenerator *array_generator, float **data,
		double **timestamp, size_t *num_data, size_t *num_row) {
	casa::Table table = GetSelectedTable(table_name, ifno, polno, false);
	*num_row = table.nrow();
	if (*num_row == 0) {
		throw casa::AipsError("Selected sky table is empty.");
	}
	casa::ROArrayColumn<float> caldata_column(table, data_name);
	casa::ROScalarColumn<double> time_column(table, "TIME");
	*num_data = caldata_column(0).nelements();
	*data = array_generator->GetAlignedArray<float>((*num_data) * (*num_row));
	*timestamp = array_generator->GetAlignedArray<double>(*num_row);
	GetArrayColumn(*data, caldata_column,
			casa::IPosition(2, *num_data, *num_row));
	GetScalarColumn(*timestamp, time_column, *num_row);
}
}

#endif /* _SAKURA_E2E_UTILS_H_ */
