#ifndef _SAKURA_E2E_UTILS_H_
#define _SAKURA_E2E_UTILS_H_

#include <memory>

#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/Vector.h>

#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/TableParse.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <libsakura/sakura.h>

namespace {

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
	template<class T> inline T *GetAlignedArray(size_t num_elements) throw (std::bad_alloc) {
		size_t num_arena = num_elements * sizeof(T) + alignment_ - 1;
		if (index_ >= pointer_holder_.size()) {
			pointer_holder_.resize(pointer_holder_.size() + 128);
		}
		void *ptr = malloc(num_arena);
		if (ptr == nullptr) {
			throw std::bad_alloc();
		}
		pointer_holder_[index_++].reset(ptr);
		return reinterpret_cast<T *>(sakura_AlignAny(num_arena, ptr,
				num_elements * sizeof(T)));
	}
	inline size_t index() {
		return index_;
	}
	inline void *release(size_t index) {
		return pointer_holder_[index].release();
	}
private:
	size_t alignment_;
	size_t index_;
	std::vector<std::unique_ptr<void, Deleter> > pointer_holder_;
};

/**
 * CASA related utility functions
 */
std::string GetTaqlString(std::string table_name, unsigned int ifno) {
	std::ostringstream oss;
	oss << "SELECT FROM \"" << table_name << "\" WHERE IFNO == " << ifno
			<< " ORDER BY TIME";
	return oss.str();
}

casa::Table GetSelectedTable(std::string table_name, unsigned int ifno) {
	casa::String taql(GetTaqlString(table_name, ifno));
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
		std::string const data_name, AlignedArrayGenerator *array_generator,
		float **data, double **timestamp, size_t *num_data, size_t *num_row) {
	casa::Table table = GetSelectedTable(table_name, ifno);
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
