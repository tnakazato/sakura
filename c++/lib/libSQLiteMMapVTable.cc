#include <new>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <cassert>
#include <cstdarg>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <vector>
#include <fcntl.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <pthread.h>
#include <sqlite3ext.h>

SQLITE_EXTENSION_INIT1

using namespace std;

namespace a {

#define unique_ptr auto_ptr
#define elementsof(x) (sizeof(x) / sizeof(*(x)))
#define sizeofMember(t,m) (sizeof(reinterpret_cast<t*>(0)->m))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define PERROR(x) perror(__FILE__  ":" #x)
#define enter() do { \
    if (false) {\
      cout << "Enter: " << __PRETTY_FUNCTION__ << endl;\
    }\
}while(0)
#define leave() do{ \
    if (false) {\
      cout << "Leave: " << __PRETTY_FUNCTION__ << endl;\
    }\
}while(0)

#define LOG if (false)

  template<typename ArrayElementType>
  class ArrayReleaser {
    ArrayElementType *str;
  public:
    ArrayReleaser(ArrayElementType *array)
      : str(array) {}
    virtual ~ArrayReleaser() {
      if (str) {
	delete[] str;
      }
    }
    ArrayElementType *get() const {
      return str;
    }
    ArrayElementType &operator [] (size_t idx) const {
      return str[idx];
    }
    ArrayElementType *release() {
      ArrayElementType *tmp = str;
      str = NULL;
      return str;
    }
  };

  char const *StrDup(char const *str) throw(bad_alloc) {
    assert(str != NULL);
    size_t len = strlen(str);
    char *p = new char[len + 1];
    strcpy(p, str);
    return p;
  }

  char const *StrNDup(char const *str, size_t len_) throw(bad_alloc) {
    assert(str != NULL);
    size_t len = min(strlen(str), len_);
    char *p = new char[len + 1];
    strncpy(p, str, len);
    p[len] = '\0';
    return p;
  }

  char *vmprintf(char const format[], va_list ap) {
    char dummy[1];
    va_list ap2;
    va_copy(ap2, ap);
    int len = vsnprintf(dummy, 0, format, ap);
    // va_end(ap); caller is responsible.
    char *result = NULL;
    if (len >= 0) {
      result = new char[len + 1];
      int len2 = vsnprintf(result, len + 1, format, ap2);
      assert(len == len2);
    }
    va_end(ap2);
    return result;
  }

  char *mprintf(char const format[], ...) {
    va_list ap;
    va_start(ap, format);
    char *result = vmprintf(format, ap);
    va_end(ap);
    return result;
  }

  string strprintf(char const format[], ...) {
    va_list ap;
    va_start(ap, format);
    ArrayReleaser<char const> tmp(vmprintf(format, ap));
    va_end(ap);
    string result(tmp.get());
    return result;
  }

  template <typename T>
  void fillZero(T data[], size_t n) {
    for (size_t i = 0; i < n; i++) {
      data[i] = 0;
    }
  }

  template <typename T>
  void deleteArray(void *data) {
    delete[] (T*)data;
  }

  class RTException {
    char const *msg;
  public:
    RTException(char const *staticString = "") : msg(staticString) {}
    RTException(RTException const &other) : msg(other.msg) {}
    char const *getMessage() { return msg; }
    virtual ~RTException() {}
  };

  RTException ASSERTION_ERROR("Assertion error");

  class Openable {
  public:
    virtual void open() throw(RTException) = 0;
    virtual ~Openable() {}
  };

  class Closable {
  public:
    virtual void close() throw(RTException) = 0;
    virtual ~Closable() {}
  };

  class Closer {
    Closable *closable_;
    Closer(Closer const &other) : closable_(NULL) {
      assert(false);
      throw ASSERTION_ERROR;
    }
    void close() throw(RTException) {
      if (closable_) {
	closable_->close();
	closable_ = NULL;
      }
    }
  public:
    Closer(Closable *closable) : closable_(closable) {
      assert(closable_ != NULL);
    }
    Closer() : closable_(NULL) {
    }
    virtual ~Closer() throw(RTException) {
      close();
    }

    void reset(Closable *closable) throw(RTException) {
      try {
	close();
      } catch (...) {
	if (closable) {
	  closable->close();
	}
	throw;
      }
      closable_ = closable;
    }

    void release() {
      closable_ = NULL;
    }
  };

  class Mutex;
  class LockHolder {
    Mutex *mutex;
    friend class Mutex;
    void lock(Mutex *m) throw(RTException);
  public:
    LockHolder(): mutex(NULL) {}
    virtual ~LockHolder() throw(RTException);
  };

  class Mutex: public virtual Closable {
    pthread_mutex_t mutex;
    Mutex(Mutex const &other) {
      assert(false);
      throw ASSERTION_ERROR;
    }
    Mutex &operator =(Mutex const &other) {
      assert(false);
      throw ASSERTION_ERROR;
    }
  public:
    Mutex() throw(RTException) {
      int result = pthread_mutex_init(&mutex, NULL);
      if (result != 0) {
	static RTException ex("pthread_mutex_init error");
	throw ex;
      }
    }

    virtual ~Mutex() throw(RTException) {
      int result = pthread_mutex_destroy(&mutex);
      if (result != 0) {
	static RTException ex("pthread_mutex_destroy error");
	throw ex;
      }
    }

    void lock() throw(RTException) {
      enter();
      int result = pthread_mutex_lock(&mutex);
      if (result != 0) {
	static RTException ex("pthread_mutex_lock error");
	throw ex;
      }
    }
    /**
     * Returns true if this thread could lock.
     * Returns false if already locked.
     */
    bool try_lock() throw(RTException) {
      enter();
      int result = pthread_mutex_trylock(&mutex);
      if (result == 0) {
	return true;
      }
      if (result == EBUSY) {
	return false;
      }
      static RTException ex("pthread_mutex_trylock error");
      throw ex;
    }

    void unlock() throw(RTException) {
      enter();
      int result = pthread_mutex_unlock(&mutex);
      if (result != 0) {
	static RTException ex("pthread_mutex_unlock error");
	throw ex;
      }
    }

    void takeLock(LockHolder &holder) throw(RTException) {
      holder.lock(this);
    }

    virtual void close() throw(RTException) {
      unlock();
    }
  };

  void LockHolder::lock(Mutex *m) throw(RTException) {
    mutex = m;
    if (mutex) {
      mutex->lock();
    }
  }

  LockHolder::~LockHolder() throw(RTException) {
    if (mutex) {
      mutex->unlock();
    }
  }

  size_t PAGE_SIZE;

  void init_module() throw(RTException) {
    long pageSize = sysconf(_SC_PAGE_SIZE);
    if (pageSize <= 0) {
      static RTException ex("Failed to get page size.");
      throw ex;
    }
    PAGE_SIZE = pageSize;
  }

  class File: public virtual Openable, public virtual Closable {
    friend class MMap;
    Mutex lock;
    string filename_;
    int mode_;
    int omode;
    int fd;
    size_t refCount;

    void ref() { 
      lock.lock();
      Closer autoUnlock(&lock);
      refCount++;
    }
    void unref() {
      lock.lock();
      Closer autoUnlock(&lock);
      assert(refCount > 0);
      refCount--;
    }
  public:
    enum Mode {
      Mode_Read = PROT_READ,
      Mode_Write = PROT_WRITE,
      Mode_None = PROT_NONE
    };
    File(char const filename[], int mode) throw(RTException)
      : filename_(filename), mode_(mode), omode(0), fd(-1), refCount(0) {
      enter();
      assert(filename != NULL);
      assert(mode_ == Mode_None || mode_ == Mode_Read || mode_ == Mode_Write
	     || mode_ == (Mode_Read|Mode_Write));
      switch (mode_) {
      case Mode_None:
	omode = O_RDONLY;
	break;
      case Mode_Read:
	omode = O_RDONLY;
	break;
      case Mode_Write:
	omode = O_WRONLY|O_CREAT;
	break;
      case Mode_Read|Mode_Write:
	omode = O_RDWR|O_CREAT;
	break;
      }
    }
    virtual ~File() throw(RTException) {
      enter();
      // without lock there should be no referer.
      if (fd >= 0) {
	close();
      }
      assert(refCount == 0);
    }
    int getMode() const {
      return mode_;
    }
    virtual void open() throw(RTException) {
      enter();
      lock.lock();
      Closer autoUnlock(&lock);
      if (fd == -1) {
	LOG cout << "omode: " << omode << endl;
	fd = ::open(filename_.c_str(), omode, 0666);
	if (fd < 0) {
	  PERROR(open);
	  static RTException ex("open error");
	  throw ex;
	}
      }
    }
    virtual void close() throw(RTException) {
      enter();
      lock.lock();
      Closer autoUnlock(&lock);
      if (fd >= 0) {
	int result = ::close(fd);
	if (result < 0) {
	  PERROR(close);
	  static RTException ex("close error");
	  throw ex;
	}
	fd = -1;
      }
    }
  };

  class MMap: public virtual Openable, public virtual Closable {
    Mutex lock;
    File *file_;
    void *mappedAddr;
    size_t size_;
    off_t fileOffset_;

    void extend() throw(RTException) {
      enter();
      static char const zero[4096] = {0};
      off_t fend = lseek(file_->fd, 0, SEEK_END);
      if (fend == static_cast<off_t>(-1)) {
	static RTException ex("failed to lseek");
	throw ex;
      }
      off_t tail = fileOffset_ + size_;
      if (fend < tail && tail > 1) {
	off_t newFend = lseek(file_->fd, tail - 1, SEEK_SET);
	if (newFend == static_cast<off_t>(-1)) {
	  static RTException ex("failed to lseek");
	  throw ex;
	}
	int wrote = write(file_->fd, zero, 1);
	if (wrote != 1) {
	  static RTException ex("failed to write");
	  throw ex;
	}
      }
    }

  public:
    MMap(File *file) throw(RTException)
      : file_(file), mappedAddr(MAP_FAILED), size_(0), fileOffset_(0) {
      enter();
      assert(file_ != NULL);
      file_->ref();
    }
    virtual ~MMap() throw(RTException) {
      enter();
      file_->unref();
      if (mappedAddr != MAP_FAILED) { // without lock because there should be no referer.
	unmap();
      }
    }

    static size_t getPageSize() {
      assert(PAGE_SIZE > 0);
      return PAGE_SIZE;
    }

    void setMapRegion(off_t offset, size_t size) throw(RTException) {
      enter();
      assert(offset >= 0 && offset % getPageSize() == 0);
      if (mappedAddr != MAP_FAILED) {
	static RTException ex("unmap before changing map region.");
	throw ex;
      }
      fileOffset_ = offset;
      size_ = size;
    }

    bool inRange(off_t offset, size_t size) const {
      assert(offset >= 0);
      if (fileOffset_ <= offset && offset + size <= fileOffset_ + size_) {
	return true;
      }
      return false;
    }

    template <typename T>
    T *remapForRegion(off_t offset, size_t size) throw(RTException) {
      if (mappedAddr == MAP_FAILED || inRange(offset, size)) {
	unmap();
	off_t newOffset = (offset / getPageSize()) * getPageSize();
	size_t newSize = offset + size - newOffset;
	newSize = (newSize + getPageSize() - 1)
	  / getPageSize() * getPageSize();
	setMapRegion(newOffset, newSize);
	map<void>();
      }
      char *p = reinterpret_cast<char *>(mappedAddr);
      p += offset - fileOffset_;
      return reinterpret_cast<T *>(p);
    }

    template <typename T>
    T *getMappedAddr() throw(RTException) {
      if (mappedAddr == MAP_FAILED) {
	static RTException ex("trial of getting mapped address w/o mapping.");
	throw ex;
      }
      return reinterpret_cast<T *>(mappedAddr);
    }

    template <typename T>
    T *map() throw(RTException) {
      enter();
      if (mappedAddr == MAP_FAILED) {
	extend();
	mappedAddr = mmap(NULL, size_, file_->getMode(),
			  MAP_SHARED, file_->fd, fileOffset_);
	if (mappedAddr == MAP_FAILED) {
	  PERROR(mmap);
	  static RTException ex("mmap error");
	  throw ex;
	}
      }
      return reinterpret_cast<T *>(mappedAddr);
    }

    void unmap() throw(RTException) {
      enter();
      if (mappedAddr != MAP_FAILED) {
	if (munmap(mappedAddr, size_) != 0) {
	  PERROR(munmap);
	  static RTException ex("munmap error");
	  throw ex;
	}
	mappedAddr = MAP_FAILED;
      }
    }

    virtual void open() throw(RTException) {
      map<void>();
    }
    virtual void close() throw(RTException) {
      unmap();
    }

    void takeLock(LockHolder &holder) throw(RTException) {
      lock.takeLock(holder);
    }
  };

  enum SQLType {
    SQLTYPE_INTEGER = SQLITE_INTEGER,
    SQLTYPE_FLOAT = SQLITE_FLOAT,
    SQLTYPE_BLOB = SQLITE_BLOB,
    SQLTYPE_NULL = SQLITE_NULL,
    SQLTYPE_TEXT = SQLITE_TEXT
  };

  char const *toSQLTypeName(SQLType type) {
    switch (type) {
    case SQLTYPE_INTEGER:
      return "INTEGER";
    case SQLTYPE_FLOAT:
      return "REAL";
    case SQLTYPE_BLOB:
      return "BLOB";
    case SQLTYPE_NULL:
      return "NULL";
    case SQLTYPE_TEXT:
      return "TEXT";
    }
    assert(false);
    return "";
  }

  class Table;
  class Column: public virtual Openable, public virtual Closable {
  public:
    struct ColumnDesc {
      SQLType type;
      bool isPK;
      bool isNotNull;
      char name[32];
    };

  private:
    Table *table;
    string name_;
  protected:
    SQLType type_;
  private:
    size_t size;
    bool isPK_;
    bool notNull_;
  protected:
    unique_ptr<File> file;
    unique_ptr<MMap> map;

    File *openFile(char const suffix[]) throw(RTException);

  private:
    friend class Table;
    void init(char const name[]) {
      assert(strlen(name) < sizeofMember(ColumnDesc, name));
      size = 0;
      switch (type_) {
      case SQLTYPE_INTEGER:
	size = sizeof(sqlite_int64);
	break;
      case SQLTYPE_FLOAT:
	size = sizeof(double);
	break;
      case SQLTYPE_BLOB:
      case SQLTYPE_TEXT:
	size = 0;
	break;
      default:
	assert(false);
	throw ASSERTION_ERROR;
      }
    }

    static void setInteger(sqlite3_context *pCtx,
			   void const*ptr, size_t size) {
      sqlite3_result_int64(pCtx, *reinterpret_cast<sqlite_int64 const*>(ptr));
    }

    static void setReal(sqlite3_context *pCtx,
			void const*ptr, size_t size) {
      sqlite3_result_double(pCtx, *reinterpret_cast<double const*>(ptr));
    }

    static void setBlob(sqlite3_context *pCtx,
			void const*ptr, size_t size) {
      sqlite3_result_blob(pCtx, ptr, size, SQLITE_TRANSIENT);
    }

    static void setText(sqlite3_context *pCtx,
			void const*ptr, size_t size) {
      sqlite3_result_text(pCtx, reinterpret_cast<char const*>(ptr), size,
			  SQLITE_TRANSIENT);
    }

  public:
    typedef void (*ResultSetter)(sqlite3_context *pCtx,
				 void const*ptr, size_t size);
    Column(ColumnDesc const &colDesc)
      : table(NULL), name_(colDesc.name), type_(colDesc.type),
	isPK_(colDesc.isPK), notNull_(colDesc.isNotNull),
	file(), map() {
      init(colDesc.name);
    }
    Column(char const name[], SQLType type, bool isPk, bool notNull)
      : table(NULL), name_(name), type_(type),
	isPK_(isPk), notNull_(notNull),
	file(), map() {
      init(name);
    }
    virtual ~Column() {
    }

    void save(ColumnDesc *colDesc) const {
      assert(colDesc != NULL);
      colDesc->type = type_;
      colDesc->isPK = isPK_;
      colDesc->isNotNull = notNull_;
      memset(colDesc->name, '\0', sizeof(colDesc->name));
      strncpy(colDesc->name, name_.c_str(), sizeof(colDesc->name) - 1);
    }

    virtual void open() throw(RTException) {
      enter();
      if (isPK()) {
	return;
      }
      file.reset(openFile(".0"));
      try {
	file->open();
	try {
	  map.reset(new MMap(file.get()));
	} catch (...) {
	  file->close();
	  throw;
	}
      } catch (...) {
	file.reset(NULL);
	throw;
      }
    }
    virtual void close() throw(RTException) {
      enter();
      if (isPK()) {
	return;
      }
      try {
	if (map.get() != NULL) {
	  map->close();
	}
      } catch (...) {
	if (file.get() != NULL) {
	  file->close();
	}
	throw;
      }
      if (file.get() != NULL) {
	file->close();
      }
    }

    virtual void insert(sqlite_int64 rowId, sqlite3_value *value) throw(RTException) {
      enter();
      assert(notNull_);
      assert(rowId > 0);
      assert(size > 0);
      assert(sqlite3_value_type(value) != SQLITE_NULL);
      {
	LockHolder protect;
	map->takeLock(protect);
	off_t offset = (rowId - 1) * size;
	switch (type_) {
	case SQLTYPE_INTEGER: {
	  assert(sqlite3_value_type(value) == SQLITE_INTEGER);
	  sqlite_int64 *p = map->remapForRegion<sqlite_int64>(offset, size);
	  *p = sqlite3_value_int64(value);
	  break;
	}
	case SQLTYPE_FLOAT: {
	  assert(sqlite3_value_type(value) == SQLITE_FLOAT);
	  double *p = map->remapForRegion<double>(offset, size);
	  *p = sqlite3_value_double(value);
	  break;
	}
	default:
	  assert(false);
	  throw ASSERTION_ERROR;
	}
      }
    }

    void fetch(sqlite3_context *pCtx, sqlite_int64 rowId) throw(RTException) {
      switch (type_) {
      case SQLTYPE_INTEGER:
	fetch(pCtx, rowId, setInteger);
	break;
      case SQLTYPE_FLOAT:
	fetch(pCtx, rowId, setReal);
	break;
      case SQLTYPE_BLOB:
	fetch(pCtx, rowId, setBlob);
	break;
      case SQLTYPE_TEXT:
	fetch(pCtx, rowId, setText);
	break;
      default:
	assert(false);
	throw ASSERTION_ERROR;
      }
    }
    virtual void fetch(sqlite3_context *pCtx, sqlite_int64 rowId,
		       ResultSetter setter)
      throw(RTException) {
      enter();
      assert(rowId > 0);
      assert(setter > 0);
      {
	LockHolder protect;
	map->takeLock(protect);
	off_t offset = (rowId - 1) * size;
	switch (type_) {
	case SQLTYPE_INTEGER:
	case SQLTYPE_FLOAT:
	  {
	    void const *p = map->remapForRegion<void const>(offset, size);
	    setter(pCtx, p, size);
	  }
	  break;
	default:
	  assert(false);
	  throw ASSERTION_ERROR;
	}
      }
    }

    string const &getName() const {
      return name_;
    }
    SQLType getType() const {
      return type_;
    }
    size_t getSize() const {
      return size;
    }
    bool isPK() const {
      return isPK_;
    }
    bool isNotNull() const {
      return notNull_;
    }
  };

  class VariableSizeColumn: public Column {
    unique_ptr<File> dataFile;
    unique_ptr<MMap> headerMap;
    unique_ptr<MMap> dataMap;
    struct Entry {
      off_t offset;
      size_t size;
    };
    struct DataHeader {
      off_t tail;
      char padding[256 - sizeof(off_t)];
      char data[1];
    };
    off_t alignUp(off_t offset) {
      off_t addr = offset;
      off_t mask = 32 - 1;
      addr += mask;
      addr &= ~mask;
      return addr;
    }
  public:
    VariableSizeColumn(ColumnDesc const &colDesc)
      :Column(colDesc), dataFile(), headerMap(), dataMap() {
    }
    VariableSizeColumn(char const name[], SQLType type, bool isPk, bool notNull)
      : Column(name, type, isPk, notNull), dataFile(), headerMap(), dataMap() {
    }
    virtual ~VariableSizeColumn() {
    }

    virtual void open() throw(RTException) {
      enter();
      Column::open();
      dataFile.reset(openFile(".data"));
      try {
	dataFile->open();
	Closer dataFileCloser(dataFile.get());
	dataMap.reset(new MMap(dataFile.get()));
	try {
	  headerMap.reset(new MMap(dataFile.get()));
	  try {
	    headerMap->setMapRegion(0, MMap::getPageSize());
	    DataHeader *header = headerMap->map<DataHeader>();
	    if (header->tail < static_cast<off_t>(offsetof(DataHeader, data))) {
	      header->tail = offsetof(DataHeader, data);
	      assert(alignUp(offsetof(DataHeader, data))
		     == offsetof(DataHeader, data));
	    }
	    dataFileCloser.release();
	  } catch (...) {
	    headerMap.reset(NULL);
	  }
	} catch (...) {
	  dataMap.reset(NULL);
	}
      } catch (...) {
	dataFile.reset(NULL);
	throw;
      }
    }
    virtual void close() throw(RTException) {
      enter();
      Closer dataFileCloser(dataFile.get());
      try {
	Closer headerMapCloser;
	Closer dataMapCloser;
	if (headerMap.get() != NULL) {
	  headerMapCloser.reset(headerMap.get());
	}
	if (dataMap.get() != NULL) {
	  dataMapCloser.reset(dataMap.get());
	}
	if (dataFile.get() != NULL) {
	  dataFile->close();
	}
      } catch (...) {
	Column::close();
	throw;
      }
      Column::close();
    }

    virtual void insert(sqlite_int64 rowId, sqlite3_value *value) throw(RTException) {
      enter();
      {
	LockHolder protect;
	LockHolder protectHeader;
	LockHolder protectData;
	map->takeLock(protect);
	headerMap->takeLock(protectHeader);
	dataMap->takeLock(protectData);

	assert(rowId > 0);
	off_t offset = (rowId - 1) * sizeof(Entry);
	Entry *entry = map->remapForRegion<Entry>(offset, sizeof(Entry));
	int type = sqlite3_value_type(value);
	if (type == SQLITE_NULL) {
	  entry->offset = 0; // means NULL SQL value
	  entry->size = 0;
	} else {
	  LOG cout << "Type: " << type << ", " << type_ << endl;
	  //assert(type == type_);
	  DataHeader *header = headerMap->getMappedAddr<DataHeader>();
	  off_t dataLocation = alignUp(header->tail);
	  int size = sqlite3_value_bytes(value);
	  void *data = dataMap->remapForRegion<void>(dataLocation, size);
	  void const *src = NULL;
	  switch (type) {
	  case SQLITE_TEXT:
	    src = sqlite3_value_text(value);
	    break;
	  case SQLITE_BLOB:
	    src = sqlite3_value_blob(value);
	    break;
	  default:
	    assert(false);
	    throw ASSERTION_ERROR;
	  }
	  memcpy(data, src, size);
	  header->tail = dataLocation + size;
	  entry->offset = dataLocation;
	  entry->size = size;
	}
      }
    }

    virtual void fetch(sqlite3_context *pCtx, sqlite_int64 rowId,
		       ResultSetter setter)
      throw(RTException) {
      enter();
      {
	LockHolder protect;
	LockHolder protectHeader;
	LockHolder protectData;
	map->takeLock(protect);
	headerMap->takeLock(protectHeader);
	dataMap->takeLock(protectData);

	assert(rowId > 0);
	off_t offset = (rowId - 1) * sizeof(Entry);
	Entry *entry = map->remapForRegion<Entry>(offset, sizeof(Entry));
	if (entry->offset == 0) { // NULL SQL value
	  sqlite3_result_null(pCtx);
	} else {
	  void *data = dataMap->remapForRegion<void>(entry->offset, entry->size);
	  setter(pCtx, data, entry->size);
	}
      }
    }
  };

  Column *createColumn(char const name[], SQLType type,
		       bool isPk, bool notNull) throw(RTException) {
    switch (type) {
    case SQLITE_TEXT:
    case SQLITE_BLOB:
      return new VariableSizeColumn(name, type, isPk, notNull);
    case SQLTYPE_INTEGER:
    case SQLTYPE_FLOAT:
      return new Column(name, type, isPk, notNull);
    }
    assert(false);
    throw ASSERTION_ERROR;
  }

  class Table: public virtual Openable, public virtual Closable  {
    string name_;
    string path_;
    vector<Column *> cols;
    sqlite_int64 nRows;
    int pkIdx;
    unique_ptr<File> tableDescFile;
    unique_ptr<MMap> tableDescMap;

    struct TableDescPage {
      char name[32];
      sqlite_int64 nRows;
      uint16_t nColumns;
      Column::ColumnDesc columns[1];
    };

    void selfOpen() throw(RTException) {
      enter();
      string tablePath = getTablePath();
      string tableDesc = strprintf("%s/%s", tablePath.c_str(), "_.desc");
      assert(tableDescFile.get() == NULL);
      tableDescFile.reset(new File(tableDesc.c_str(),
				   File::Mode_Read | File::Mode_Write));
      try {
	tableDescFile->open();
	try {
	  tableDescMap.reset(new MMap(tableDescFile.get()));
	  try {
	    tableDescMap->setMapRegion(0, MMap::getPageSize());
	    TableDescPage *desc = tableDescMap->map<TableDescPage>();
	    try {
	      if (strcmp(desc->name, name_.c_str()) == 0) {
		load();
	      }
	    } catch (...) {
	      tableDescMap->close();
	      throw;
	    }
	  } catch (...) {
	    tableDescMap.reset(NULL);
	    throw;
	  }
	} catch (...) {
	  tableDescFile->close();
	  throw;
	}
      } catch (...) {
	tableDescFile.reset(NULL);
	leave();
	throw;
      }
      leave();
    }
    void selfClose() throw(RTException) {
      if (tableDescMap.get() != NULL) {
	tableDescMap->close();
      }
      if (tableDescFile.get() != NULL) {
	tableDescFile->close();
      }
    }
    void closeCols(uint16_t startCol)
      throw(RTException) {
      if (startCol < cols.size()) {
	try {
	  cols[startCol]->close();
	} catch (...) {
	  try {
	    closeCols(startCol + 1);
	  } catch (...) {
	  }
	  throw;
	}
	closeCols(startCol + 1);
      }
    }
  public:
    Table(char const name[], char const path[]) throw(RTException)
      : name_(name), path_(path), cols(), nRows(0), pkIdx(-1),
	tableDescFile(), tableDescMap() {
      enter();
      if (strlen(name) >= sizeofMember(TableDescPage, name)) {
	static RTException ex("too long table name");
	throw ex;
      }
    }
    virtual ~Table() {
      enter();
      vector<Column *>::size_type end = cols.size();
      for (vector<Column *>::size_type i = 0; i < end; i++) {
	delete cols[i];
	cols[i] = NULL;
      }
    }

    string getTablePath() const {
      string tableDir = strprintf("%s/%s", path_.c_str(), name_.c_str());
      return tableDir;
    }

    virtual void open() throw(RTException) {
      enter();
      selfOpen();

      uint16_t closeUpTo = 0;
      try {
	for (uint16_t i = 0; i < cols.size(); ) {
	  cols[i]->open();
	  i++;
	  closeUpTo = i;
	}
      } catch (...) {
	for (uint16_t i = 0; i < closeUpTo; i++) {
	  try {
	    cols[i]->close();
	  } catch (...) {
	  }
	}
	leave();
	throw;
      }
      leave();
    }

    virtual void close() throw(RTException) {
      enter();
      try {
	closeCols(0);
      } catch (...) {
	selfClose();
	throw;
      }
      selfClose();
    }

    void load() throw(RTException) {
      enter();
      TableDescPage *desc = tableDescMap->getMappedAddr<TableDescPage>();
      nRows = desc->nRows;
      if (cols.size() != desc->nColumns) {
	static RTException ex("table mismatch");
	throw ex;
      }
    }
    void save() throw(RTException) {
      enter();
      TableDescPage *desc = tableDescMap->getMappedAddr<TableDescPage>();
      strcpy(desc->name, name_.c_str());
      desc->nRows = nRows;
      desc->nColumns = cols.size();
      for (uint16_t i = 0; i < desc->nColumns; i++) {
	cols[i]->save(&desc->columns[i]);
      }
      leave();
    }
    void create() throw(RTException) {
      enter();
      string tableDir = getTablePath();
      cout << "Parent dir of table: " << tableDir << endl;
      int result = mkdir(tableDir.c_str(), 0775);
      if (result != 0) {
	struct stat info;
	result = stat(tableDir.c_str(), &info);
	if (result != 0 || ! S_ISDIR(info.st_mode)) {
	  static RTException ex("mkdir error");
	  throw ex;
	}
      }
    }
    void addColumn(Column *col) throw(RTException) {
      assert(col != NULL);
      vector<Column *>::size_type idx = cols.size();
      if (&reinterpret_cast<TableDescPage*>(0)->columns[idx + 1] >
	  reinterpret_cast<void *>(MMap::getPageSize())) {
	static RTException ex("too much columns to store in a page.");
	throw ex;
      }
      if (col->isPK()) {
	if (pkIdx != -1) {
	  static RTException ex("duplicate primary key error");
	  throw ex;
	}
	pkIdx = (int)idx;
      }
      cols.push_back(col);
      col->table = this;
    }
    int getPKColumnIndex() const {
      return pkIdx;
    }
    Column *getColumn(int idx) const {
      return cols[idx];
    }
    string const &getName() const {
      return name_;
    }
    string const &getPath() const {
      return path_;
    }
    sqlite_int64 getNumberOfRows() const {
      return nRows;
    }
    void remove(sqlite_int64 rowId) throw(RTException) {
      enter();
      if (1 <= rowId && rowId == nRows) {
	// TODO delete
	nRows--;
      }
      static RTException ex("cannot delete any record other than last record.");
      throw ex;
    }
    sqlite_int64 insert(sqlite_int64 rowId,
			int argc, sqlite3_value **argv) throw(RTException) {
      enter();
      if (rowId < 0) {
	rowId = nRows + 1;
      }
      if (rowId != nRows + 1) {
	static RTException ex("cannot insert record with rowId which is not count(*) + 1.");
	throw ex;
      }
      for (int i = 0; i < argc; i++) {
	LOG cout << "Inserting column[" << i << "]\n";
	if (! cols[i]->isPK()) {
	  cols[i]->insert(rowId, argv[i]);
	}
      }
      nRows++;
      LOG cout << "Table::inserted" << nRows << endl;
      return rowId;
    }
    virtual void drop() throw(RTException) {
      enter();
      string tableDir = getTablePath();
      string command = strprintf("/bin/rm -rf %s", tableDir.c_str());
      int result = system(command.c_str());
      if (result != 0) {
	static RTException ex("failed to remove dir");
	throw ex;
      }
    }
  };

  File *Column::openFile(char const suffix[]) throw(RTException) {
    string dir = table->getTablePath();
    string filepath = strprintf("%s/%s%s", dir.c_str(), name_.c_str(), suffix);
    return new File(filepath.c_str(),
		    File::Mode_Read | File::Mode_Write);
  }

  class PredicateBase {
  public:
    virtual ~PredicateBase() {}
    //virtual bool isMatch(void const *lhsValue) = 0;
    /**
     * rowId starts with 1.
     */
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow) const = 0;
  };

  template <typename T>
  class Predicate: public PredicateBase {
  protected:
    T rValue;
  public:
    Predicate(T const &rhsValue)
      : rValue(rhsValue) {
    }
  };

  template <typename T>
  class PredEQ: public Predicate<T> {
  public:
    PredEQ(T const &rValue) : Predicate<T>(rValue) {}
    virtual bool isMatch(void const *lhsValue) {
      return *reinterpret_cast<T const *>(lhsValue) == this->rValue;
    }
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow)
      const {
      sqlite_int64 rValue = static_cast<sqlite_int64>(this->rValue);
      if (rValue <= nRow) {
	return rValue;
      }
      return nRow + 1; // not found
    }

  };

  template <typename T>
  class PredGT: public Predicate<T> {
  public:
    PredGT(T const &rValue) : Predicate<T>(rValue) {}
    virtual bool isMatch(void const *lhsValue) {
      return *reinterpret_cast<T const *>(lhsValue) > this->rValue;
    }
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow)
      const {
      sqlite_int64 rValue = static_cast<sqlite_int64>(this->rValue);
      if (rowId > rValue) {
	return rowId;
      }
      return rValue + 1;
    }
  };

  template <typename T>
  class PredGE: public Predicate<T> {
  public:
    PredGE(T const &rValue) : Predicate<T>(rValue) {}
    virtual bool isMatch(void const *lhsValue) {
      return *reinterpret_cast<T const *>(lhsValue) >= this->rValue;
    }
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow)
      const {
      sqlite_int64 rValue = static_cast<sqlite_int64>(this->rValue);
      if (rowId >= rValue) {
	return rowId;
      }
      return rValue;
    }
  };

  template <typename T>
  class PredLT: public Predicate<T> {
  public:
    PredLT(T const &rValue) : Predicate<T>(rValue) {}
    virtual bool isMatch(void const *lhsValue) {
      return *reinterpret_cast<T const *>(lhsValue) < this->rValue;
    }
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow)
      const {
      sqlite_int64 rValue = static_cast<sqlite_int64>(this->rValue);
      if (rowId < rValue) {
	return rowId;
      }
      return nRow + 1; // not found
    }
  };

  template <typename T>
  class PredLE: public Predicate<T> {
  public:
    PredLE(T const &rValue) : Predicate<T>(rValue) {}
    virtual bool isMatch(void const *lhsValue) {
      return *reinterpret_cast<T const *>(lhsValue) <= this->rValue;
    }
    virtual sqlite_int64 nextMatch(sqlite_int64 rowId, sqlite_int64 nRow)
      const {
      sqlite_int64 rValue = static_cast<sqlite_int64>(this->rValue);
      if (rowId <= rValue) {
	return rowId;
      }
      return nRow + 1; // not found
    }
  };

  class Cursor: public virtual Closable  {
    Mutex lock;
    enum State {
      NotSearchedYet,
      Scanning,
      EOFReached
    } state;
    Table *table_;
    sqlite_int64 rowId;
    vector<PredicateBase *> preds;
  public:
    Cursor(Table *table)
      : lock(), state(NotSearchedYet), table_(table), rowId(0), preds() {
      enter();
    }
    virtual ~Cursor() {
      enter();
      vector<PredicateBase *>::size_type end = preds.size();
      for (vector<PredicateBase *>::size_type i = 0; i < end; i++) {
	delete preds[i];
	preds[i] = NULL;
      }
    }

    Table *getTable() const {
      return table_;
    }

    /**
     * pred is owned by this instance.
     */
    void addPredicate(PredicateBase *pred) {
      assert(state == NotSearchedYet);
      preds.push_back(pred);
    }

    void find() {
      enter();
      assert(state == NotSearchedYet);
      assert(rowId == 0);
      state = Scanning;
      next();
    }

    void next() {
      enter();
      assert(state == Scanning);
      sqlite_int64 const nRows = table_->getNumberOfRows();
      //cout << "Cursor::next: nRows: " << nRows << endl;

      rowId++;
      //cout << "Cursor::next: rowId: " << rowId << endl;
      for (sqlite_int64 lastRowId = rowId; rowId <= nRows; lastRowId = rowId) {
	vector<PredicateBase *>::size_type end = preds.size();
	//cout << "Cursor::next: preds: " << end << endl;
	for (vector<PredicateBase *>::size_type i = 0; i < end; i++) {
	  rowId = preds[i]->nextMatch(rowId, nRows);
	}
	if (lastRowId == rowId) { // matches all preds or not matches at all.
	  break;
	}
      }
      if (rowId > nRows) {
	//cout << "Cursor::next: eof reached: " << endl;
	state = EOFReached;
      }
    }

    bool isEOF() const {
      assert(state == Scanning || state == EOFReached);
      return state == EOFReached;
    }

    sqlite_int64 getRowId() const {
      if (state == Scanning) {
	return rowId;
      }
      static RTException ex("attempting to get RowId against a cursor which doesn't point a valid record.");
      throw ex;
    }

    virtual void close() throw(RTException) {
      enter();
      // TODO implement here
    }
  };

#if 0


template<typename T>
T getParam(T (*func)(sqlite3_value *), sqlite3_value *value, T valueIfNull) {
  if (sqlite3_value_type(value) == SQLITE_NULL) {
    return valueIfNull;
  }
  return func(value);
}

extern "C" {
}

#endif

  template <typename E>
  class Scanner {
    E const *tokenId_;
    char const * const *tokenStr_;
    size_t const num_;
    int (*cmp_)(const char *, const char *, size_t);
    char const *separators_;
  public:
    Scanner(E const tokenId[], char const *tokenStr[], size_t num,
	    int (*cmp)(const char *, const char *, size_t),
	    char const separators[])
      : tokenId_(tokenId), tokenStr_(tokenStr), num_(num), cmp_(cmp),
	separators_(separators) { }
    virtual ~Scanner() {}

    void init(char const str[], char const **ctx) const {
      assert(ctx != NULL);
      *ctx = str;
    }
    E nextToken(char const **ctx, E tokenOfEOF, E invalidToken,
		char const **invalidTokenStr, size_t *invalidTokenStrLen) const throw(RTException) {
      assert(ctx != NULL);
      assert(invalidTokenStr != NULL);
      assert(invalidTokenStrLen != NULL);
      *invalidTokenStr = NULL;
      *invalidTokenStrLen = 0;
      char const *p = *ctx;
      while (*p && strchr(separators_, *p) != NULL) { // skip seps
	p++;
      }
      char const *tokStart = p;
      while (*p && strchr(separators_, *p) == NULL) { // find end of token.
	p++;
      }
      *ctx = p;
      size_t len = p - tokStart;
      if (len == 0) {
	return tokenOfEOF;
      }
      for (size_t i = 0; i < num_; i++) {
	if (cmp_(tokStart, tokenStr_[i], len) == 0) {
	  return tokenId_[i];
	}
      }
      *invalidTokenStr = tokStart;
      *invalidTokenStrLen = len;
      return invalidToken;
    }
  };

  class NewTableParser {
    enum Token {
      T_int,
      T_text,
      T_blob,
      T_none,
      T_real,
      T_primary,
      T_key,
      T_not,
      T_null,
      T_NUM,
      T_id,
      T_EOF
    };
    static Token const tokens[T_NUM];
    static char const *tokenStrs[T_NUM];
    static Scanner<Token> const scanner;

    NewTableParser() {
      assert(false);
      throw ASSERTION_ERROR;
    }
  public:
    static SQLType parse(char const *str,
			 char const **name, size_t *nameLen,
			 bool *isPK, bool *notNull) throw(RTException) {
      assert(str != NULL);
      assert(name != NULL);
      assert(nameLen != NULL);
      assert(isPK != NULL);
      assert(notNull != NULL);
      *name = NULL;
      *nameLen = 0;
      *isPK = false;
      *notNull = false;

      static RTException ex("parse error");

      char const *ctx = NULL;
      scanner.init(str, &ctx);

      char const *colName = NULL;
      size_t colNameLen = 0;
      {
	Token t = scanner.nextToken(&ctx, T_EOF, T_id, &colName, &colNameLen);
	if (t != T_id) {
	  throw ex;
	}
      }

      SQLType type = SQLTYPE_NULL;
      {
	char const *invTok = NULL;
	size_t invTokLen = 0;
	Token t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	switch (t) {
	case T_int:
	  type = SQLTYPE_INTEGER;
	  break;
	case T_text:
	  type = SQLTYPE_TEXT;
	  break;
	case T_blob:
	case T_none:
	  type = SQLTYPE_BLOB;
	  break;
	case T_real:
	  type = SQLTYPE_FLOAT;
	  break;
	default:
	  assert(false);
	  throw ASSERTION_ERROR;
	}
      }

      do {
	char const *invTok = NULL;
	size_t invTokLen = 0;
	Token t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	if (t == T_EOF) {
	  *isPK = false;
	  break;
	}
	if (t == T_primary) {
	  t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	  if (t == T_key) {
	    t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	    if (t == T_EOF) {
	      *isPK = true;
	      *notNull = true;
	      break;
	    }
	  }
	} else if (t == T_not) {
	  t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	  if (t == T_null) {
	    t = scanner.nextToken(&ctx, T_EOF, T_id, &invTok, &invTokLen);
	    if (t == T_EOF) {
	      *notNull = true;
	      break;
	    }
	  }
	}
	throw ex;
      } while (true);
      *name = colName;
      *nameLen = colNameLen;
      return type;
    }
  };
  NewTableParser::Token const NewTableParser::tokens[] = {
    NewTableParser::T_int,
    NewTableParser::T_text,
    NewTableParser::T_blob,
    NewTableParser::T_none,
    NewTableParser::T_real,
    NewTableParser::T_primary,
    NewTableParser::T_key,
    NewTableParser::T_not,
    NewTableParser::T_null
  };
  char const *NewTableParser::tokenStrs[] = {
      "integer",
      "text",
      "blob",
      "none",
      "real",
      "primary",
      "key",
      "not",
      "null"
    };
  Scanner<NewTableParser::Token> const
  NewTableParser::scanner(NewTableParser::tokens,
			  NewTableParser::tokenStrs,
			  elementsof(NewTableParser::tokens),
			  strncasecmp, " \t\n\r");

  struct MMapVTab {
    sqlite3_vtab base;
    Table *table;
  };
  struct MMapCursor {
    sqlite3_vtab_cursor base;
    Cursor *cursor;
  };

  int mmapCreateOrConnectMethod(sqlite3 *db, void *pAux,
				int argc, char const*const*argv,
				sqlite3_vtab **ppVTab,
				char **pzErr,
				bool doCreate) throw() {
    if (argc < 5) {
      *pzErr = sqlite3_mprintf("Too less arguments for a virtual table.");
      return SQLITE_ERROR;
    }
    char const *modName = argv[0];
    char const *dbName = argv[1];
    char const *vtName = argv[2];
    string sqlStr(argv[3]);
    string::size_type start = 0;
    string::size_type len =  sqlStr.size();
    if (len > 0 && sqlStr[0] == '\'') {
      start = 1;
      len--;
    }
    if (len > 0 && sqlStr[start + len - 1] == '\'') {
      len--;
    }
    sqlStr = sqlStr.substr(start, len);
    char const *vtPath = sqlStr.c_str();
    unique_ptr<Table> vt(new Table(vtName, vtPath));
    string ddl = "create table ";
    ddl += vtName;
    char const *sep = "(";
    try {
      for (int i = 4; i < argc; i++) {
	char const *colName = NULL;
	size_t colNameLen = 0;
	bool isPK = false;
	bool notNull = false;
	SQLType type =
	  NewTableParser::parse(argv[i], &colName, &colNameLen, &isPK, &notNull);
	string colname(colName, colNameLen);
	unique_ptr<Column> col(createColumn(colname.c_str(), type, isPK, notNull));
	vt->addColumn(col.get());
	col.release();
	ddl += sep;
	ddl += colname + " " + toSQLTypeName(type);
	sep = ", ";
      }
      ddl += ")";
    } catch (...) {
      *pzErr = sqlite3_mprintf("Failed to create virtual table.");
      return SQLITE_ERROR;
    }
    LOG cout << ddl << endl;

    try {
      if (doCreate) {
	vt->create();
      }
    } catch (...) {
      *pzErr = sqlite3_mprintf("Failed to create virtual table on disk.");
      return SQLITE_ERROR;
    }
    try {
      vt->open();
      try {
	MMapVTab *mmapVtab = (MMapVTab *)sqlite3_malloc(sizeof(MMapVTab));
	if (mmapVtab == NULL) {
	  vt->close();
	  return SQLITE_NOMEM;
	}
	mmapVtab->base.zErrMsg = NULL;
	mmapVtab->table = vt.get();
	int result = sqlite3_declare_vtab(db, ddl.c_str());
	if (result != SQLITE_OK) {
	  sqlite3_free(mmapVtab);
	  *pzErr = sqlite3_mprintf("sqlite3_declare_vtab failed.");
	  vt->close();
	  return SQLITE_ERROR;
	}
	vt.release();
	*ppVTab = reinterpret_cast<sqlite3_vtab *>(mmapVtab);
	return SQLITE_OK;
      } catch (...) {
	vt->close();
	throw;
      }
    } catch (...) {
    }
    *pzErr = sqlite3_mprintf("failed to open table.");
    return SQLITE_ERROR;
  }
  int mmapCreateMethod(sqlite3 *db, void *pAux,
		       int argc, char const*const*argv,
		       sqlite3_vtab **ppVTab,
		       char **pzErr) throw() {
    enter();
    return mmapCreateOrConnectMethod(db, pAux, argc, argv, ppVTab,
				     pzErr, true);
  }

  int mmapConnectMethod(sqlite3 *db, void *pAux,
		       int argc, char const*const*argv,
		       sqlite3_vtab **ppVTab,
		       char **pzErr) {
    enter();
    return mmapCreateOrConnectMethod(db, pAux, argc, argv, ppVTab,
				     pzErr, false);
  }

  int mmapDisconnectMethod(sqlite3_vtab *pVTab) throw() {
    enter();
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    int result = SQLITE_OK;
    try {
      tab->save();
    } catch (...) {
      result = SQLITE_ERROR;
    }
    try {
      tab->close();
    } catch (...) {
      result = SQLITE_ERROR;
    }
    delete tab;
    sqlite3_free(mmapVtab);
    return result;
  }

  int mmapDestroyMethod(sqlite3_vtab *pVTab) throw() {
    enter();
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    int result = SQLITE_OK;
    try {
      tab->close();
    } catch (...) {
      result = SQLITE_ERROR;
    }
    try {
      tab->drop();
    } catch (...) {
      result = SQLITE_ERROR;
    }
    delete tab;
    sqlite3_free(mmapVtab);
    return result;
  }

  int mmapOpenMethod(sqlite3_vtab *pVTab, sqlite3_vtab_cursor **ppCursor)
    throw() {
    enter();
    assert(pVTab != NULL);
    assert(ppCursor != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;

    unique_ptr<Cursor> cursor(new Cursor(tab));
    MMapCursor *mmapCursor = (MMapCursor *)sqlite3_malloc(sizeof(MMapCursor));
    mmapCursor->cursor = cursor.get();
    cursor.release();
    *ppCursor = reinterpret_cast<sqlite3_vtab_cursor *>(mmapCursor);
    return SQLITE_OK;
  }

  int mmapCloseMethod(sqlite3_vtab_cursor *pCursor) throw() {
    enter();
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    cursor->close();
    delete cursor;
    sqlite3_free(mmapCursor);
    return SQLITE_OK;
  }

  int mmapBestIndexMethod(sqlite3_vtab *pVTab, sqlite3_index_info *pIdxInfo)
    throw() {
    enter();
    assert(pVTab != NULL);
    assert(pIdxInfo != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    char idxStr[pIdxInfo->nConstraint * 2 + 1];
    int idx = 0;

    LOG cout << "best index\n";
    LOG cout << "best index: constraint: " << pIdxInfo->nConstraint << "\n";
    int pkCol = tab->getPKColumnIndex(); // may be -1
    for (int i = 0; i < pIdxInfo->nConstraint; i++) {
      if (pIdxInfo->aConstraint[i].usable) {
	if (pIdxInfo->aConstraint[i].iColumn != pkCol
	    && pIdxInfo->aConstraint[i].iColumn != -1) {
	  sqlite3_free(pVTab->zErrMsg);
	  pVTab->zErrMsg = sqlite3_mprintf("Any column other than PK can't be a constraint.");
	  leave();
	  return SQLITE_ERROR;
	}
	switch (pIdxInfo->aConstraint[i].op) {
	case SQLITE_INDEX_CONSTRAINT_MATCH:
	  sqlite3_free(pVTab->zErrMsg);
	  pVTab->zErrMsg = sqlite3_mprintf("LIKE operater is not allowed for this column.");
	  leave();
	  return SQLITE_ERROR;
	}
	assert(pIdxInfo->aConstraintUsage != NULL);
	idxStr[idx++] = pIdxInfo->aConstraint[i].op;
	idxStr[idx++] = pIdxInfo->aConstraint[i].iColumn + 1 + 'A';
	pIdxInfo->aConstraintUsage[i].argvIndex = idx / 2;
	pIdxInfo->aConstraintUsage[i].omit = 0; // set to 1 after debug
      }
    }
    LOG cout << "best index loop end\n";
    pIdxInfo->idxNum = idx;
    idxStr[idx++] = '\0';
    assert(idx <= pIdxInfo->nConstraint * 2 + 1);
    pIdxInfo->needToFreeIdxStr = 1;
    pIdxInfo->idxStr = sqlite3_mprintf("%s", idxStr);
    if (pIdxInfo->idxStr == NULL) {
      leave();
      return SQLITE_NOMEM;
    }
    assert(! pIdxInfo->orderByConsumed);
    pIdxInfo->estimatedCost = 200000.0;
    leave();
    return SQLITE_OK;
  }

  int mmapFilterMethod(sqlite3_vtab_cursor *pCursor,
		       int idxNum, char const*idxStr,
		       int argc, sqlite3_value **argv) throw() {
    enter();
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    Table *table = cursor->getTable();
    const int nConstraint = idxNum / 2;
    assert(argc == nConstraint);
    for (int i = 0; i < nConstraint; i++) {
      int op = idxStr[i*2];
      int iColumn = idxStr[i*2 + 1] - 'A' - 1;
      sqlite3_value *value = argv[i];
      // TODO set to Cursor
      Column *col = NULL;
      if (iColumn >= 0) {
	col = table->getColumn(iColumn);
      }
      if (iColumn < 0
	  || (col && col->isPK() && col->getType() == SQLTYPE_INTEGER)) {
	// RowId
	switch (op) {
	case SQLITE_INDEX_CONSTRAINT_EQ:
	  {
	    PredicateBase *pred =
	      new PredEQ<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_GT:
	  {
	    PredicateBase *pred =
	      new PredGT<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_GE:
	  {
	    PredicateBase *pred =
	      new PredGE<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_LT:
	  {
	    PredicateBase *pred =
	      new PredLT<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_LE:
	  {
	    PredicateBase *pred =
	      new PredLE<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;

	case SQLITE_INDEX_CONSTRAINT_MATCH:
	default:
	  assert(false);
	  throw ASSERTION_ERROR;
	}
      } else {
#if 0
	switch (op) {
	case SQLITE_INDEX_CONSTRAINT_EQ:
	  {
	    PredicateBase *pred =
	      new PredEQ<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_GT:
	  {
	    PredicateBase *pred =
	      new PredGT<sqlite_int64>(static_cast<sqlite_int64>(sqlite3_value_int(argv[i])));
	    cursor->addPredicate(pred);
	  }
	  break;
	case SQLITE_INDEX_CONSTRAINT_LE:
	case SQLITE_INDEX_CONSTRAINT_LT:
	case SQLITE_INDEX_CONSTRAINT_GE:
	  break;
	case SQLITE_INDEX_CONSTRAINT_MATCH:
	default:
	  assert(false);
	  throw ASSERTION_ERROR;
	}
#else
	assert(false);
	throw ASSERTION_ERROR;
#endif
      }
    }
    cursor->find();
    return SQLITE_OK;
  }

  int mmapEofMethod(sqlite3_vtab_cursor *pCursor) throw() {
    enter();
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    return cursor->isEOF();
  }

  int mmapNextMethod(sqlite3_vtab_cursor *pCursor) throw() {
    enter();
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    cursor->next();
    return SQLITE_OK;
  }

  int mmapRowidMethod(sqlite3_vtab_cursor *pCursor,
		      sqlite_int64 *pRowid) throw() {
    enter();
    assert(pCursor != NULL);
    assert(pRowid != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    sqlite_int64 rowId = -1;
    try {
      rowId = cursor->getRowId();
    } catch (...) {
      return SQLITE_ERROR;
    }
    *pRowid = rowId;
    return SQLITE_OK;
  }

  int mmapColumnMethod(sqlite3_vtab_cursor *pCursor, sqlite3_context *pCtx,
		       int pos) throw() {
    enter();
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    //cout << "column: " << pos << "\n";
    if (false) { // error
      sqlite3_result_error_code(pCtx, SQLITE_IOERR);
      return SQLITE_IOERR;
    }
    try {
      Column *col = cursor->getTable()->getColumn(pos);
      if (col->isPK()) { // PK should be rowid because of restriction.
	sqlite3_result_int64(pCtx, cursor->getRowId());
      } else {
	col->fetch(pCtx, cursor->getRowId());
      }
    } catch (...) {
      return SQLITE_ERROR;
    }
    return SQLITE_OK;
  }

  int mmapUpdateMethod(sqlite3_vtab *pVTab,
		       int argc, sqlite3_value **argv,
		       sqlite_int64 *pRowid) throw() {
    // TODO lock table
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    int i = 0;
    assert(argc > 0);
    LOG cout << "update: argc: " << argc << endl;
    if (sqlite3_value_type(argv[i]) != SQLITE_NULL) {
      if (sqlite3_value_type(argv[i]) == SQLITE_INTEGER) {
	sqlite_int64 rowId = sqlite3_value_int(argv[i]);
	try {
	  tab->remove(rowId);
	} catch (...) {
	  return SQLITE_ERROR;
	}
      } else {
	assert(false);
	throw ASSERTION_ERROR;
      }
    }
    if (argc > 1) {
      i++;
      sqlite_int64 rowId = -1;
      if (sqlite3_value_type(argv[i]) != SQLITE_NULL) {
	if (sqlite3_value_type(argv[i]) == SQLITE_INTEGER) {
	  rowId = sqlite3_value_int(argv[i]);
	} else {
	  assert(false);
	  throw ASSERTION_ERROR;
	}
      }
      try {
	rowId = tab->insert(rowId, argc - 2, &argv[2]);
	*pRowid = rowId;
      } catch (...) {
	return SQLITE_ERROR;
      }
    } else {
      assert(false);
      throw ASSERTION_ERROR;
    }
    return SQLITE_OK;
  }

  int mmapRenameMethod(sqlite3_vtab *pVTab, const char *zNew) throw() {
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    // TODO support rename
    return SQLITE_ERROR;
  }

  static sqlite3_module const mmapModule = {
    /* iVersion      */ 1,
    /* xCreate       */ mmapCreateMethod,
    /* xConnect      */ mmapConnectMethod,
    /* xBestIndex    */ mmapBestIndexMethod,
    /* xDisconnect   */ mmapDisconnectMethod,
    /* xDestroy      */ mmapDestroyMethod,
    /* xOpen         */ mmapOpenMethod,
    /* xClose        */ mmapCloseMethod,
    /* xFilter       */ mmapFilterMethod,
    /* xNext         */ mmapNextMethod,
    /* xEof          */ mmapEofMethod,
    /* xColumn       */ mmapColumnMethod,
    /* xRowid        */ mmapRowidMethod,
    /* xUpdate       */ mmapUpdateMethod,
    /* xBegin        */ NULL, // mmapBeginMethod,
    /* xSync         */ NULL, // mmapSyncMethod,
    /* xCommit       */ NULL, // mmapCommitMethod,
    /* xRollback     */ NULL, // mmapRollbackMethod,
    /* xFindFunction */ NULL, // mmapFindFunctionMethod,
    /* xRename */       mmapRenameMethod,
    /* xSavepoint    */ NULL, // mmapSavepointMethod,
    /* xRelease      */ NULL, // mmapReleaseMethod,
    /* xRollbackTo   */ NULL, // mmapRollbackToMethod,
  };

} // namespace

using namespace a;
extern "C" {
  int sqlite3_extension_init(sqlite3 *db,
			     char **pzErrMsg,
			     const sqlite3_api_routines *pApi) {
    SQLITE_EXTENSION_INIT2(pApi)
      ;
    init_module();
    void *pClientData = NULL;         /* Client data for xCreate/xConnect */
    void(*xDestroy)(void*) = NULL;     /* Module destructor function */
    int result = sqlite3_create_module_v2(db,
					  "MMapVTable",
					  &mmapModule, pClientData, xDestroy);
    if (result != 0) {
      *pzErrMsg = sqlite3_mprintf("Can't create MMapVTable module.");
      return result;
    }
    return SQLITE_OK;
  }
}
