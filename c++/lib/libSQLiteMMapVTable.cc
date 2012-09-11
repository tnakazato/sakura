#include <new>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <cassert>
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
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define PERROR(x) perror(__FILE__  ":" #x)

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
    Closer(Closer const &other) : closable_(NULL) { assert(false); }
  public:
    Closer(Closable *closable) : closable_(closable) {
      assert(closable_ != NULL);
    }
    virtual ~Closer() throw(RTException) { closable_->close(); }
  };

  class Mutex: public virtual Closable {
    pthread_mutex_t mutex;
    Mutex(Mutex const &other) { assert(false); }
    Mutex &operator =(Mutex const &other) { assert(false); }
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
      int result = pthread_mutex_unlock(&mutex);
      if (result != 0) {
	static RTException ex("pthread_mutex_unlock error");
	throw ex;
      }
    }

    virtual void close() throw(RTException) {
      unlock();
    }
  };

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
    char const *filename_;
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
      : filename_(NULL), mode_(mode), omode(0), fd(-1), refCount(0) {
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
      filename_ = StrDup(filename);
    }
    virtual ~File() throw(RTException) {
      // without lock there should be no referer.
      if (fd >= 0) {
	close();
      }
      delete[] filename_;
      assert(refCount == 0);
    }
    int getMode() const {
      return mode_;
    }
    virtual void open() throw(RTException) {
      lock.lock();
      Closer autoUnlock(&lock);
      if (fd == -1) {
	fd = ::open(filename_, omode);
	if (fd < 0) {
	  PERROR(open);
	  static RTException ex("open error");
	  throw ex;
	}
      }
    }
    virtual void close() throw(RTException) {
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
  public:
    MMap(File *file) throw(RTException)
      : file_(file), mappedAddr(MAP_FAILED), size_(0), fileOffset_(0) {
      assert(file_ != NULL);
      file_->ref();
    }
    virtual ~MMap() throw(RTException) {
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
      assert(offset >= 0 && offset % getPageSize() == 0);
      lock.lock();
      Closer autoUnlock(&lock);
      if (mappedAddr != MAP_FAILED) {
	static RTException ex("unmap before changing map region.");
	throw ex;
      }
      fileOffset_ = offset;
      size_ = size;
    }
      
    void map() throw(RTException) {
      lock.lock();
      Closer autoUnlock(&lock);
      if (mappedAddr == MAP_FAILED) {
	mappedAddr = mmap(NULL, size_, file_->getMode(),
			  MAP_SHARED, file_->fd, fileOffset_);
	if (mappedAddr == MAP_FAILED) {
	  PERROR(mmap);
	  static RTException ex("mmap error");
	  throw ex;
	}
      }
    }
    void unmap() throw(RTException) {
      lock.lock();
      Closer autoUnlock(&lock);
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
      map();
    }
    virtual void close() throw(RTException) {
      unmap();
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
  class Column {
    Table *table;
    char const *name_;
    SQLType type_;
    size_t size;
    bool isPK_;
    bool notNull_;
    friend class Table;
  public:
    Column(char const name[], SQLType type, bool isPk, bool notNull)
      : table(NULL), name_(NULL), type_(type), isPK_(isPk), notNull_(notNull) {
      size = 0;
      switch (type_) {
      case SQLTYPE_INTEGER:
	size = sizeof(int64_t);
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
	break;
      }
      name_ = StrDup(name);
    }
    virtual ~Column() {
      delete[] name_;
    }
    char const *getName() const {
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

  class Table: public virtual Closable  {
    char const *name_;
    char const *path_;
    vector<Column *> cols;
    sqlite_int64 nRows;
    int pkIdx;
  public:
    Table(char const name[], char const path[])
      : name_(NULL), path_(NULL), cols(), nRows(0), pkIdx(-1) {
      name_ = StrDup(name);
      try {
	path_ = StrDup(path);
      } catch (...) {
	delete[] name_;
      }
    }
    virtual ~Table() {
      delete[] path_;
      delete[] name_;
      vector<Column *>::size_type end = cols.size();
      for (vector<Column *>::size_type i = 0; i < end; i++) {
	delete cols[i];
	cols[i] = NULL;
      }
    }
    void addColumn(Column *col) throw(RTException) {
      assert(col != NULL);
      vector<Column *>::size_type idx = cols.size();
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
    char const *getName() const {
      return name_;
    }
    char const *getPath() const {
      return path_;
    }
    sqlite_int64 getNumberOfRows() const {
      return nRows;
    }
    void remove(sqlite_int64 rowId) throw(RTException) {
      if (1 <= rowId && rowId == nRows) {
	// TODO delete
	nRows--;
      }
      static RTException ex("cannot delete any record other than last record.");
      throw ex;
    }
    sqlite_int64 insert(sqlite_int64 rowId,
			int argc, sqlite3_value **argv) throw(RTException) {
      if (rowId < 0) {
	rowId = nRows + 1;
      }
      if (rowId != nRows + 1) {
	static RTException ex("cannot insert record with rowId which is not count(*) + 1.");
	throw ex;
      }
      // TODO insert
      nRows++;
      cout << "Table::inserted" << nRows << endl;
      return rowId;
    }
    virtual void close() throw(RTException) {
      // TODO implement here
    }
    virtual void drop() throw(RTException) {
      // TODO implement here
    }
  };

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
    }
    virtual ~Cursor() {
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
      assert(state == NotSearchedYet);
      assert(rowId == 0);
      state = Scanning;
      next();
    }

    void next() {
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

    NewTableParser() { assert(false); }
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
	  break;
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

  int mmapCreateMethod(sqlite3 *db, void *pAux,
		       int argc, char const*const*argv,
		       sqlite3_vtab **ppVTab,
		       char **pzErr) throw() {
    if (argc < 5) {
      *pzErr = sqlite3_mprintf("Too less arguments for a virtual table.");
      return SQLITE_ERROR;
    }
    char const *modName = argv[0];
    char const *dbName = argv[1];
    char const *vtName = argv[2];
    char const *vtPath = argv[3];
    auto_ptr<Table> vt(new Table(vtName, vtPath));
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
	auto_ptr<Column> col(new Column(colname.c_str(), type, isPK, notNull));
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
    cout << ddl << endl;
    int result = sqlite3_declare_vtab(db, ddl.c_str());
    if (result != SQLITE_OK) {
      *pzErr = sqlite3_mprintf("sqlite3_declare_vtab failed.");
      return SQLITE_ERROR;
    }
    MMapVTab *mmapVtab = (MMapVTab *)sqlite3_malloc(sizeof(MMapVTab));
    mmapVtab->base.zErrMsg = NULL;
    mmapVtab->table = vt.get();
    vt.release();
    *ppVTab = reinterpret_cast<sqlite3_vtab *>(mmapVtab);
    return SQLITE_OK;
  }

  int mmapConnectMethod(sqlite3 *db, void *pAux,
		       int argc, char const*const*argv,
		       sqlite3_vtab **ppVTab,
		       char **pzErr) {
    return mmapCreateMethod(db, pAux, argc, argv, ppVTab, pzErr);
  }

  int mmapDisconnectMethod(sqlite3_vtab *pVTab) throw() {
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    tab->close();
    delete tab;
    sqlite3_free(mmapVtab);
    return SQLITE_OK;
  }

  int mmapDestroyMethod(sqlite3_vtab *pVTab) throw() {
    assert(pVTab != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    tab->close();
    tab->drop();
    delete tab;
    sqlite3_free(mmapVtab);
    return SQLITE_OK;
  }

  int mmapOpenMethod(sqlite3_vtab *pVTab, sqlite3_vtab_cursor **ppCursor)
    throw() {
    assert(pVTab != NULL);
    assert(ppCursor != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;

    auto_ptr<Cursor> cursor(new Cursor(tab));
    MMapCursor *mmapCursor = (MMapCursor *)sqlite3_malloc(sizeof(MMapCursor));
    mmapCursor->cursor = cursor.get();
    cursor.release();
    *ppCursor = reinterpret_cast<sqlite3_vtab_cursor *>(mmapCursor);
    return SQLITE_OK;
  }

  int mmapCloseMethod(sqlite3_vtab_cursor *pCursor) throw() {
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
    assert(pVTab != NULL);
    assert(pIdxInfo != NULL);
    MMapVTab *mmapVtab = reinterpret_cast<MMapVTab *>(pVTab);
    Table *tab = mmapVtab->table;
    char idxStr[pIdxInfo->nConstraint * 2 + 1];
    int idx = 0;

    cout << "best index\n";
    cout << "best index: constraint: " << pIdxInfo->nConstraint << "\n";
    int pkCol = tab->getPKColumnIndex(); // may be -1
    for (int i = 0; i < pIdxInfo->nConstraint; i++) {
      if (pIdxInfo->aConstraint[i].usable) {
	if (pIdxInfo->aConstraint[i].iColumn != pkCol
	    && pIdxInfo->aConstraint[i].iColumn != -1) {
	  sqlite3_free(pVTab->zErrMsg);
	  pVTab->zErrMsg = sqlite3_mprintf("Any column other than PK can't be a constraint.");
	  return SQLITE_ERROR;
	}
	switch (pIdxInfo->aConstraint[i].op) {
	case SQLITE_INDEX_CONSTRAINT_MATCH:
	  sqlite3_free(pVTab->zErrMsg);
	  pVTab->zErrMsg = sqlite3_mprintf("LIKE operater is not allowed for this column.");
	  return SQLITE_ERROR;
	}
	assert(pIdxInfo->aConstraintUsage != NULL);
	idxStr[idx++] = pIdxInfo->aConstraint[i].op;
	idxStr[idx++] = pIdxInfo->aConstraint[i].iColumn + 1 + 'A';
	pIdxInfo->aConstraintUsage[i].argvIndex = idx / 2;
	pIdxInfo->aConstraintUsage[i].omit = 0; // set to 1 after debug
      }
    }
    cout << "best index loop end\n";
    pIdxInfo->idxNum = idx;
    idxStr[idx++] = '\0';
    assert(idx <= pIdxInfo->nConstraint * 2 + 1);
    pIdxInfo->needToFreeIdxStr = 1;
    pIdxInfo->idxStr = sqlite3_mprintf("%s", idxStr);
    if (pIdxInfo->idxStr == NULL) {
      return SQLITE_NOMEM;
    }
    assert(! pIdxInfo->orderByConsumed);
    pIdxInfo->estimatedCost = 200000.0;
    cout << "best index exit\n";
    return SQLITE_OK;
  }

  int mmapFilterMethod(sqlite3_vtab_cursor *pCursor,
		       int idxNum, char const*idxStr,
		       int argc, sqlite3_value **argv) throw() {
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    Table *table = cursor->getTable();
    const int nConstraint = idxNum / 2;
    assert(argc == nConstraint);
    cout << "filter\n";
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
	  break;
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
	  break;
	}
#else
	assert(false);
#endif
      }
    }
    cursor->find();
    return SQLITE_OK;
  }

  int mmapEofMethod(sqlite3_vtab_cursor *pCursor) throw() {
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    //cout << "eof\n";
    return cursor->isEOF();
  }

  int mmapNextMethod(sqlite3_vtab_cursor *pCursor) throw() {
    assert(pCursor != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    cursor->next();
    //cout << "next\n";
    return SQLITE_OK;
  }

  int mmapRowidMethod(sqlite3_vtab_cursor *pCursor,
		      sqlite_int64 *pRowid) throw() {
    assert(pCursor != NULL);
    assert(pRowid != NULL);
    MMapCursor *mmapCursor = reinterpret_cast<MMapCursor *>(pCursor);
    Cursor *cursor = mmapCursor->cursor;
    sqlite_int64 rowId = -1;
    //cout << "rowid\n";
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
	// TODO implement here
	//sqlite3_result_blob(pCtx, "", 0, SQLITE_TRANSIENT);
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
    cout << "update: argc: " << argc << endl;
    if (sqlite3_value_type(argv[i]) != SQLITE_NULL) {
      if (sqlite3_value_type(argv[i]) == SQLITE_INTEGER) {
	sqlite_int64 rowId = sqlite3_value_int(argv[i]);
	try {
	  tab->remove(rowId);
	} catch (...) {
      cout << "update: error del" << endl;
	  return SQLITE_ERROR;
	}
      } else {
	assert(false);
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
	}
      }
      try {
	rowId = tab->insert(rowId, argc - 2, &argv[2]);
	*pRowid = rowId;
      } catch (...) {
	cout << "update: error" << endl;
	return SQLITE_ERROR;
      }
    } else {
      assert(false);
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
