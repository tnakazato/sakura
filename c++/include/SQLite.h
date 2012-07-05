#ifndef SQLITE_H_
#define SQLITE_H_

#include <string>
#include <sqlite3.h>

#ifdef DEBUG
#define SQLDEBUG 1
#endif

/**
 * C++ wrapper library for SQLite.
 * All strings handled by this library should be UTF-8.
 *
 * Lifetime of following objects should be:
 *
 * ResultSet within PreparedStatement within Connection.
 * @author Kohji Nakamura
 * @version $Revision$
 */
namespace sqlite {

/**
 * This exception will be raised when there is an error while using this API.
 */
class SQLException {
 public:
  virtual ~SQLException();
};

/**
 * This is an abstract base class for all statement classes.
 */
class Statement {
 public:
  virtual ~Statement();
};

class PreparedStatement;
class Connection;

/**
 * This class represents a result of a query.
 * You need to call {@link #next()} before fetching the first row
 * since this class initially points a pseudo row
 * just before the first result row.
 *
 */
class ResultSet {
  PreparedStatement *const stmt;

  ResultSet(PreparedStatement *stmt);
  friend class PreparedStatement;
 public:
  virtual ~ResultSet() throw (SQLException);
  /**
   * Go to the next row if it exists.
   * @return true if there was a next row, otherwise false.
   */
  virtual bool next() throw (SQLException);
  /**
   * Returns a number of columns within the result.
   * @return a number of columns.
   */
  virtual int getColumnCount() throw (SQLException);
  /**
   * Get a column name for the `pos'-th column.
   * @param pos A position of the column. It starts with 1.
   * @return name of the column.
   */
  virtual std::string getColumnName(int pos) throw (SQLException);
  /**
   * Returns whether the specified column is null or not.
   * This method must be called before any invocations of
   * get&lt;Type&gt;() methods for the same column for each row.
   * Generally, you should call this method before calling get&lt;Type&gt;()
   * methods for the same column when accessing nullable column.
   * @param pos A position of the column. It starts with 1.
   * @return true if the value of the column is null, otherwise false.
   */
  virtual bool isNull(int pos) throw (SQLException);
  /**
   * Fetches an int value from the specified column.
   * @param pos A position of the column. It starts with 1.
   * @return an int value.
   */
  virtual int getInt(int pos) throw (SQLException);
  /**
   * Fetches a string value from the specified column.
   * @param pos A position of the column. It starts with 1.
   * @param size A pointer to an integer in which size of the returned value
   * (including trailing '\\0') will be stored.
   * @return '\\0' terminated UTF-8 string value. The returned value is valid
   * until calling {@link #next()} or get&lt;Type&gt;() methods for the same column
   * or destructing the instance of this class. You must copy the value
   * if you need to access it beyond such cases.
   * You must not release the returned area.
   */
  virtual char const *getTransientString(int pos, int *size) throw (SQLException);
  /**
   * Fetches a float value from the specified column.
   * @param pos A position of the column. It starts with 1.
   * @return a float value.
   */
  virtual float getFloat(int pos) throw (SQLException);
  /**
   * Fetches a double value from the specified column.
   * @param pos A position of the column. It starts with 1.
   * @return a double value.
   */
  virtual double getDouble(int pos) throw (SQLException);
  /**
   * Fetches a BLOB value from the specified column.
   * @param pos A position of the column. It starts with 1.
   * @param size A pointer to an integer in which size of the returned value
   * will be stored.
   * @return a blob value. The returned value is valid
   * until calling {@link #next()} or get&lt;Type&gt;() methods for the same column
   * or destructing the instance of this class. You must copy the value
   * if you need to access it beyond such cases.
   * You must not release the returned area.
   */
  virtual void const *getTransientBlob(int pos, int *size) throw (SQLException);
};

/**
 * This class represents a prepared statement.
 * An instance of this class must exists while all instances created from
 * the instance exist.
 */
class PreparedStatement: public Statement {
  Connection *const con;
  sqlite3_stmt *const stmt;
  ResultSet *owner;
#ifdef SQLDEBUG
  std::string sql;
#endif

  PreparedStatement(Connection *con, sqlite3_stmt *stmt
#ifdef SQLDEBUG
		    , char const *sql
#endif
		    );
  friend class ResultSet;
  friend class Connection;
  virtual void reset(ResultSet *owner) throw (SQLException);
 public:
  ~PreparedStatement() throw (SQLException);
  /**
   * Returns a number of parameters in the %SQL statement supplied for
   * {@link Connection::prepare(char const *sql)}.
   * @return a number of parameters.
   */
  virtual int getParameterCount() throw (SQLException);
  /**
   * Clears all values bound to parameters.
   */
  virtual void clearParameters() throw (SQLException);
  /**
   * Binds a null value to the specified column.
   * @param pos A position of the column. It starts with 1.
   */
  virtual void setNull(int pos) throw (SQLException);
  /**
   * Binds an int value to the specified column.
   * @param pos A position of the column. It starts with 1.
   * @param value an int value.
   */
  virtual void setInt(int pos, int value) throw (SQLException);
  /**
   * Binds a UTF-8 string value to the specified column.
   * If `value' is a static data, use this method because this method is
   * more efficient than {@link #setTransientString(int pos, char const *value)}.
   * @param pos A position of the column. It starts with 1.
   * @param value a UTF-8 string value.
   */
  virtual void setStaticString(int pos, char const *value) throw (SQLException);
  /**
   * Binds a UTF-8 string value to the specified column.
   * `value' is copied in this method, thus you are free to change `value'
   * immediately after this method returns.
   * @param pos A position of the column. It starts with 1.
   * @param value a UTF-8 string value.
   */
  virtual void setTransientString(int pos, char const *value) throw (SQLException);
  /**
   * Binds a double value to the specified column.
   * @param pos A position of the column. It starts with 1.
   * @param value a double value.
   */
  virtual void setDouble(int pos, double value) throw (SQLException);
  /**
   * Binds a float value to the specified column.
   * @param pos A position of the column. It starts with 1.
   * @param value a float value.
   */
  virtual void setFloat(int pos, float value) throw (SQLException);
  /**
   * Binds a BLOB value to the specified column.
   * If `value' is a static data, use this method because this method is
   * more efficient than {@link #setTransientBlob(int pos, void const *value, int size)}.
   * @param pos A position of the column. It starts with 1.
   * @param value a BLOB value.
   * @param size a size of the `value' in bytes.
   */
  virtual void setStaticBlob(int pos, void const *value, int size) throw (SQLException);
  /**
   * Binds a BLOB value to the specified column.
   * `value' is copied in this method, thus you are free to change `value'
   * immediately after this method returns.
   * @param pos A position of the column. It starts with 1.
   * @param value a BLOB value.
   * @param size a size of the `value' in bytes.
   */
  virtual void setTransientBlob(int pos, void const *value, int size) throw (SQLException);
  /**
   * Executes a non-query statement.
   * Before calling this method, bind all parameters by calling set&lt;Type&gt;() methods.
   * @return a number of rows affected.
   */
  virtual int executeUpdate() throw (SQLException);
  /**
   * Executes a query.
   * Before calling this method, bind all parameters by calling set&lt;Type&gt;() methods.
   * @return a {@link ResultSet} instance to be used to access the result of the query.
   */
  virtual ResultSet *executeQuery() throw (SQLException);
};

/**
 * This class represents a connection to SQLite DBMS.
 * An instance of this class must exists while all instances created from
 * the instance exist.
 */
class Connection {
  sqlite3 *const db;
  Connection(sqlite3 *db);
  friend class ResultSet;
  friend class PreparedStatement;
 public:
  virtual ~Connection() throw (SQLException);
  /**
   * This is a factory method for this class.
   * @param dburl DB URL to be opened.
   * @return a connection opened.
   * @see http://www.sqlite.org/c3ref/open.html#urifilenamesinsqlite3open
   */
  static Connection *open(char const *dburl) throw (SQLException);
    
  /**
   * Executes one or more %SQL statements separated by ';'.
   */
  virtual void execute(char const *sql) throw (SQLException);
  /**
   * Creates a prepared statement.
   * @param sql an %SQL statement.
   * @return a prepared statement.
   */
  virtual PreparedStatement *prepare(char const *sql) throw (SQLException);
  /**
   * Returns the most recent successful inserted RowId.
   * This method is convenient when you want to know a INTEGER AUTOINCREMENTed
   * primary key for the last inserted row.
   * @return the most recent successful inserted RowId,
   * or 0 if no successful INSERTs have ever occurred on this connection.
   * @see http://www.sqlite.org/c3ref/last_insert_rowid.html
   * @see last_insert_rowid() described in http://www.sqlite.org/lang_corefunc.html
   */
  virtual int64_t getLastInsertRowId() throw (SQLException);
  /**
   * Creates an %SQL function.
   * @see http://www.sqlite.org/c3ref/create_function.html
   */
  virtual void createFunction(
      char const *functionName,
      int nArg,
      int eTextRep,
      void *pApp,
      void (*xFunc)(sqlite3_context *, int, sqlite3_value **),
      void (*xStep)(sqlite3_context *, int, sqlite3_value **),
      void (*xFinal)(sqlite3_context *),
      void(*xDestroy)(void *)) throw (SQLException);
};

/**
 * This class provides some useful utility functions.
 */
class SQL {
 public:
  /**
   * Returns a string which contains `values' with its order separated by `delimiter'.
   * @param values a NULL terminated array of strings.
   * @param delimiter a delimiter string.
   * @return a string which contains `values' with its order separated by `delimiter'.
  */
  static std::string join(char const * const values[], char const *delimiter=", ");

  /**
   * Returns a string which contains same number of bind characters('?') as elements in `values'
   * (excluding a trailing NULL) separated by `delimiter'.
   * @param values a NULL terminated array of pointers.
   * This is used just to figure out a number of bind characters to be listed.
   * @param delimiter a delimiter string.
   * @return a bind character list string.
   */
  template <typename T>
  static std::string bindChars(T *values[], char const *delimiter=", ") {
    return repeat<T>(values, "?", delimiter);
  }

  /**
   * Returns a string which contains same number of `rept' as elements in `values'
   * (excluding a trailing NULL) separated by `delimiter'.
   * @param values a NULL terminated array of pointers.
   * This is used just to figure out a number of `rept' to be listed.
   * @param rept a string to be listed.
   * @param delimiter a delimiter string.
   * @return a `rept' list string.
   */
  template <typename T>
  static std::string repeat(T *values[], char const *rept, char const *delimiter=", ") {
    std::string result;
    char const *sep = "";
    for (size_t i = 0; values[i]; i++) {
      result += sep;
      result += rept;
      sep = delimiter;
    }
    return result;
  }
};

}

#endif
