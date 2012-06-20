#ifndef SQLITE_H_
#define SQLITE_H_

#include <string>
#include <sqlite3.h>

#ifdef DEBUG
#define SQLDEBUG 1
#endif

namespace sqlite {

// All strings should be UTF-8
/**
 * Lifetime of following objects should be:
 * ResultSet within PreparedStatement within Connection.
 */

class SQLException {
};

class Statement {
 public:
  virtual ~Statement();
};

class PreparedStatement;
class Connection;

class ResultSet {
  PreparedStatement *const stmt;

  ResultSet(PreparedStatement *stmt);
  friend class PreparedStatement;
 public:
  virtual ~ResultSet() throw (SQLException);
  virtual bool next() throw (SQLException);
  virtual int getColumnCount() throw (SQLException);
  // pos starts with 1.
  virtual std::string getColumnName(int pos) throw (SQLException);
  virtual bool wasNull(int pos) throw (SQLException);
  virtual int getInt(int pos) throw (SQLException);
  // Don't release returned area
  virtual char const *getTransientString(int pos, int *size) throw (SQLException);
  virtual float getFloat(int pos) throw (SQLException);
  virtual double getDouble(int pos) throw (SQLException);
  // Don't release returned area
  virtual void const *getTransientBlob(int pos, int *size) throw (SQLException);
};

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
  virtual int getParameterCount() throw (SQLException);
  virtual void clearParameters() throw (SQLException);
  // pos starts with 1.
  virtual void setNull(int pos) throw (SQLException);
  virtual void setInt(int pos, int value) throw (SQLException);
  virtual void setStaticString(int pos, char const *value) throw (SQLException);
  virtual void setTransientString(int pos, char const *value) throw (SQLException);
  virtual void setDouble(int pos, double value) throw (SQLException);
  virtual void setFloat(int pos, float value) throw (SQLException);
  virtual void setStaticBlob(int pos, void const *value, int size) throw (SQLException);
  virtual void setTransientBlob(int pos, void const *value, int size) throw (SQLException);
  virtual int executeUpdate() throw (SQLException);
  virtual ResultSet *executeQuery() throw (SQLException);
  ~PreparedStatement() throw (SQLException);
};

class Connection {
  sqlite3 *const db;
  Connection(sqlite3 *db);
  friend class ResultSet;
  friend class PreparedStatement;
 public:
  virtual ~Connection() throw (SQLException);
  static Connection *open(char const *dburl) throw (SQLException);
    
  virtual void execute(char const *sql) throw (SQLException);
  virtual PreparedStatement *prepare(char const *sql) throw (SQLException);
  // See the document of sqlite3_create_function_v2().
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

class SQL {
 public:
  static std::string join(char const * const values[], char const *delimiter=", ");

  template <typename T>
  static std::string bindChars(T *values[], char const *delimiter=", ") {
    return repeat<T>(values, "?", delimiter);
  }

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
