#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <cassert>
#include "SQLite.h"
#ifdef SQLDEBUG
# include <iostream>
#endif

namespace sqlite {

SQLException::~SQLException() {}

Connection::Connection(sqlite3 *db): db(db) {}

Connection::~Connection() throw (SQLException) {
  int result = sqlite3_close(db); 
  if (SQLITE_OK != result) { 
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(db)); 
    throw SQLException();
  } 
}

Connection *Connection::open(char const *dburl) throw (SQLException) {
  sqlite3 *db = NULL;
  int result =
    sqlite3_open_v2(dburl, &db,
		    SQLITE_OPEN_URI|SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE|SQLITE_OPEN_PRIVATECACHE,
		    NULL); 
  if (SQLITE_OK != result) { 
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(db)); 
    throw SQLException();
  } 
  return new Connection(db);
}

void Connection::execute(char const *sql) throw (SQLException) {
#ifdef SQLDEBUG
  std::cerr << "sqlite::" << __FUNCTION__ << ": Executing " << sql << std::endl;
#endif
  char *errmsg = NULL;
  int result = sqlite3_exec(db, sql, NULL, NULL, &errmsg);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, errmsg); 
    sqlite3_free(errmsg);
    throw SQLException();
  } 
  sqlite3_free(errmsg);
}

int64_t Connection::getLastInsertRowId() throw (SQLException) {
  sqlite3_int64 result = sqlite3_last_insert_rowid(db);
  return result;
}

void Connection::createFunction(char const *functionName,
				int nArg,
				int eTextRep,
				void *pApp,
				void (*xFunc)(sqlite3_context *, int, sqlite3_value **),
				void (*xStep)(sqlite3_context *, int, sqlite3_value **),
				void (*xFinal)(sqlite3_context *),
				void(*xDestroy)(void *)) throw (SQLException) {
  int result =
    sqlite3_create_function_v2(db,
			       functionName,
			       nArg,
			       eTextRep,
			       pApp,
			       xFunc, xStep, xFinal, xDestroy);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(db)); 
    throw SQLException();
  } 
}

PreparedStatement *Connection::prepare(char const *sql) throw (SQLException) {
  sqlite3_stmt *stmt = NULL;
  const char *nextSql = NULL;
  int result = sqlite3_prepare_v2(db, sql, -1, &stmt, &nextSql);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(db)); 
    throw SQLException();
  } 
  return new PreparedStatement(this, stmt
#ifdef SQLDEBUG
			       , sql
#endif
			       );
}

Statement::~Statement() {}

PreparedStatement::PreparedStatement(Connection *con, sqlite3_stmt *stmt
#ifdef SQLDEBUG
				     , char const *sql
#endif
				     )
  : con(con), stmt(stmt), owner(NULL)
#ifdef SQLDEBUG
  , sql(sql)
#endif
{}

PreparedStatement::~PreparedStatement() throw (SQLException) {
  int result = sqlite3_finalize(stmt);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
  if (owner != NULL) {
    throw SQLException();
  }
}

int PreparedStatement::getParameterCount() throw (SQLException) {
  int result = sqlite3_bind_parameter_count(stmt);
  return result;
}

void PreparedStatement::clearParameters() throw (SQLException) {
  int result = sqlite3_clear_bindings(stmt);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setNull(int pos) throw (SQLException) {
  int result = sqlite3_bind_null(stmt, pos);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setInt(int pos, int64_t value) throw (SQLException) {
  int result = sqlite3_bind_int64(stmt, pos, value);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setFloat(int pos, float value) throw (SQLException) {
  int result = sqlite3_bind_double(stmt, pos, value);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setDouble(int pos, double value) throw (SQLException) {
  int result = sqlite3_bind_double(stmt, pos, value);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setStaticString(int pos, char const *value) throw (SQLException) {
  int result = sqlite3_bind_text(stmt, pos, value, strlen(value), SQLITE_STATIC);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setTransientString(int pos, char const *value) throw (SQLException) {
  int result = sqlite3_bind_text(stmt, pos, value, strlen(value), SQLITE_TRANSIENT);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::setStaticBlob(int pos, void const *value, int size) throw (SQLException) {
  assert(size >= 0);
  int result = sqlite3_bind_blob(stmt, pos, value, size, SQLITE_STATIC);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}
void PreparedStatement::setTransientBlob(int pos, void const *value, int size) throw (SQLException) {
  assert(size >= 0);
  int result = sqlite3_bind_blob(stmt, pos, value, size, SQLITE_TRANSIENT);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
}

void PreparedStatement::reset(ResultSet *owner) throw (SQLException) {
  if (this->owner != owner) {
    throw SQLException();
  }
  this->owner = NULL;
  int result = sqlite3_reset(stmt);
  if (SQLITE_OK != result) {
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    throw SQLException();
  } 
  //clearParameters(stmt);
}

int PreparedStatement::executeUpdate() throw (SQLException) {
#ifdef SQLDEBUG
  std::cerr << "sqlite::" << __FUNCTION__ << ": Executing " << sql << std::endl;
#endif
  int result = sqlite3_step(stmt);
  int count = -1;
  switch (result) {
  case SQLITE_DONE:
    count = sqlite3_changes(con->db);
    break;
  default:
    fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(con->db)); 
    break;
  } 
  try {
    // reset statement;
    ResultSet rs(this);
    owner = &rs;
  } catch (...) {
    throw;
  }
  switch (result) {
  case SQLITE_DONE:
    return count;
  } 
  throw SQLException();
}

ResultSet *PreparedStatement::executeQuery() throw (SQLException) {
#ifdef SQLDEBUG
  std::cerr << "sqlite::" << __FUNCTION__ << ": Executing " << sql << std::endl;
#endif
  if (owner != NULL) {
    throw SQLException();
  }

  ResultSet *rs = new ResultSet(this);
  owner = rs;
  return rs;
}

ResultSet::ResultSet(PreparedStatement *stmt) : stmt(stmt) {}

ResultSet::~ResultSet() throw (SQLException) {
  stmt->reset(this);
}

bool ResultSet::next() throw (SQLException) {
  int result = sqlite3_step(stmt->stmt);
  switch (result) {
  case SQLITE_DONE:
    return false;
  case SQLITE_ROW:
    return true;
  } 
  fprintf(stderr, "%s: %d %s\n", __FUNCTION__, result, sqlite3_errmsg(stmt->con->db)); 
  throw SQLException();
}

int ResultSet::getColumnCount() throw (SQLException) {
  int columns = sqlite3_column_count(stmt->stmt);
  return columns;
}

std::string ResultSet::getColumnName(int pos) throw (SQLException) {
  char const *name = sqlite3_column_name(stmt->stmt, pos - 1);
  if (name == NULL) {
    fprintf(stderr, "%s: %s\n", __FUNCTION__, sqlite3_errmsg(stmt->con->db)); 
    throw SQLException();
  }
  std::string result(name);
  return result;
}

bool ResultSet::isNull(int pos) throw (SQLException) {
  int type = sqlite3_column_type(stmt->stmt, pos - 1);
  if (SQLITE_NULL == type) {
    return true;
  }
  return false;
}

int64_t ResultSet::getInt(int pos) throw (SQLException) {
  sqlite3_int64 value = sqlite3_column_int64(stmt->stmt, pos - 1);
  return value;
}

float ResultSet::getFloat(int pos) throw (SQLException) {
  double value = sqlite3_column_double(stmt->stmt, pos - 1);
  return value;
}

double ResultSet::getDouble(int pos) throw (SQLException) {
  double value = sqlite3_column_double(stmt->stmt, pos - 1);
  return value;
}

char const *ResultSet::getTransientString(int pos, int *size) throw (SQLException) {
  unsigned char const *value = sqlite3_column_text(stmt->stmt, pos - 1);
  *size = sqlite3_column_bytes(stmt->stmt, pos - 1);
  return (char const *)value;
}

void const *ResultSet::getTransientBlob(int pos, int *size) throw (SQLException) {
  void const *value = sqlite3_column_blob(stmt->stmt, pos - 1);
  *size = sqlite3_column_bytes(stmt->stmt, pos - 1);
  return value;
}

std::string SQL::join(char const * const values[], char const *delimiter) {
  std::string result;
  char const *sep = "";
  for (size_t i = 0; values[i]; i++) {
    result += sep;
    result += values[i];
    sep = delimiter;
  }
  return result;
}

}
