/*
 * logger.h
 *
 *  Created on: 2013/09/13
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_LOGGER_H_
#define LIBSAKURA_LIBSAKURA_LOGGER_H_

#include <libsakura/config.h>

#if LIBSAKURA_HAS_LOG4CXX
# include <log4cxx/logger.h>
#endif

namespace LIBSAKURA_PREFIX {

#if LIBSAKURA_HAS_LOG4CXX
typedef ::log4cxx::Logger LoggerType;
typedef ::log4cxx::LoggerPtr LoggerPtrType;
struct Logger {
	static ::log4cxx::LoggerPtr GetLogger(char const *suffix = "");
	static bool IsFatalEnabled(::log4cxx::Logger * const logger) {
		return logger->isFatalEnabled();
	}
	static bool IsErrorEnabled(::log4cxx::Logger * const logger) {
		return logger->isErrorEnabled();
	}
	static bool IsWarnEnabled(::log4cxx::Logger * const logger) {
		return logger->isWarnEnabled();
	}
	static bool IsInfoEnabled(::log4cxx::Logger * const logger) {
		return logger->isInfoEnabled();
	}
	static bool IsDebugEnabled(::log4cxx::Logger * const logger) {
		return logger->isDebugEnabled();
	}
	static bool IsTraceEnabled(::log4cxx::Logger * const logger) {
		return logger->isTraceEnabled();
	}
	static void Fatal(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_FATAL(logger, message);
	}
	static void Error(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_ERROR(logger, message);
	}
	static void Warn(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_WARN(logger, message);
	}
	static void Info(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_INFO(logger, message);
	}
	static void Debug(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_DEBUG(logger, message);
	}
	static void Trace(::log4cxx::LoggerPtr logger, char const *message) {
		LOG4CXX_TRACE(logger, message);
	}
};

#else	/* LIBSAKURA_HAS_LOG4CXX */
class LoggerType;
typedef void *LoggerPtrType;
struct NullLogger;
typedef NullLogger Logger;

# define LOG4CXX_ERROR(logger, msg) ((void)0)
# define LOG4CXX_WARN(logger, msg) ((void)0)
# define LOG4CXX_INFO(logger, msg) ((void)0)
# define LOG4CXX_DEBUG(logger, msg) ((void)0)

#endif	/* LIBSAKURA_HAS_LOG4CXX */

struct NullLogger {
	static LoggerPtrType GetLogger(char const *suffix) {
		return nullptr;
	}
	static constexpr bool IsFatalEnabled(LoggerType * const logger) {
		return false;
	}
	static constexpr bool IsErrorEnabled(LoggerType * const logger) {
		return false;
	}
	static constexpr bool IsWarnEnabled(LoggerType * const logger) {
		return false;
	}
	static constexpr bool IsInfoEnabled(LoggerType * const logger) {
		return false;
	}
	static constexpr bool IsDebugEnabled(LoggerType * const logger) {
		return false;
	}
	static constexpr bool IsTraceEnabled(LoggerType * const logger) {
		return false;
	}
	static void Fatal(LoggerPtrType logger, char const *message) {
	}
	static void Error(LoggerPtrType logger, char const *message) {
	}
	static void Warn(LoggerPtrType logger, char const *message) {
	}
	static void Info(LoggerPtrType logger, char const *message) {
	}
	static void Debug(LoggerPtrType logger, char const *message) {
	}
	static void Trace(LoggerPtrType logger, char const *message) {
	}
};

} /* namespace LIBSAKURA_PREFIX */

#endif /* LIBSAKURA_LIBSAKURA_LOGGER_H_ */
