/*
 * logger.h
 *
 *  Created on: 2013/09/13
 *      Author: kohji
 */

#ifndef LIBSAKURA_LIBSAKURA_LOGGER_H_
#define LIBSAKURA_LIBSAKURA_LOGGER_H_

#include <log4cxx/logger.h>
#include <libsakura/sakura.h>

namespace LIBSAKURA_PREFIX {

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

struct NullLogger {
	static ::log4cxx::LoggerPtr GetLogger(char const *suffix) {
		return nullptr;
	}
	static constexpr bool IsFatalEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static constexpr bool IsErrorEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static constexpr bool IsWarnEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static constexpr bool IsInfoEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static constexpr bool IsDebugEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static constexpr bool IsTraceEnabled(::log4cxx::Logger * const logger) {
		return false;
	}
	static void Fatal(::log4cxx::LoggerPtr logger, char const *message) {
	}
	static void Error(::log4cxx::LoggerPtr logger, char const *message) {
	}
	static void Warn(::log4cxx::LoggerPtr logger, char const *message) {
	}
	static void Info(::log4cxx::LoggerPtr logger, char const *message) {
	}
	static void Debug(::log4cxx::LoggerPtr logger, char const *message) {
	}
	static void Trace(::log4cxx::LoggerPtr logger, char const *message) {
	}
};

} /* namespace LIBSAKURA_PREFIX */

#endif /* LIBSAKURA_LIBSAKURA_LOGGER_H_ */
