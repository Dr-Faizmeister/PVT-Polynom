#pragma once
#include "stdafx.h"
#include "common.h"

#define ERROR_STRING "ERROR : "
#define WARNING_STRING "WARNING : "
#define INFO_STRING "MESSAGE : "
#define DEBUG_INFO_STRING "DEBUG_MESSAGE : "
#define TIME_MESSAGE_STRING "TIME_MESSAGE : "

using std::string;

// Типы сообщений
enum MessageType {
	LOG_ERROR,			// Ошибка. Отображать в консоль и лог.
	LOG_WARNING,		// Предупреждение. В зависимости от настроек выводится в консоль.
	LOG_INFO,			// Сообщение. Выводится в консоль.
	LOG_DEBUG_INFO,		// Отладочный вывод. Работает при включенной опции отладки.
	LOG_TIME_MESSAGE	// Лог в ходе расчета, отображает процесс решения
};

class MessageLogger {
public:
	static MessageLogger* Inst() {
		if (inst == 0) {
			inst = new MessageLogger();
		}
		return inst;
	}

	void setOutputFile(string Filename);

	MessageLogger& operator<<(MessageType MsgType);
	MessageLogger& operator<<(char* Str);
	MessageLogger& operator<<(string Str);
	MessageLogger& operator<<(double Value);
	MessageLogger& operator<<(int Value);
private:
	MessageLogger() {};
	static MessageLogger* inst;

	std::ofstream fOut;
	MessageType currentMsgType;

	bool isPrintToConsole(MessageType MsgType);
};

#define WMSG (*MessageLogger::Inst())

//Returns the last Win32 error, in string format. Returns an empty string if there is no error.
std::string GetLastErrorAsString();