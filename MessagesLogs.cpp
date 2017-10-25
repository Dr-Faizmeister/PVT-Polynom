#include "stdafx.h"
#include "MessagesLogs.h"

MessageLogger* MessageLogger::inst = 0;

/*string MessageType2String(MessageType MsgType) {
	switch (MsgType) {
	case LOG_ERROR:
		return "ERROR";
	case LOG_WARNING:
		return "WARNING";
	case LOG_INFO:
		return "MESSAGE";
	case LOG_DEBUG_INFO:
		return "DEBUG MESSAGE";
	case LOG_TIME_MESSAGE:
		return "TIME_MESSAGE";
	default:
		return "Unknow message type";
	}
}*/

string MessageType2String(MessageType MsgType) {
	switch (MsgType) {
	case LOG_ERROR:
		return ERROR_STRING;
	case LOG_WARNING:
		return WARNING_STRING;
	case LOG_INFO:
		return INFO_STRING;
	case LOG_DEBUG_INFO:
		return DEBUG_INFO_STRING;
	case LOG_TIME_MESSAGE:
		return TIME_MESSAGE_STRING;
	default:
		return "Unknow message type";
	}
}

void MessageLogger::setOutputFile(string Filename) {
	fOut = std::ofstream(Filename.c_str());
}


bool MessageLogger::isPrintToConsole(MessageType MsgType) {
#ifdef _DEBUG
	return true;
#endif
	return MsgType == LOG_ERROR || MsgType == LOG_INFO;
}

MessageLogger& MessageLogger::operator<<(MessageType MsgType) {
	currentMsgType = MsgType;
	fOut << MessageType2String(MsgType);
	fOut.flush();
	return *this;
}

MessageLogger& MessageLogger::operator<<(char* Str) {
	fOut << Str;
	if (isPrintToConsole(currentMsgType)) {
		std::cout << Str;
	}
	fOut.flush();
	return *this;
}

MessageLogger& MessageLogger::operator<<(string Str) {
	fOut << Str;
	if (isPrintToConsole(currentMsgType)) {
		std::cout << Str;
	}
	fOut.flush();
	return *this;
}

MessageLogger& MessageLogger::operator<<(double Value) {
	fOut << Value;
	if (isPrintToConsole(currentMsgType)) {
		std::cout << Value;
	}
	fOut.flush();
	return *this;
}

MessageLogger& MessageLogger::operator<<(int Value) {
	fOut << Value;
	if (isPrintToConsole(currentMsgType)) {
		std::cout << Value;
	}
	fOut.flush();
	return *this;
}



std::string GetLastErrorAsString() {
    //Get the error message, if any.
    DWORD errorMessageID = GetLastError();
    if(errorMessageID == 0)
        return std::string(); //No error message has been recorded

    LPSTR messageBuffer = nullptr;
    size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                                 NULL, errorMessageID, MAKELANGID(LANG_ENGLISH, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);

    std::string message(messageBuffer, size);

    //Free the buffer.
    LocalFree(messageBuffer);

    return message;
}