#ifndef LOGGER_HPP
#define LOGGER_HPP
#include <string>
#include <vector>
#include <mutex>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <iostream>

/**
    @brief Simple logger class to log the actions realized throughout an execution onto a file. 
*/
class Logger {
public:

    /**
        @brief Calls an instance of this logger. 
        @return    An instance of this object.
    */
    static Logger& instance() {
        static Logger logger;
        return logger;
    }

    /**
        @brief Logs a string message to the logger. 
        @param message    The message to be logged. 
    */
    void log(const std::string& message) {
        std::lock_guard<std::mutex> lock (mtx);
        std::string timestamped = "[" + current_time() + "] " + message;
        std::cout << timestamped << "\n";
        messages.push_back (timestamped);
    }

    /**
        @brief Getter method to get all the messages stored in logger, as a vector of strings.
        @return    Vector of strings containing messages. 
    */
    std::vector<std::string> get_messages () const {
        std::lock_guard<std::mutex> lock(mtx);
        return messages;
    }

    /**
        @brief Clears this logger. 
    */
    void clear() {
        std::lock_guard<std::mutex> lock (mtx);
        messages.clear ();
    }

private:
    std::vector<std::string> messages;
    mutable std::mutex mtx;

    /**
        @brief Gets the current time. 
        @return    A string with the current time. 
    */
    std::string current_time () const {
        auto now = std::chrono::system_clock::now ();
        std::time_t t = std::chrono::system_clock::to_time_t (now);
        std::tm tm_buf {};
#if defined(_WIN32)
        localtime_s(&tm_buf, &t);
#else
        localtime_r(&t, &tm_buf);
#endif
        std::ostringstream oss;
        oss << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }
};
#endif