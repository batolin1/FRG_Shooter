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

class Logger {
public:
    static Logger& instance() {
        static Logger logger;
        return logger;
    }

    void log(const std::string& message) {
        std::lock_guard<std::mutex> lock(mtx);
        std::string timestamped = "[" + current_time() + "] " + message;
        std::cout << timestamped << "\n";
        messages.push_back(timestamped);
    }

    std::vector<std::string> get_messages() const {
        std::lock_guard<std::mutex> lock(mtx);
        return messages;
    }

    void clear() {
        std::lock_guard<std::mutex> lock(mtx);
        messages.clear();
    }

private:
    std::vector<std::string> messages;
    mutable std::mutex mtx;

    // Not thread-safe on its own (uses a thread-safe localtime variant),
    // but is only ever called while `mtx` is held by log(), so concurrent
    // calls from multiple threads never overlap.
    std::string current_time() const {
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        std::tm tm_buf{};
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