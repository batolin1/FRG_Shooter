#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>
#include <vector>
#include <mutex>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>

class Logger {
public:
    static Logger& instance() {
        static Logger logger;
        return logger;
    }

    void log(const std::string& message) {
        std::string timestamped = "[" + current_time() + "] " + message;
        std::cout << timestamped << std::endl;
        std::lock_guard<std::mutex> lock(mtx);
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

    std::string current_time() const {
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        std::ostringstream oss;
        oss << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }
};
#endif