#ifndef FILE_PARSING_EXCEPTION_HPP
#define FILE_PARSING_EXCEPTION_HPP

#include <exception>
#include <string>

/**
    @brief a file parsing exception to be thrown if parsing of file fails. 
*/
class File_Parsing_Exception : public std::exception {
        
    private:

        std::string message;

    public:
    
        explicit File_Parsing_Exception (const std::string& msg) : message (msg) {}

        const char* what () const noexcept override {
            return message.c_str ();
        }
};

#endif