#ifndef HOOKS_H_
#define HOOKS_H_

#include <iostream>
#include <chrono>
#include <unordered_map>
#include <string>
#include <iomanip>
#include <sstream>
#include <stack>
#include <utility>

// Global timer variable
static std::stack< std::pair<std::string, std::chrono::steady_clock::time_point> > Fields;

void HooksAddField(std::string field, std::chrono::steady_clock::time_point cur_time) {
    Fields.push(std::make_pair(field, cur_time));
}

void HooksRegionBegin(std::string name)
{
    // Add region name to the result string
    HooksAddField(name, std::chrono::steady_clock::now());
}

bool HooksRegionEnd(std::string name, size_t str_len, size_t num_len, bool newline=false)
{
    if (Fields.size() == 0 || Fields.top().first != name) {
        std::cout << "Timing function usage wrong" << std::endl;
        return false;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << std::setw(str_len) << name + " (secs)"<< " : " << std::setw(num_len) << std::fixed << std::setprecision(6) << std::chrono::duration_cast<std::chrono::microseconds>(
            end - Fields.top().second).count()/1000000.0 << std::endl;
    if (newline) std::cout << std::endl;

    Fields.pop();

    return true;
}

#endif
