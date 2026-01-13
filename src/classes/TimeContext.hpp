//  TimeContext.hpp
//  Created by SHUD-up contributors.
//
#ifndef TimeContext_hpp
#define TimeContext_hpp

#include <string>

class TimeContext {
public:
    TimeContext();

    void setBaseDate(long yyyymmdd);
    long baseDate() const;

    int julianDay(double t_min) const;           /* 1-366 */
    std::string formatISO(double t_min) const;   /* YYYY-MM-DD HH:MM */

    std::string formatDate(double t_min) const;  /* YYYY-MM-DD */

private:
    long base_yyyymmdd;
    int base_year;
    unsigned base_month;
    unsigned base_day;
    long long base_days_since_epoch;

    static bool isLeapYear(int year);
    static int daysInMonth(int year, unsigned month);
    static int dayOfYear(int year, unsigned month, unsigned day);

    static long long daysFromCivil(int y, unsigned m, unsigned d);
    static void civilFromDays(long long z, int& y, unsigned& m, unsigned& d);
    static bool parseYYYYMMDD(long yyyymmdd, int& y, unsigned& m, unsigned& d);

    void toCivil(double t_min, int& y, unsigned& m, unsigned& d, int& hour, int& minute) const;
};

#endif  /* TimeContext_hpp */

