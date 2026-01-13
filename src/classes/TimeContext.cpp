#include "TimeContext.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

TimeContext::TimeContext()
    : base_yyyymmdd(0),
      base_year(0),
      base_month(0),
      base_day(0),
      base_days_since_epoch(0) {
}

void TimeContext::setBaseDate(long yyyymmdd) {
    int y = 0;
    unsigned m = 0;
    unsigned d = 0;
    if (!parseYYYYMMDD(yyyymmdd, y, m, d)) {
        base_yyyymmdd = 0;
        base_year = 0;
        base_month = 0;
        base_day = 0;
        base_days_since_epoch = 0;
        return;
    }

    base_yyyymmdd = yyyymmdd;
    base_year = y;
    base_month = m;
    base_day = d;
    base_days_since_epoch = daysFromCivil(base_year, base_month, base_day);
}

long TimeContext::baseDate() const {
    return base_yyyymmdd;
}

int TimeContext::julianDay(double t_min) const {
    int y = 0;
    unsigned m = 0;
    unsigned d = 0;
    int hour = 0;
    int minute = 0;
    toCivil(t_min, y, m, d, hour, minute);
    if (y == 0 || m == 0 || d == 0) {
        return 0;
    }
    return dayOfYear(y, m, d);
}

std::string TimeContext::formatISO(double t_min) const {
    int y = 0;
    unsigned m = 0;
    unsigned d = 0;
    int hour = 0;
    int minute = 0;
    toCivil(t_min, y, m, d, hour, minute);

    std::ostringstream oss;
    oss << std::setw(4) << std::setfill('0') << y << "-"
        << std::setw(2) << std::setfill('0') << m << "-"
        << std::setw(2) << std::setfill('0') << d << " "
        << std::setw(2) << std::setfill('0') << hour << ":"
        << std::setw(2) << std::setfill('0') << minute;
    return oss.str();
}

std::string TimeContext::formatDate(double t_min) const {
    int y = 0;
    unsigned m = 0;
    unsigned d = 0;
    int hour = 0;
    int minute = 0;
    toCivil(t_min, y, m, d, hour, minute);

    std::ostringstream oss;
    oss << std::setw(4) << std::setfill('0') << y << "-"
        << std::setw(2) << std::setfill('0') << m << "-"
        << std::setw(2) << std::setfill('0') << d;
    return oss.str();
}

bool TimeContext::isLeapYear(int year) {
    if ((year % 4) != 0) {
        return false;
    }
    if ((year % 100) != 0) {
        return true;
    }
    return (year % 400) == 0;
}

int TimeContext::daysInMonth(int year, unsigned month) {
    static const int dim[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month < 1 || month > 12) {
        return 0;
    }
    if (month == 2 && isLeapYear(year)) {
        return 29;
    }
    return dim[month - 1];
}

int TimeContext::dayOfYear(int year, unsigned month, unsigned day) {
    static const int cum[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    if (month < 1 || month > 12) {
        return 0;
    }
    int doy = cum[month - 1] + (int)day;
    if (month > 2 && isLeapYear(year)) {
        doy += 1;
    }
    return doy;
}

long long TimeContext::daysFromCivil(int y, unsigned m, unsigned d) {
    y -= m <= 2;
    const int era = (y >= 0 ? y : y - 399) / 400;
    const unsigned yoe = (unsigned)(y - era * 400);
    const unsigned doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return (long long)era * 146097 + (long long)doe - 719468;
}

void TimeContext::civilFromDays(long long z, int& y, unsigned& m, unsigned& d) {
    z += 719468;
    const long long era = (z >= 0 ? z : z - 146096) / 146097;
    const unsigned doe = (unsigned)(z - era * 146097);
    const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    y = (int)(yoe) + (int)(era * 400);
    const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const unsigned mp = (5 * doy + 2) / 153;
    d = doy - (153 * mp + 2) / 5 + 1;
    m = mp + (mp < 10 ? 3 : -9);
    y += (m <= 2);
}

bool TimeContext::parseYYYYMMDD(long yyyymmdd, int& y, unsigned& m, unsigned& d) {
    if (yyyymmdd <= 0) {
        return false;
    }
    y = (int)(yyyymmdd / 10000);
    long md = yyyymmdd % 10000;
    m = (unsigned)(md / 100);
    d = (unsigned)(md % 100);
    if (m < 1 || m > 12) {
        return false;
    }
    int dim = daysInMonth(y, m);
    if (dim <= 0 || d < 1 || (int)d > dim) {
        return false;
    }
    return true;
}

void TimeContext::toCivil(double t_min, int& y, unsigned& m, unsigned& d, int& hour, int& minute) const {
    if (base_yyyymmdd <= 0) {
        y = 0;
        m = 0;
        d = 0;
        hour = 0;
        minute = 0;
        return;
    }

    long long total_minutes = 0;
    if (!std::isnan(t_min) && !std::isinf(t_min)) {
        total_minutes = (long long)t_min;
    }

    long long day_offset = total_minutes / 1440;
    long long minute_of_day = total_minutes % 1440;
    if (minute_of_day < 0) {
        minute_of_day += 1440;
        day_offset -= 1;
    }

    long long days = base_days_since_epoch + day_offset;
    civilFromDays(days, y, m, d);

    hour = (int)(minute_of_day / 60);
    minute = (int)(minute_of_day % 60);
}

