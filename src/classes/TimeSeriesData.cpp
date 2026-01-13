#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits>
#include "TimeSeriesData.hpp"

void CheckFile(std::ifstream * fp, const char *s)
{
    if (fp == NULL || !fp->is_open()) {
        fprintf(stderr, "\n  Fatal Error: \n %s is in use or does not exist!\n", s);
        myexit(ERRFileIO);
    }
}

_TimeSeriesData::_TimeSeriesData(){
    timeRangeCached = 0;
    minTime = NA_VALUE;
    maxTime = NA_VALUE;
}

_TimeSeriesData::~_TimeSeriesData(){
    for (int i = 0; i < MAXQUE; i++) {
        delete[] ts[i];
    }
#ifdef DEBUG
    printf("TSD to %s Destructed.\n", fn.c_str());
#endif
}

void _TimeSeriesData::initialize(int n)
{
    ncol = n; /* Include the Time and Data columns */
    
    iNow = 0;
    nQue = 0;
    eof = 0;
    timeRangeCached = 0;
    minTime = NA_VALUE;
    maxTime = NA_VALUE;
    for (int i = 0; i < MAXQUE + 1; i++) {
        pRing[i] = i + 1;
    }
    pRing[MAXQUE] = 0;
    //for repeating last value;
    iNow = 0;
    iNext = pRing[iNow];
    
    for (int i = 0; i < MAXQUE + 1; i++) {
        ts[i] = new double[n];
    }
}

void _TimeSeriesData::computeTimeRange() const
{
    if (timeRangeCached) {
        return;
    }

    std::ifstream file(fn);
    CheckFile(&file, fn.c_str());

    std::string line;
    if (!std::getline(file, line)) {
        fprintf(stderr, "\n  Fatal Error: Empty time-series file.\n");
        fprintf(stderr, "  File: %s\n", fn.c_str());
        myexit(ERRFileIO);
    }
    if (!std::getline(file, line)) {
        fprintf(stderr, "\n  Fatal Error: Missing header line in time-series file.\n");
        fprintf(stderr, "  File: %s\n", fn.c_str());
        myexit(ERRFileIO);
    }

    bool hasData = false;
    double prevTimeMin = 0.0;
    long prevLineNo = 0;
    long lineNo = 2; // already consumed 2 lines

    double localMin = std::numeric_limits<double>::infinity();
    double localMax = -std::numeric_limits<double>::infinity();

    while (std::getline(file, line)) {
        lineNo++;
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#') {
            continue;
        }

        std::istringstream iss(line);
        double timeDay = 0.0;
        if (!(iss >> timeDay)) {
            fprintf(stderr, "\n  Fatal Error: Failed to parse time value.\n");
            fprintf(stderr, "  File: %s\n", fn.c_str());
            fprintf(stderr, "  Line: %ld\n", lineNo);
            fprintf(stderr, "  Content: %s\n", line.c_str());
            myexit(ERRDATAIN);
        }
        const double timeMin = timeDay * 1440.0; // Day -> minute

        if (!hasData) {
            hasData = true;
            prevTimeMin = timeMin;
            prevLineNo = lineNo;
            localMin = timeMin;
            localMax = timeMin;
            continue;
        }

        if (timeMin + 1e-12 < prevTimeMin) {
            fprintf(stderr, "\n  Fatal Error: Time column is not monotonic non-decreasing.\n");
            fprintf(stderr, "  File: %s\n", fn.c_str());
            fprintf(stderr, "  Previous line: %ld, time = %.15g day (%.3f min)\n",
                    prevLineNo, prevTimeMin / 1440.0, prevTimeMin);
            fprintf(stderr, "  Current  line: %ld, time = %.15g day (%.3f min)\n",
                    lineNo, timeMin / 1440.0, timeMin);
            fprintf(stderr, "  Fix: sort the time column ascending (allow equal), or regenerate the file.\n");
            myexit(ERRDATAIN);
        }

        if (timeMin < localMin) localMin = timeMin;
        if (timeMin > localMax) localMax = timeMin;
        prevTimeMin = timeMin;
        prevLineNo = lineNo;
    }

    if (!hasData) {
        fprintf(stderr, "\n  Fatal Error: No data rows in time-series file.\n");
        fprintf(stderr, "  File: %s\n", fn.c_str());
        fprintf(stderr, "  Fix: ensure there are numeric data rows after the 2-line header.\n");
        myexit(ERRFileIO);
    }

    minTime = localMin;
    maxTime = localMax;
    timeRangeCached = 1;
}

double _TimeSeriesData::getMinTime() const
{
    computeTimeRange();
    return minTime;
}

double _TimeSeriesData::getMaxTime() const
{
    computeTimeRange();
    return maxTime;
}

void _TimeSeriesData::read_csv()
{
    if (!eof) {
        computeTimeRange(); // validates monotonic time column and caches min/max time
        std::ifstream file(fn);
        CheckFile(&file, fn.c_str());
        std::string str;
        Length = 0;
#ifdef DEBUG
        if (nQue > 0) {
            std::cout << "No of Queue = " << nQue << std::endl;;
        }
#endif
        for (int i = 0; i < MAXQUE * nQue + 2; i++) {
            /* Line 1= size of table; Line 2= Head of table */
            getline(file, str);
        }
        long lineNo = (long)(MAXQUE * nQue + 2); // 1-based line number of last skipped line
        bool hasPrev = false;
        double prevTimeMin = 0.0;
        long prevLineNo = 0;
        for (int i = 0; i < MAXQUE && getline(file, str); ) {
            lineNo++;

            const size_t first = str.find_first_not_of(" \t\r\n");
            if (first == std::string::npos || str[first] == '#') {
                continue;
            }

            std::istringstream iss(str);
            double timeDay = 0.0;
            if (!(iss >> timeDay)) {
                fprintf(stderr, "\n  Fatal Error: Failed to parse time value.\n");
                fprintf(stderr, "  File: %s\n", fn.c_str());
                fprintf(stderr, "  Line: %ld\n", lineNo);
                fprintf(stderr, "  Content: %s\n", str.c_str());
                myexit(ERRDATAIN);
            }
            const double timeMin = timeDay * 1440.0; // Day -> minute

            if (hasPrev && timeMin + 1e-12 < prevTimeMin) {
                fprintf(stderr, "\n  Fatal Error: Time column is not monotonic non-decreasing.\n");
                fprintf(stderr, "  File: %s\n", fn.c_str());
                fprintf(stderr, "  Previous line: %ld, time = %.15g day (%.3f min)\n",
                        prevLineNo, prevTimeMin / 1440.0, prevTimeMin);
                fprintf(stderr, "  Current  line: %ld, time = %.15g day (%.3f min)\n",
                        lineNo, timeMin / 1440.0, timeMin);
                fprintf(stderr, "  Fix: sort the time column ascending (allow equal), or regenerate the file.\n");
                myexit(ERRDATAIN);
            }

            ts[i][0] = timeMin;
            for (int j = 1; j < ncol; j++) {
                if (!(iss >> ts[i][j])) {
                    fprintf(stderr, "\n  Fatal Error: Failed to parse numeric column %d.\n", j + 1);
                    fprintf(stderr, "  File: %s\n", fn.c_str());
                    fprintf(stderr, "  Line: %ld\n", lineNo);
                    fprintf(stderr, "  Content: %s\n", str.c_str());
                    myexit(ERRDATAIN);
                }
            }

            hasPrev = true;
            prevTimeMin = timeMin;
            prevLineNo = lineNo;
            Length++;
            i++;
        }
        nQue++;            /* Number of the Queue was reload */
        if (!file.eof()) {
            eof = 0;
        } else {
            eof = 1;
        }
#ifdef DEBUG
        std::cout << fn << "\n\tUpdate queue. Length = " << Length << std::endl;
#endif
        if(Length <= 0){
            fprintf(stderr, "Reading fail, file = %s\n", fn.c_str());
            myexit(ERRFileIO);
        }
        file.close();
    }
}
void _TimeSeriesData::readDimensions()
{
    int tmp, nc;
    FILE *fp = fopen(fn.c_str(), "r");
    CheckFile(fp, fn.c_str());
    fscanf(fp, "%d %d %ld", &tmp, &nc, &StartTime);
#ifdef DEBUG
    fprintf(stdout, "Header of %s : %d\t%d\t%ld ",fn.c_str(),  tmp, nc, StartTime);
#endif
    fclose(fp);
    initialize(nc);
}
double _TimeSeriesData::getX(double t, int col)
{
    return ts[iNow][col];
}
int _TimeSeriesData::get_Ncol(){
    return ncol;
}
long _TimeSeriesData::getStartTime() const
{
    return StartTime;
}
void _TimeSeriesData::applyCalib(double prcp, double temp)
{
    for (int i = 0; i < Length; i++) {
        ts[i][1] *= prcp;    /* Calibration of prcp */
        ts[i][2] += temp;    /* Calibration of temp */
    }
}
void _TimeSeriesData::movePointer(double t){
    while (t >= ts[iNext][0] && ts[iNext][0] >= ts[iNow][0]) {
        if (iNow == 0) {
            for (int i = 0; i < ncol; i++) {
                ts[Length][i] = ts[Length - 1][i];
                //repeat last item;
            }
        }
        iNow = iNext;
        iNext = pRing[iNow];
        if (iNext == 0) {
            read_csv();
        }
        //printf("%s,  %.1f \t [%.4f] \t %.1f\n", fn.c_str(), ts[iNow][0], t, ts[iNext][0]);
    }
    if (ts[iNext][0] < ts[iNow][0] && t - ts[iNow][0] > 1 && ts[iNow][0] + 1440 < t) {
        fprintf(stderr, "\nError in reading file: %s\n", fn.c_str());
        fprintf(stderr, "\nError: missing forcing data after t=%.3lf\n\n", ts[iNow][0] / 1440. + 1);
        myexit(ERRFileIO);
    }
}

void _TimeSeriesData::checkValue(int icol, double xmin, double xmax, const char *varname){
    for(int i = 0; i < MAXQUE & i < Length; i++){
        if( ts[i][icol] < xmin || ts[i][icol] > xmax){
            fprintf(stderr, "Warning: value of %s(t=%g min) = %g is out of range (%f, %f).\n", varname, ts[i][0], ts[i][icol], xmin, xmax);
        }
    }
}
