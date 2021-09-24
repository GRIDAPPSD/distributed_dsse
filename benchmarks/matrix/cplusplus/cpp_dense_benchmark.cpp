
#include <string>
using std::string;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
#include <sstream>

#include <unistd.h>
#include <sys/time.h>

#include <unordered_map>
using std::unordered_map;

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using Eigen::MatrixXd;

void process_mem_usage(string& vm_used) {
    using std::ios_base;
    using std::ifstream;
    using std::string;

    // 'file' stat seems to give the most reliable results
    ifstream stat_stream("/proc/self/stat",ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    unsigned long vsize;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize; // don't care about the rest
    stat_stream.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    uint vm_usage = (int)(vsize / (1024.0*1024.0));
    string units = " MB";
    vm_used = std::to_string(vm_usage) + units;
}


double getWallTime() {
    struct timeval time;

    if (gettimeofday(&time, NULL))
        return 0;

    return (double)time.tv_sec + (double)time.tv_usec*0.000001;
}


string getMinSec(double seconds) {
    uint useconds = (uint)seconds;
    uint sec = useconds % 60;
    string secstr = (sec<10)? "0"+std::to_string(sec): std::to_string(sec);
    return (std::to_string(useconds/60) + ":" + secstr + std::to_string(seconds - useconds).substr(1));
}


int main(int argc, char** argv) {

    double startTime, elapsedTime;
    string vm_used;

    string line, istr, jstr, valstr;
    uint i, j;
    double val;

    std::ifstream DimFile("../mat/dimensions.csv");
    std::getline(DimFile, line);
    std::istringstream stream(line);
    std::getline(stream, istr, ','); const uint xqty = std::stoi(istr);
    std::getline(stream, jstr, ','); const uint zqty = std::stoi(jstr);

    {
        MatrixXd Supd(zqty,zqty);
        std::ifstream SupdFile("../mat/Supd_trip.csv");
        while (std::getline(SupdFile, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            Supd.coeffRef(i,j) = val;
        }

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();
            MatrixXd SupdInv = Supd.inverse();
            elapsedTime = getWallTime() - startTime;
            cout << "C++ Dense Supd inverse time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Dense Supd inverse memory: " << vm_used << "\n" << std::flush;
    }

    //process_mem_usage(vm_used);
    //cout << "C++ Dense After delete memory: " << vm_used << "\n" << std::flush;

    {
        MatrixXd K2(xqty,zqty);
        std::ifstream K2File("../mat/K2_trip.csv");
        while (std::getline(K2File, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            K2.coeffRef(i,j) = val;
        }

        MatrixXd K3(zqty,zqty);
        std::ifstream K3File("../mat/K3_trip.csv");
        while (std::getline(K3File, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            K3.coeffRef(i,j) = val;
        }

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();
            MatrixXd Kupd = K2 * K3;
            elapsedTime = getWallTime() - startTime;
            cout << "C++ Dense Kupd multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Dense Kupd multiply memory: " << vm_used << "\n" << std::flush;
    }

    //process_mem_usage(vm_used);
    //cout << "C++ Dense After delete memory: " << vm_used << "\n" << std::flush;

    {
        MatrixXd K2(xqty,zqty);
        std::ifstream K2File("../mat/K2_trip.csv");
        while (std::getline(K2File, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            K2.coeffRef(i,j) = val;
        }

        MatrixXd K3sparse(zqty,zqty);
        std::ifstream K3File("../mat/K3sparse_trip.csv");
        while (std::getline(K3File, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            K3sparse.coeffRef(i,j) = val;
        }

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();
            MatrixXd Kupd_sparse = K2 * K3sparse;
            elapsedTime = getWallTime() - startTime;
            cout << "C++ Dense Kupd_sparse multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Dense Kupd_sparse multiply memory: " << vm_used << "\n" << std::flush;
    }

    //process_mem_usage(vm_used);
    //cout << "C++ Dense After delete memory: " << vm_used << "\n" << std::flush;

}

