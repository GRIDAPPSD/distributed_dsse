
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
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::SparseLU;


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
    std::istringstream stream1(line);
    std::getline(stream1, istr, ','); const uint xqty = std::stoi(istr);
    std::getline(stream1, jstr, ','); const uint zqty = std::stoi(jstr);

    std::getline(DimFile, line);
    std::istringstream stream2(line);
    std::getline(stream2, istr, ','); const uint Supd_entries = std::stoi(istr);

    std::getline(DimFile, line);
    std::istringstream stream3(line);
    std::getline(stream3, istr, ','); const uint K2_entries = std::stoi(istr);

    std::getline(DimFile, line);
    std::istringstream stream4(line);
    std::getline(stream4, istr, ','); const uint K3_entries = std::stoi(istr);

    std::getline(DimFile, line);
    std::istringstream stream5(line);
    std::getline(stream5, istr, ','); const uint K3sparse_entries = std::stoi(istr);

    {
        SparseMatrix<double> Supd(zqty,zqty);
        Supd.reserve(Supd_entries);

        std::ifstream SupdFile("../mat/Supd_trip.csv");
        while (std::getline(SupdFile, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            Supd.coeffRef(i,j) = val;
        }
        Supd.makeCompressed();

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();

            SparseLU<SparseMatrix<double> > solver;
            solver.compute(Supd);
            SparseMatrix<double> I(zqty,zqty);
            I.setIdentity();
            auto SupdInv = solver.solve(I);

            elapsedTime = getWallTime() - startTime;
            cout << "C++ Sparse Supd SparseI inverse time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Sparse Supd SparseI inverse memory: " << vm_used << "\n" << std::flush;

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();

            SparseLU<SparseMatrix<double> > solver;
            solver.compute(Supd);
            MatrixXd I(zqty,zqty);
            I.setIdentity();
            auto SupdInv = solver.solve(I);

            elapsedTime = getWallTime() - startTime;
            cout << "C++ Sparse Supd DenseI inverse time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Sparse Supd DenseI inverse memory: " << vm_used << "\n" << std::flush;
    }

    exit(0);
    //process_mem_usage(vm_used);
    //cout << "C++ Sparse After delete memory: " << vm_used << "\n" << std::flush;

    {
        cout << "C++ Sparse reading " << K2_entries << " K2 entries..." << "\n" << std::flush;
        SparseMatrix<double> K2(xqty,zqty);
        K2.reserve(K2_entries);
        std::ifstream K2file("../mat/K2_trip.csv");
        uint entries = 0;
        while (std::getline(K2file, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            //cout << "K2 i: " << i << ", j: " << j << ", val: " << val << endl;
            K2.coeffRef(i,j) = val;
            if (++entries % 10000 == 0)
                cout << "C++ Sparse read " << entries << " K2 entries (" << (int)(100.0*((double)entries/(double)K2_entries)) << "%)\n" << std::flush;
        }
        K2.makeCompressed();

#if 111
        MatrixXd K3dense(zqty,zqty);
        std::ifstream K3File("../mat/K3_trip.csv");
        while (std::getline(K3File, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            K3dense.coeffRef(i,j) = val;
        }
        SparseMatrix<double> K3 = K3dense.sparseView();
        K3.makeCompressed();
#else
        cout << "C++ Sparse reading " << K3_entries << " K3 entries..." << "\n" << std::flush;
        SparseMatrix<double> K3(zqty,zqty);
        K3.reserve(K3_entries);
        std::ifstream K3file("../mat/K3_trip.csv");
        entries = 0;
        while (std::getline(K3file, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            //cout << "K3 i: " << i << ", j: " << j << ", val: " << val << endl;
            K3.coeffRef(i,j) = val;
            if (++entries % 10000 == 0)
                cout << "C++ Sparse read " << entries << " K3 entries (" << (int)(100.0*((double)entries/(double)K3_entries)) << "%)\n" << std::flush;
        }
        K3.makeCompressed();
#endif
        cout << "C++ Sparse Kupd about to start multiply..." << "\n" << std::flush;

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();

            SparseMatrix<double> Kupd = K2 * K3;

            elapsedTime = getWallTime() - startTime;
            cout << "C++ Sparse Kupd multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;
       }

        process_mem_usage(vm_used);
        cout << "C++ Sparse Kupd multiply memory: " << vm_used << "\n" << std::flush;
    }

    //process_mem_usage(vm_used);
    //cout << "C++ Sparse After delete memory: " << vm_used << "\n" << std::flush;

    {
        cout << "C++ Sparse reading " << K2_entries << " K2 entries..." << "\n" << std::flush;
        SparseMatrix<double> K2(xqty,zqty);
        K2.reserve(K2_entries);
        std::ifstream K2file("../mat/K2_trip.csv");
        uint entries = 0;
        while (std::getline(K2file, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            //cout << "K2 i: " << i << ", j: " << j << ", val: " << val << endl;
            K2.coeffRef(i,j) = val;
            if (++entries % 10000 == 0)
                cout << "C++ Sparse read " << entries << " K2 entries (" << (int)(100.0*((double)entries/(double)K2_entries)) << "%)\n" << std::flush;
        }
        K2.makeCompressed();

        cout << "C++ Sparse reading " << K3sparse_entries << " K3sparse entries..." << "\n" << std::flush;
        SparseMatrix<double> K3sparse(zqty,zqty);
        K3sparse.reserve(K3sparse_entries);
        std::ifstream K3file("../mat/K3sparse_trip.csv");
        entries = 0;
        while (std::getline(K3file, line)) {
            std::istringstream stream(line);
            std::getline(stream, istr, ','); i = std::stoi(istr);
            std::getline(stream, jstr, ','); j = std::stoi(jstr);
            std::getline(stream, valstr, ','); val = std::stod(valstr);
            //cout << "K3sparse i: " << i << ", j: " << j << ", val: " << val << endl;
            K3sparse.coeffRef(i,j) = val;
            if (++entries % 10000 == 0)
                cout << "C++ Sparse read " << entries << " K3sparse entries (" << (int)(100.0*((double)entries/(double)K3sparse_entries)) << "%)\n" << std::flush;
        }
        K3sparse.makeCompressed();
        cout << "C++ Sparse Kupd_sparse about to start multiply..." << "\n" << std::flush;

        for (uint i=0; i<5; i++) {
            startTime = getWallTime();

            SparseMatrix<double> Kupd_sparse = K2 * K3sparse;

            elapsedTime = getWallTime() - startTime;
            cout << "C++ Sparse Kupd_sparse multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;
        }

        process_mem_usage(vm_used);
        cout << "C++ Sparse Kupd_sparse multiply memory: " << vm_used << "\n" << std::flush;
    }

    //process_mem_usage(vm_used);
    //cout << "C++ Sparse After delete memory: " << vm_used << "\n" << std::flush;

}

