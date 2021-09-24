
#include <string>
using std::string;

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
#include <sstream>

#include "cs.h"
#include "klu.h"

#include <unistd.h>
#include <sys/time.h>

#include <unordered_map>
using std::unordered_map;


void stdout_cs_compress(cs *&a, const uint &precision=16) {
    // First copy into a map
    unordered_map<uint,unordered_map<uint,double>> mat;
    for ( uint i = 0 ; i < a->n ; i++ ) {
        for ( uint j = a->p[i] ; j < a->p[i+1] ; j++ ) {
            mat[a->i[j]][i] = a->x[j];
        }
    }
    // write to stdout
    cout << "printing CS compressed matrix:\n" << std::flush;
    for ( uint i = 0 ; i < a->m ; i++ )
        for ( uint j = 0 ; j < a->n ; j++ )
            cout << mat[i][j] << ( j == a->n-1 ? "\n" : "," );
    cout << "done printing CS compressed matrix\n" << std::flush;
}



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

    double startTime, elapsedTime;
    string vm_used;

#if 100
    cs* SupdRaw = cs_spalloc(zqty, zqty, Supd_entries, 1, 1);
    std::ifstream SupdFile("../mat/Supd_trip.csv");
    while (std::getline(SupdFile, line)) {
        std::istringstream stream(line);
        std::getline(stream, istr, ','); i = std::stoi(istr);
        std::getline(stream, jstr, ','); j = std::stoi(jstr);
        std::getline(stream, valstr, ','); val = std::stod(valstr);
        cs_entry(SupdRaw, i, j, val);
    }
    cs* Supd = cs_compress(SupdRaw); cs_spfree(SupdRaw);

    double *rhs = (double *)calloc(zqty*zqty, sizeof(double));

    try {
        for (uint i=0; i<5; i++) {
            startTime = getWallTime();

            // Initialize klusolve variables
            klu_symbolic *klusym;
            klu_numeric *klunum;
            klu_common klucom;
            if (!klu_defaults(&klucom)) throw "klu_defaults failed";

            klusym = klu_analyze(Supd->m,Supd->p,Supd->i,&klucom);
            if (!klusym) throw "klu_analyze failed";

            klunum = klu_factor(Supd->p,Supd->i,Supd->x,klusym,&klucom);
            if (!klunum) {
                cout << "Common->status is: " << klucom.status << "\n" << std::flush;
                if ( klucom.status == 1 ) cout << "\tKLU_SINGULAR\n" << std::flush;
                throw "klu_factor failed";
            }

            // KLU condition number estimation
            //(void)klu_condest(Supd->p,Supd->x,klusym,klunum,&klucom);
            //cout << "klu_condest Supd condition number estimate: " << klucom.condest << "\n" << std::flush;

            // initialize an identity right-hand side
            // assumes all non-diagonal entries are 0 through calloc call
            for ( uint ii = 0 ; ii < zqty ; ii++ )
                rhs[ii*zqty + ii] = 1.0;

            klu_solve(klusym,klunum,Supd->m,Supd->n,rhs,&klucom);
            if (klucom.status) {
                cout << "Common->status is: " << klucom.status << "\n" << std::flush;
                throw "klu_solve failed";
            }

            elapsedTime = getWallTime() - startTime;
            cout << "C++ SuiteSparse Supd inverse time: " << getMinSec(elapsedTime) << "\n" << std::flush;

            // free klusym and klunum or memory leak results
            klu_free_symbolic(&klusym, &klucom);
            klu_free_numeric(&klunum, &klucom);
        }

        process_mem_usage(vm_used);
        cout << "C++ SuiteSparse Supd inverse memory: " << vm_used << "\n" << std::flush;
    } catch (const char *msg) {
        cout << "KLU ERROR: " << msg << "\n" << std::flush;
        throw "klu_error";
    }
    cs_spfree(Supd);
    free(rhs);

    //process_mem_usage(vm_used);
    //cout << "C++ SuiteSparse After delete memory: " << vm_used << "\n" << std::flush;
#endif

    cs* K2Raw = cs_spalloc(xqty, zqty, K2_entries, 1, 1);
    std::ifstream K2file("../mat/K2_trip.csv");
    while (std::getline(K2file, line)) {
        std::istringstream stream(line);
        std::getline(stream, istr, ','); i = std::stoi(istr);
        std::getline(stream, jstr, ','); j = std::stoi(jstr);
        std::getline(stream, valstr, ','); val = std::stod(valstr);
        //cout << "K2 i: " << i << ", j: " << j << ", val: " << val << endl;
        cs_entry(K2Raw, i, j, val);
    }
    cs* K2 = cs_compress(K2Raw); cs_spfree(K2Raw);

    cs* K3Raw = cs_spalloc(zqty, zqty, K3_entries, 1, 1);
    std::ifstream K3file("../mat/K3_trip.csv");
    while (std::getline(K3file, line)) {
        std::istringstream stream(line);
        std::getline(stream, istr, ','); i = std::stoi(istr);
        std::getline(stream, jstr, ','); j = std::stoi(jstr);
        std::getline(stream, valstr, ','); val = std::stod(valstr);
        //cout << "K3 i: " << i << ", j: " << j << ", val: " << val << endl;
        cs_entry(K3Raw, i, j, val);
    }
    cs* K3 = cs_compress(K3Raw); cs_spfree(K3Raw);

    for (uint i=0; i<5; i++) {
        startTime = getWallTime();

        cs *Kupd = cs_multiply(K2,K3);

        elapsedTime = getWallTime() - startTime;
        cout << "C++ SuiteSparse Kupd multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;

        cs_spfree(Kupd);
    }

    process_mem_usage(vm_used);
    cout << "C++ SuiteSparse Kupd multiply memory: " << vm_used << "\n" << std::flush;

    //stdout_cs_compress(Kupd);

    cs_spfree(K2); cs_spfree(K3);

    //process_mem_usage(vm_used);
    //cout << "C++ SuiteSparse After delete memory: " << vm_used << "\n" << std::flush;

    K2Raw = cs_spalloc(xqty, zqty, K2_entries, 1, 1);
    std::ifstream K2file2("../mat/K2_trip.csv");
    while (std::getline(K2file2, line)) {
        std::istringstream stream(line);
        std::getline(stream, istr, ','); i = std::stoi(istr);
        std::getline(stream, jstr, ','); j = std::stoi(jstr);
        std::getline(stream, valstr, ','); val = std::stod(valstr);
        //cout << "K2 i: " << i << ", j: " << j << ", val: " << val << endl;
        cs_entry(K2Raw, i, j, val);
    }
    K2 = cs_compress(K2Raw); cs_spfree(K2Raw);

    cs* K3sparseRaw = cs_spalloc(zqty, zqty, K3sparse_entries, 1, 1);
    std::ifstream K3file2("../mat/K3sparse_trip.csv");
    while (std::getline(K3file2, line)) {
        std::istringstream stream(line);
        std::getline(stream, istr, ','); i = std::stoi(istr);
        std::getline(stream, jstr, ','); j = std::stoi(jstr);
        std::getline(stream, valstr, ','); val = std::stod(valstr);
        //cout << "K3 i: " << i << ", j: " << j << ", val: " << val << endl;
        cs_entry(K3sparseRaw, i, j, val);
    }
    cs* K3sparse = cs_compress(K3sparseRaw); cs_spfree(K3sparseRaw);

    for (uint i=0; i<5; i++) {
        startTime = getWallTime();

        cs *Kupd_sparse = cs_multiply(K2,K3sparse);

        elapsedTime = getWallTime() - startTime;
        cout << "C++ SuiteSparse Kupd_sparse multiply time: " << getMinSec(elapsedTime) << "\n" << std::flush;

        cs_spfree(Kupd_sparse);
    }

    process_mem_usage(vm_used);
    cout << "C++ SuiteSparse Kupd_sparse multiply memory: " << vm_used << "\n" << std::flush;

    //stdout_cs_compress(Kupd_sparse);

    cs_spfree(K2); cs_spfree(K3sparse);

    //process_mem_usage(vm_used);
    //cout << "C++ SuiteSparse After delete memory: " << vm_used << "\n" << std::flush;

}

