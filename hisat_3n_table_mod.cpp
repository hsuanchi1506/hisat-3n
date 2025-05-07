#include <iostream>
#include <getopt.h>
#include <fstream>
#include <thread>
#include <vector>
#include <chrono>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "position_3n_table.h"

using namespace std;

string alignmentFileName;
string refFileName;
string outputFileName;
bool standardInMode = false;
bool uniqueOnly = false;
bool multipleOnly = false;
bool CG_only = false;
int nThreads = 1;
long long int loadingBlockSize = 1000000;
char convertFrom = '0';
char convertTo = '0';
char convertFromComplement;
char convertToComplement;
bool addedChrName = false;
bool removedChrName = false;

Positions *positions;

bool fileExist(string &filename)
{
    ifstream file(filename);
    return file.good();
}

enum
{
    ARG_ADDED_CHRNAME = 256,
    ARG_REMOVED_CHRNAME
};

static const char *short_options = "p:umr:a:b:o:h";
static struct option long_options[] = {
    {"alignments", required_argument, 0, 'a'},
    {"ref", required_argument, 0, 'r'},
    {"output-name", required_argument, 0, 'o'},
    {"base-change", required_argument, 0, 'b'},
    {"unique-only", no_argument, 0, 'u'},
    {"multiple-only", no_argument, 0, 'm'},
    {"threads", required_argument, 0, 'p'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}};

static void printHelp(ostream &out)
{
    out << "hisat-3n-table developed by Yun (Leo) Zhang" << endl;
    out << "Usage:" << endl
        << "hisat-3n-table [options]* --alignments <alignmentFile> --ref <refFile> --output-name <outputFile> --base-change <char1,char2>" << endl
        << "  <alignmentFile>           SORTED SAM filename. Please enter '-' for standard input." << endl
        << "  <refFile>                 reference file (should be FASTA format)." << endl
        << "  <outputFile>              file name to save the 3n table (tsv format)." << endl
        << "  <chr1,chr2>               the char1 is the nucleotide converted from, the char2 is the nucleotide converted to." << endl;
    out << "Options (defaults in parentheses):" << endl
        << "  -u/--unique-only          only count the base which is in unique mapped reads." << endl
        << "  -m/--multiple-only        only count the base which is in multiple mapped reads." << endl
        << "  -p/--threads <int>        number of threads to launch (1)." << endl
        << "  -h/--help                 print this usage message." << endl;
}

static void parseOption(int next_option, const char *optarg)
{
    switch (next_option)
    {
    case 'a':
        alignmentFileName = optarg;
        if (alignmentFileName == "-")
        {
            standardInMode = true;
        }
        else if (!fileExist(alignmentFileName))
        {
            cerr << "The alignment file does not exist." << endl;
            throw(1);
        }
        break;
    case 'r':
        refFileName = optarg;
        if (!fileExist(refFileName))
        {
            cerr << "Reference (FASTA) file does not exist." << endl;
            throw(1);
        }
        break;
    case 'o':
        outputFileName = optarg;
        break;
    case 'b':
    {
        string arg = optarg;
        if (arg.size() != 3 || arg[1] != ',')
        {
            cerr << "Error: expected 2 comma-separated arguments to --base-change (e.g. C,T)." << endl;
            throw 1;
        }
        convertFrom = toupper(arg.front());
        convertTo = toupper(arg.back());
        break;
    }
    case 'u':
        uniqueOnly = true;
        break;
    case 'm':
        multipleOnly = true;
        break;
    case 'p':
        nThreads = stoi(optarg);
        if (nThreads < 1)
            nThreads = 1;
        break;
    case 'h':
        printHelp(cerr);
        throw 0;
    default:
        printHelp(cerr);
        throw 1;
    }
}

static void parseOptions(int argc, const char **argv)
{
    int option_index = 0;
    int next_option;
    while ((next_option = getopt_long(argc, const_cast<char **>(argv), short_options, long_options, &option_index)) != -1)
        parseOption(next_option, optarg);

    if (refFileName.empty() || alignmentFileName.empty())
    {
        cerr << "No reference or SAM file specified!" << endl;
        printHelp(cerr);
        throw 1;
    }
    if (convertFrom == '0' || convertTo == '0')
    {
        cerr << "The --base-change argument is required." << endl;
        throw 1;
    }
    convertFromComplement = 'G';
    convertToComplement = 'A';
}

bool getSAMChromosomePos(string *line, string &chr, long long int &pos)
{
    size_t start = 0;
    int field = 0;
    while (field < 4)
    {
        size_t tab = line->find('\t', start);
        if (tab == string::npos)
            return false;
        if (field == 2)
            chr = line->substr(start, tab - start);
        else if (field == 3)
            pos = stoll(line->substr(start, tab - start));
        start = tab + 1;
        field++;
    }
    return chr != "*";
}

int hisat_3n_table()
{
    positions = new Positions(refFileName, nThreads, addedChrName, removedChrName);

    vector<thread *> workers;
    for (int i = 0; i < nThreads; i++)
        workers.push_back(new thread(&Positions::append, positions, i));

    thread outputThread(&Positions::outputFunction, positions, outputFileName);

    if (standardInMode)
    {
        cerr << "Standard input is not supported with memory-mapped mode." << endl;
        exit(EXIT_FAILURE);
    }

    int fd = open(alignmentFileName.c_str(), O_RDONLY);
    if (fd == -1)
    {
        perror("open");
        exit(EXIT_FAILURE);
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1)
    {
        perror("fstat");
        close(fd);
        exit(EXIT_FAILURE);
    }
    size_t fileSize = sb.st_size;

    char *data = static_cast<char *>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (data == MAP_FAILED)
    {
        perror("mmap");
        close(fd);
        exit(EXIT_FAILURE);
    }

    char *ptr = data;
    char *end = data + fileSize;
    string samChromosome;
    long long int samPos, reloadPos = 0, lastPos = 0;

    while (ptr < end)
    {
        char *lineStart = ptr;
        while (ptr < end && *ptr != '\n')
            ptr++;
        string lineStr(lineStart, ptr - lineStart);
        ptr++;

        if (lineStr.empty() || lineStr[0] == '@')
            continue;

        string *line;
        positions->getFreeStringPointer(line);
        *line = lineStr;

        if (!getSAMChromosomePos(line, samChromosome, samPos))
        {
            positions->returnLine(line);
            continue;
        }

        while (positions->linePool.size() > 1000 * nThreads)
            this_thread::sleep_for(std::chrono::microseconds(1));

        if (samChromosome != positions->chromosome)
        {
            while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000)
                this_thread::sleep_for(std::chrono::microseconds(1));
            positions->appendingFinished();
            positions->moveAllToOutput();
            positions->loadNewChromosome(samChromosome);
            reloadPos = loadingBlockSize;
            lastPos = 0;
        }

        while (samPos > reloadPos)
        {
            while (!positions->linePool.empty() || positions->outputPositionPool.size() > 100000)
                this_thread::sleep_for(std::chrono::microseconds(1));
            positions->appendingFinished();
            positions->moveBlockToOutput();
            positions->loadMore();
            reloadPos += loadingBlockSize;
        }

        if (lastPos > samPos)
        {
            cerr << "Input alignment file is not sorted." << endl;
            throw 1;
        }

        positions->linePool.push(line);
        lastPos = samPos;
    }

    if (munmap(data, fileSize) == -1)
        perror("munmap");
    close(fd);

    while (!positions->linePool.empty())
        this_thread::sleep_for(std::chrono::microseconds(100));

    positions->appendingFinished();
    positions->moveAllToOutput();
    while (!positions->outputPositionPool.empty())
        this_thread::sleep_for(std::chrono::microseconds(100));

    string *line;
    while (positions->freeLinePool.popFront(line))
        delete line;

    positions->working = false;
    for (auto &worker : workers)
    {
        worker->join();
        delete worker;
    }
    outputThread.join();

    delete positions;
    return 0;
}

int main(int argc, const char **argv)
{
    int ret = 0;
    try
    {
        parseOptions(argc, argv);
        ret = hisat_3n_table();
    }
    catch (std::exception &e)
    {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    catch (int e)
    {
        return e;
    }
    return ret;
}