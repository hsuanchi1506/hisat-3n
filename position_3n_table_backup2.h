/*
* Copyright 2020, Yun (Leo) Zhang <imzhangyun@gmail.com>
*
* This file is part of HISAT-3N.
*
* HISAT-3N is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* HISAT-3N is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with HISAT-3N.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef POSITION_3N_TABLE_H
#define POSITION_3N_TABLE_H

#include <string>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>
#include <cassert>
#include "alignment_3n_table.h"

// Add mmap related headers
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

extern bool CG_only;
extern long long int loadingBlockSize;

/**
 * store unique information for one base information with readID, and the quality.
 */
class uniqueID
{
public:
    unsigned long long readNameID;
    bool isConverted;
    char quality;
    bool removed;

    uniqueID(unsigned long long InReadNameID,
            bool InIsConverted,
            char& InQual){
        readNameID = InReadNameID;
        isConverted = InIsConverted;
        quality = InQual;
        removed = false;
    }
};

/**
 * basic class to store reference position information
 */
class Position{
    mutex mutex_;
public:
    string chromosome; // reference chromosome name
    long long int location; // 1-based position
    char strand; // +(REF) or -(REF-RC)
    string convertedQualities; // each char is a mapping quality on this position for converted base.
    string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.
    vector<uniqueID> uniqueIDs; // each value represent a readName which contributed the base information.
                                // readNameIDs is to make sure no read contribute 2 times in same position.

    void initialize() {
        chromosome.clear();
        location = -1;
        strand = '?';
        convertedQualities.clear();
        unconvertedQualities.clear();
        vector<uniqueID>().swap(uniqueIDs);
    }

    Position(){
        initialize();
    };

    bool empty() {
        return convertedQualities.empty() && unconvertedQualities.empty();
    }

    void set (string& inputChr, long long int inputLoc) {
        chromosome = inputChr;
        location = inputLoc + 1;
    }

    void set(char inputStrand) {
        strand = inputStrand;
    }

    // [UPDATE] Iteratively
    int searchReadNameID(unsigned long long& readNameID, int start, int end) {
        while (start <= end) {
            int middle = start + (end - start) / 2;
            if (uniqueIDs[middle].readNameID == readNameID) {
                return middle;
            }
            if (uniqueIDs[middle].readNameID > readNameID) {
                end = middle - 1;
            } else {
                start = middle + 1;
            }
        }
        return start;
    }

    bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
        int idCount = uniqueIDs.size();
        if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
            uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
        int index = searchReadNameID(InAlignment.readNameID, 0, idCount);
        if (uniqueIDs[index].readNameID == InAlignment.readNameID) {
            if (uniqueIDs[index].removed) {
                return false;
            }
            if (uniqueIDs[index].isConverted != InBase.converted) {
                uniqueIDs[index].removed = true;
                if (uniqueIDs[index].isConverted) {
                    for (int i = 0; i < convertedQualities.size(); i++) {
                        if (convertedQualities[i] == InBase.qual) {
                            convertedQualities.erase(convertedQualities.begin()+i);
                            return false;
                        }
                    }
                } else {
                    for (int i = 0; i < unconvertedQualities.size(); i++) {
                        if (unconvertedQualities[i] == InBase.qual) {
                            unconvertedQualities.erase(unconvertedQualities.begin()+i);
                            return false;
                        }
                    }
                }
            }
            return false;
        } else {
            uniqueIDs.emplace(uniqueIDs.begin()+index, InAlignment.readNameID, InBase.converted, InBase.qual);
            return true;
        }
    }

    void appendBase (PosQuality& input, Alignment& a) {
        mutex_.lock();
        if (appendReadNameID(input,a)) {
            if (input.converted) {
                convertedQualities += input.qual;
            } else {
                unconvertedQualities += input.qual;
            }
        }
        mutex_.unlock();
    }
};

/**
 * store all reference position in this class.
 * The following changes focus on reading reference genome data using mmap to accelerate data access.
 */
class Positions{
public:
    vector<Position*> refPositions;
    string chromosome;
    long long int location;
    char lastBase = 'X';
    SafeQueue<string*> linePool;
    SafeQueue<string*> freeLinePool;
    SafeQueue<Position*> freePositionPool;
    SafeQueue<Position*> outputPositionPool;
    bool working;
    mutex mutex_;
    long long int refCoveredPosition;
    // Replaced ifstream-based reference file reading with mmap
    // ifstream refFile;
    // Add members required for mmap
    char *refData;
    size_t refSize;
    size_t refOffset;

    vector<mutex*> workerLock;
    int nThreads = 1;
    ChromosomeFilePositions chromosomePos;
    bool addedChrName = false;
    bool removedChrName = false;

    // Modified constructor of Positions to open the reference file using mmap
    Positions(string inputRefFileName, int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName) {
        working = true;
        nThreads = inputNThreads;
        addedChrName = inputAddedChrName;
        removedChrName = inputRemovedChrName;
        for (int i = 0; i < nThreads; i++) {
            workerLock.push_back(new mutex);
        }
        // Load reference file using open/fstat/mmap
        int fd = open(inputRefFileName.c_str(), O_RDONLY);
        if (fd < 0) {
            perror("open");
            exit(1);
        }
        struct stat sb;
        if (fstat(fd, &sb) == -1) {
            perror("fstat");
            exit(1);
        }
        refSize = sb.st_size;
        refData = static_cast<char*>(mmap(NULL, refSize, PROT_READ, MAP_PRIVATE, fd, 0));
        if (refData == MAP_FAILED) {
            perror("mmap");
            exit(1);
        }
        close(fd);
        refOffset = 0;
        // Initialize chromosome position information from mmap data
        LoadChromosomeNamesPos();
    }

    ~Positions() {
        for (int i = 0; i < workerLock.size(); i++) {
            delete workerLock[i];
        }
        Position* pos;
        while(freePositionPool.popFront(pos)) {
            delete pos;
        }
        // Unmap mmap
        munmap(refData, refSize);
    }

    int getIndex(long long int &targetPos) {
        int firstPos = refPositions[0]->location;
        return targetPos - firstPos;
    }

    string getChrName(string& inputLine) {
        string name;
        for (int i = 1; i < inputLine.size(); i++)
        {
            char c = inputLine[i];
            if (isspace(c)){
                break;
            }
            name += c;
        }

        if(removedChrName) {
            if(name.find("chr") == 0) {
                name = name.substr(3);
            }
        } else if(addedChrName) {
            if(name.find("chr") != 0) {
                name = string("chr") + name;
            }
        }
        return name;
    }

    // Scan the entire reference file using mmap data to record the start position of each chromosome
    void LoadChromosomeNamesPos() {
        size_t tempOffset = 0;
        while (tempOffset < refSize) {
            size_t lineStart = tempOffset;
            while (tempOffset < refSize && refData[tempOffset] != '\n') {
                tempOffset++;
            }
            string line(refData + lineStart, refData + tempOffset);
            if (tempOffset < refSize && refData[tempOffset] == '\n') {
                tempOffset++;
            }
            if (!line.empty() && line.front() == '>') {
                string chrName = getChrName(line);
                // Use current offset as the position in file
                streampos currentPos = tempOffset;
                chromosomePos.append(chrName, currentPos);
            }
        }
        chromosomePos.sort();
        // Note: Do not change refOffset here since subsequent reads may start from the beginning or a specific location
    }

    // Helper function: read next line from current refData/offset
    string getNextLine() {
        if (refOffset >= refSize) return "";
        size_t lineStart = refOffset;
        while (refOffset < refSize && refData[refOffset] != '\n') {
            refOffset++;
        }
        string line(refData + lineStart, refData + refOffset);
        if (refOffset < refSize && refData[refOffset] == '\n')
            refOffset++;
        return line;
    }

    // Modified version of loadNewChromosome using mmap
    void loadNewChromosome(string targetChromosome) {
        streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
        refOffset = startPos;
        chromosome = targetChromosome;
        refCoveredPosition = 2 * loadingBlockSize;
        string line;
        lastBase = 'X';
        location = 0;
        while (refOffset < refSize) {
            line = getNextLine();
            if (line.empty()) { continue; }
            if (line.front() == '>') { // Reached next chromosome header
                break;
            } else {
                // Convert line to uppercase
                for (auto &c : line)
                    c = toupper(c);
                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
    }

    // Load more reference content via mmap
    void loadMore() {
        refCoveredPosition += loadingBlockSize;
        string line;
        while (refOffset < refSize) {
            line = getNextLine();
            if (line.empty()) { continue; }
            if (line.front() == '>') { // Next chromosome
                break;
            } else {
                for (auto &c : line)
                    c = toupper(c);
                appendRefPosition(line);
                if (location >= refCoveredPosition) {
                    return;
                }
            }
        }
    }

    void appendRefPosition(string& line) {
        Position* newPos;
        char* b;
        for (int i = 0; i < line.size(); i++) {
            getFreePosition(newPos);
            newPos->set(chromosome, location+i);
            b = &line[i];
            if (CG_only) {
                if (lastBase == 'C' && *b == 'G') {
                    refPositions.back()->set('+');
                    newPos->set('-');
                }
            } else {
                if (*b == convertFrom) {
                    newPos->set('+');
                } else if (*b == convertFromComplement) {
                    newPos->set('-');
                }
            }
            refPositions.push_back(newPos);
            lastBase = *b;
        }
        location += line.size();
    }

    void appendingFinished() {
        for (int i = 0; i < nThreads; i++) {
            workerLock[i]->lock();
            workerLock[i]->unlock();
        }
    }

    void outputFunction(string outputFileName) {
        ofstream tableFile;
        ostream* out_ = &cout;
        
        if (!outputFileName.empty()) {
            tableFile.open(outputFileName, ios_base::out | ios_base::binary); // 使用 binary 模式加快寫入
            out_ = &tableFile;
        }
        
        string outputBuffer;
        outputBuffer.reserve(1024 * 1024); // 1MB buffer
        
        *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
        
        Position* pos;
        while (working || !outputPositionPool.empty()) {
            if (outputPositionPool.popFront(pos)) {
                outputBuffer.clear();
                outputBuffer.append(pos->chromosome)
                            .append("\t")
                            .append(to_string(pos->location))
                            .append("\t")
                            .append(1, pos->strand)
                            .append("\t")
                            .append(pos->convertedQualities)
                            .append("\t")
                            .append(to_string(pos->convertedQualities.size()))
                            .append("\t")
                            .append(pos->unconvertedQualities)
                            .append("\t")
                            .append(to_string(pos->unconvertedQualities.size()))
                            .append("\n");
                
                out_->write(outputBuffer.data(), outputBuffer.size());
                returnPosition(pos);
                
                // 定期 flush
                if (outputBuffer.size() > 1024 * 512) { // 512KB
                    out_->flush();
                }
            } else {
                this_thread::sleep_for(chrono::microseconds(1));
            }
        }
        
        if (tableFile.is_open()) {
            tableFile.close();
        }
    }


    void moveBlockToOutput() {
        if (refPositions.empty()) {
            return;
        }
        int index;
        for (index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
                if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                    returnPosition(refPositions[index]);
                } else {
                    outputPositionPool.push(refPositions[index]);
                }
            } else {
                break;
            }
        }
        if (index != 0) {
            refPositions.erase(refPositions.begin(), refPositions.begin()+index);
        }
    }

    void moveAllToOutput() {
        if (refPositions.empty()) {
            return;
        }
        for (int index = 0; index < refPositions.size(); index++) {
            if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                returnPosition(refPositions[index]);
            } else {
                vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
                outputPositionPool.push(refPositions[index]);
            }
        }
        refPositions.clear();
    }

    void getFreeStringPointer(string*& newLine) {
        if (freeLinePool.popFront(newLine)) {
            return;
        } else {
            newLine = new string();
        }
    }

    void getFreePosition(Position*& newPosition) {
        while (outputPositionPool.size() >= 10000) {
            this_thread::sleep_for (std::chrono::microseconds(1));
        }
        if (freePositionPool.popFront(newPosition)) {
            return;
        } else {
            newPosition = new Position();
        }
    }

    void returnLine(string* line) {
        line->clear();
        freeLinePool.push(line);
    }

    void returnPosition(Position* pos) {
        pos->initialize();
        freePositionPool.push(pos);
    }

    void append(int threadID) {
        string* line;
        Alignment newAlignment;
        

        while (working) {
            workerLock[threadID]->lock();
            if(!linePool.popFront(line)) {
                workerLock[threadID]->unlock();
                this_thread::sleep_for (std::chrono::nanoseconds(1));
                continue;
            }
            while (refPositions.empty()) {
                this_thread::sleep_for (std::chrono::microseconds(1));
            }
            newAlignment.parse(line);
            returnLine(line);
            appendPositions(newAlignment);
            workerLock[threadID]->unlock();
        }
    }

    void appendPositions(Alignment& newAlignment) {
        if (!newAlignment.mapped || newAlignment.bases.empty()) {
            return;
        }
        long long int startPos = newAlignment.location;
        int index = getIndex(newAlignment.location);

        for (int i = 0; i < newAlignment.sequence.size(); i++) {
            PosQuality* b = &newAlignment.bases[i];
            if (b->remove) {
                continue;
            }

            Position* pos = refPositions[index+b->refPos];
            assert (pos->location == startPos + b->refPos);

            if (pos->strand == '?') {
                continue;
            }
            pos->appendBase(newAlignment.bases[i], newAlignment);
        }
    }
};

#endif //POSITION_3N_TABLE_H
