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
 #include <algorithm>
 #include <condition_variable>
 #include <cctype>
 #include <locale>
 #include "alignment_3n_table.h"
 
 using namespace std;
 
 extern bool CG_only;
 extern long long int loadingBlockSize;
 
 /**
  * store unique information for one base information with readID, and the quality.
  */
 class uniqueID {
 public:
     unsigned long long readNameID;
     bool isConverted;
     char quality;
     bool removed;
 
     uniqueID(unsigned long long InReadNameID, bool InIsConverted, char InQual)
         : readNameID(InReadNameID), isConverted(InIsConverted), quality(InQual), removed(false) {}
 };
 
 /**
  * basic class to store reference position information
  */
 class Position {
     mutex mutex_;
 public:
     string chromosome; // reference chromosome name
     long long int location; // 1-based position
     char strand; // +(REF) or -(REF-RC)
     string convertedQualities; // each char is a mapping quality on this position for converted base.
     string unconvertedQualities; // each char is a mapping quality on this position for unconverted base.
     vector<uniqueID> uniqueIDs; // each value represent a readName which contributed the base information.
 
     void initialize() {
         chromosome.clear();
         location = -1;
         strand = '?';
         convertedQualities.clear();
         unconvertedQualities.clear();
         vector<uniqueID>().swap(uniqueIDs);
     }
 
     Position() {
         initialize();
     };
 
     /**
      * return true if there is mapping information in this reference position.
      */
     bool empty() {
         return convertedQualities.empty() && unconvertedQualities.empty();
     }
 
     /**
      * set the chromosome, location (position), and strand information.
      */
     void set(string& inputChr, long long int inputLoc) {
         chromosome = inputChr;
         location = inputLoc + 1;
     }
 
     void set(char inputStrand) {
         strand = inputStrand;
     }
 
     /**
      * 以迭代方式（利用 lower_bound）取得 readNameID 插入的位置
      */
     int searchReadNameID(unsigned long long readNameID) {
         auto it = std::lower_bound(uniqueIDs.begin(), uniqueIDs.end(), readNameID,
             [](const uniqueID &u, unsigned long long id) {
                 return u.readNameID < id;
             });
         return it - uniqueIDs.begin();
     }
 
     /**
      * 將新的 readNameID 加入 uniqueIDs，若已存在則依情況處理，並回傳是否新增成功
      */
     bool appendReadNameID(PosQuality& InBase, Alignment& InAlignment) {
         int idCount = uniqueIDs.size();
         if (idCount == 0 || InAlignment.readNameID > uniqueIDs.back().readNameID) {
             uniqueIDs.emplace_back(InAlignment.readNameID, InBase.converted, InBase.qual);
             return true;
         }
         int index = searchReadNameID(InAlignment.readNameID);
         if (index < uniqueIDs.size() && uniqueIDs[index].readNameID == InAlignment.readNameID) {
             // 如果新讀取的資訊與已存在的 conversion 狀態一致，則忽略；否則刪除原有資訊
             if (uniqueIDs[index].removed) {
                 return false;
             }
             if (uniqueIDs[index].isConverted != InBase.converted) {
                 uniqueIDs[index].removed = true;
                 if (uniqueIDs[index].isConverted) {
                     for (size_t i = 0; i < convertedQualities.size(); i++) {
                         if (convertedQualities[i] == InBase.qual) {
                             convertedQualities.erase(convertedQualities.begin() + i);
                             return false;
                         }
                     }
                 } else {
                     for (size_t i = 0; i < unconvertedQualities.size(); i++) {
                         if (unconvertedQualities[i] == InBase.qual) {
                             unconvertedQualities.erase(unconvertedQualities.begin() + i);
                             return false;
                         }
                     }
                 }
             }
             return false;
         } else {
             uniqueIDs.emplace(uniqueIDs.begin() + index, InAlignment.readNameID, InBase.converted, InBase.qual);
             return true;
         }
     }
 
     /**
      * 將 SAM 中的資訊附加到此 Position 中
      */
     void appendBase(PosQuality& input, Alignment& a) {
         // 使用 lock_guard 確保 mutex 正確解鎖
         std::lock_guard<std::mutex> lock(mutex_);
         if (appendReadNameID(input, a)) {
             if (input.converted) {
                 convertedQualities.push_back(input.qual);
             } else {
                 unconvertedQualities.push_back(input.qual);
             }
         }
     }
 };
 
 /**
  * store all reference position in this class.
  */
 class Positions {
 public:
     vector<Position*> refPositions; // 當前參考序列的 pool
     string chromosome; // 當前參考染色體名稱
     long long int location; // 當前染色體上的位置（bp）
     char lastBase = 'X'; // 先前讀入的 base，用於 CG_only 模式
     SafeQueue<string*> linePool; // 存放未處理 SAM 行的 pool
     SafeQueue<string*> freeLinePool; // 可重複利用的 string* pool
     SafeQueue<Position*> freePositionPool; // 可重複利用的 Position* pool
     SafeQueue<Position*> outputPositionPool; // 待輸出的 Position* pool
     bool working;
     mutex mutex_;
     long long int refCoveredPosition; // 已載入參考序列的最後位置
     ifstream refFile;
     vector<mutex*> workerLock; // 每個 worker 執行緒各自的鎖
     int nThreads = 1;
     ChromosomeFilePositions chromosomePos; // 儲存染色體名稱與對應在檔案中的位置
     bool addedChrName = false;
     bool removedChrName = false;
 
     // 新增條件變數以改善忙等待
     condition_variable cv_freePosition;
     mutex mtx_freePosition;
     condition_variable cv_refPositions;
     mutex mtx_refPositions;
     condition_variable cv_outputPool;
     mutex mtx_outputPool;
 
     Positions(string inputRefFileName, int inputNThreads, bool inputAddedChrName, bool inputRemovedChrName) {
         working = true;
         nThreads = inputNThreads;
         addedChrName = inputAddedChrName;
         removedChrName = inputRemovedChrName;
         for (int i = 0; i < nThreads; i++) {
             workerLock.push_back(new mutex);
         }
         refFile.open(inputRefFileName, ios_base::in);
         LoadChromosomeNamesPos();
     }
 
     ~Positions() {
         for (size_t i = 0; i < workerLock.size(); i++) {
             delete workerLock[i];
         }
         Position* pos;
         while (freePositionPool.popFront(pos)) {
             delete pos;
         }
     }
 
     /**
      * 由目標位置計算在 refPositions 中的索引
      */
     int getIndex(long long int &targetPos) {
         int firstPos = refPositions[0]->location;
         return targetPos - firstPos;
     }
 
     /**
      * 由 FASTA 的 header 行（以 '>' 開頭）解析出染色體名稱
      */
     string getChrName(string& inputLine) {
         string name;
         for (size_t i = 1; i < inputLine.size(); i++) {
             char c = inputLine[i];
             if (isspace(c)) {
                 break;
             }
             name.push_back(c);
         }
 
         if (removedChrName) {
             if (name.find("chr") == 0) {
                 name = name.substr(3);
             }
         } else if (addedChrName) {
             if (name.find("chr") != 0) {
                 name = string("chr") + name;
             }
         }
         return name;
     }
 
     /**
      * 掃描參考檔案，記錄每條染色體的檔案位置
      */
     void LoadChromosomeNamesPos() {
         string line;
         while (getline(refFile, line)) {
             if (!line.empty() && line.front() == '>') { // 此行為染色體名稱
                 chromosome = getChrName(line);
                 streampos currentPos = refFile.tellg();
                 chromosomePos.append(chromosome, currentPos);
             }
         }
         chromosomePos.sort();
         chromosome.clear();
         // 將檔案指標回復開頭，方便後續依染色體載入
         refFile.clear();
         refFile.seekg(0, ios::beg);
     }
 
     /**
      * 將 FASTA 序列中的一行（非 header）解析成多個 Position
      */
     void appendRefPosition(string& line) {
         Position* newPos;
         // 逐一檢查每個 base
         for (size_t i = 0; i < line.size(); i++) {
             getFreePosition(newPos);
             newPos->set(chromosome, location + i);
             char b = line[i];
             if (CG_only) {
                 if (lastBase == 'C' && b == 'G') {
                     // 將前一個位置設定為 '+'，而本位置設定為 '-'
                     if (!refPositions.empty())
                         refPositions.back()->set('+');
                     newPos->set('-');
                 }
             } else {
                 if (b == convertFrom) {
                     newPos->set('+');
                 } else if (b == convertFromComplement) {
                     newPos->set('-');
                 }
             }
             {
                 // 將新增的 Position 推入 refPositions，並通知等待者
                 std::lock_guard<std::mutex> lock(mtx_refPositions);
                 refPositions.push_back(newPos);
             }
             cv_refPositions.notify_all();
             lastBase = b;
         }
         location += line.size();
     }
 
     /**
      * 若能通過所有 workerLock，代表暫無 worker 在寫入
      */
     void appendingFinished() {
         for (int i = 0; i < nThreads; i++) {
             workerLock[i]->lock();
             workerLock[i]->unlock();
         }
     }
 
     /**
      * 輸出執行緒的工作函式：將 outputPositionPool 內的 Position 輸出
      */
     void outputFunction(string outputFileName) {
         ostream* out_ = &cout;
         ofstream tableFile;
         if (!outputFileName.empty()) {
             tableFile.open(outputFileName, ios_base::out);
             out_ = &tableFile;
         }
 
         *out_ << "ref\tpos\tstrand\tconvertedBaseQualities\tconvertedBaseCount\tunconvertedBaseQualities\tunconvertedBaseCount\n";
         Position* pos;
         while (working) {
             if (outputPositionPool.popFront(pos)) {
                 *out_ << pos->chromosome << '\t'
                       << to_string(pos->location) << '\t'
                       << pos->strand << '\t'
                       << pos->convertedQualities << '\t'
                       << to_string(pos->convertedQualities.size()) << '\t'
                       << pos->unconvertedQualities << '\t'
                       << to_string(pos->unconvertedQualities.size()) << '\n';
                 returnPosition(pos);
                 // 當有新 Position 可用時，通知等待中的 getFreePosition
                 cv_freePosition.notify_one();
             } else {
                 // 使用條件變數等待，降低 busy-wait CPU 佔用
                 unique_lock<mutex> lock(mtx_outputPool);
                 cv_outputPool.wait_for(lock, std::chrono::microseconds(1));
             }
         }
         tableFile.close();
     }
 
     /**
      * 將位置值小於 (refCoveredPosition - loadingBlockSize) 的 Position 移入 output pool
      */
     void moveBlockToOutput() {
         if (refPositions.empty()) {
             return;
         }
         size_t index = 0;
         for (; index < refPositions.size(); index++) {
             if (refPositions[index]->location < refCoveredPosition - loadingBlockSize) {
                 if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                     returnPosition(refPositions[index]);
                 } else {
                     outputPositionPool.push(refPositions[index]);
                     // 通知 outputFunction 有新資料
                     cv_outputPool.notify_one();
                 }
             } else {
                 break;
             }
         }
         if (index != 0) {
             refPositions.erase(refPositions.begin(), refPositions.begin() + index);
         }
     }
 
     /**
      * 將所有 refPositions 移入 output pool
      */
     void moveAllToOutput() {
         if (refPositions.empty()) {
             return;
         }
         for (size_t index = 0; index < refPositions.size(); index++) {
             if (refPositions[index]->empty() || refPositions[index]->strand == '?') {
                 returnPosition(refPositions[index]);
             } else {
                 vector<uniqueID>().swap(refPositions[index]->uniqueIDs);
                 outputPositionPool.push(refPositions[index]);
                 cv_outputPool.notify_one();
             }
         }
         refPositions.clear();
     }
 
     /**
      * 載入新染色體的參考序列（初始載入約 2*loadingBlockSize bp）
      */
     void loadNewChromosome(string targetChromosome) {
         refFile.clear();
         // 根據染色體名稱找到檔案中的起始位置
         streampos startPos = chromosomePos.getChromosomePosInRefFile(targetChromosome);
         chromosome = targetChromosome;
         refFile.seekg(startPos, ios::beg);
         refCoveredPosition = 2 * loadingBlockSize;
         string line;
         lastBase = 'X';
         location = 0;
         while (getline(refFile, line)) {
             if (!line.empty() && line.front() == '>') {
                 // 遇到下一條染色體，直接返回
                 return;
             } else {
                 if (line.empty()) { continue; }
                 // 轉換為大寫（結果與原本完全一致）
                 std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                 appendRefPosition(line);
                 if (location >= refCoveredPosition) {
                     return;
                 }
             }
         }
     }
 
     /**
      * 載入更多參考序列（每次約 loadingBlockSize bp），若遇到下一條染色體則返回
      */
     void loadMore() {
         refCoveredPosition += loadingBlockSize;
         string line;
         while (getline(refFile, line)) {
             if (!line.empty() && line.front() == '>') {
                 return;
             } else {
                 if (line.empty()) { continue; }
                 std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                 appendRefPosition(line);
                 if (location >= refCoveredPosition) {
                     return;
                 }
             }
         }
     }
 
     /**
      * 將 Alignment 的資訊加入對應的參考位置
      */
     void appendPositions(Alignment& newAlignment) {
         if (!newAlignment.mapped || newAlignment.bases.empty()) {
             return;
         }
         long long int startPos = newAlignment.location; // 1-based
         int index = getIndex(newAlignment.location);
 
         for (size_t i = 0; i < newAlignment.sequence.size(); i++) {
             PosQuality* b = &newAlignment.bases[i];
             if (b->remove) {
                 continue;
             }
             Position* pos = refPositions[index + b->refPos];
             assert(pos->location == startPos + b->refPos);
             if (pos->strand == '?') {
                 continue;
             }
             pos->appendBase(newAlignment.bases[i], newAlignment);
         }
     }
 
     /**
      * 從 freeLinePool 取得一個 string*，若沒有則新建立
      */
     void getFreeStringPointer(string*& newLine) {
         if (freeLinePool.popFront(newLine)) {
             return;
         } else {
             newLine = new string();
         }
     }
 
     /**
      * 從 freePositionPool 取得一個 Position*，若沒有則新建立
      * 利用條件變數避免忙等待
      */
     void getFreePosition(Position*& newPosition) {
         {
             unique_lock<mutex> lock(mtx_freePosition);
             cv_freePosition.wait(lock, [this]() { return outputPositionPool.size() < 10000; });
         }
         if (freePositionPool.popFront(newPosition)) {
             return;
         } else {
             newPosition = new Position();
         }
     }
 
     /**
      * 將 string 回收至 freeLinePool
      */
     void returnLine(string* line) {
         line->clear();
         freeLinePool.push(line);
     }
 
     /**
      * 將 Position 回收至 freePositionPool
      */
     void returnPosition(Position* pos) {
         pos->initialize();
         freePositionPool.push(pos);
         cv_freePosition.notify_one();
     }
 
     /**
      * worker 執行緒的主要工作函式：
      * 從 linePool 取出 SAM 行、解析、並將資訊加入參考位置
      */
     void append(int threadID) {
         string* line;
         Alignment newAlignment;
 
         while (working) {
             {
                 // 使用各自的 workerLock 保護
                 std::lock_guard<std::mutex> lock(*workerLock[threadID]);
                 if (!linePool.popFront(line)) {
                     // 若沒有資料，則稍待並繼續
                     this_thread::sleep_for(std::chrono::nanoseconds(1));
                     continue;
                 }
             }
             {
                 // 等待 refPositions 有資料可用
                 unique_lock<mutex> lock(mtx_refPositions);
                 cv_refPositions.wait(lock, [this]() { return !refPositions.empty(); });
             }
             newAlignment.parse(line);
             returnLine(line);
             appendPositions(newAlignment);
         }
     }
 };
 
 #endif //POSITION_3N_TABLE_H
 