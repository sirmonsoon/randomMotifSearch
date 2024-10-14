#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <climits>


using namespace std;

//Given Motifs we can construct Profile(Motifs)
//Using pseudocounts to avoid 0 probability option
vector<vector<double>> createProfile(const vector<string>& motifs, int k){
    vector<vector<double>> profile(4, vector<double>(k, 1.0));
    int t = motifs.size();

    for(int i = 0; i < t; i++){
        for(int j = 0; j < k; j++){
            char nucleotide = motifs[i][j];
            if(nucleotide == 'A') 
            profile[0][j]++;
            else if(nucleotide == 'C')
            profile[1][j]++;
            else if(nucleotide == 'G')
            profile[2][j]++;
            else if(nucleotide == 'T')
            profile[3][j]++;
        }
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < k; j++){
            profile[i][j] /= (t + 4.0);
        }
    }
    return profile;
}

//Scoring Motifs
int scoreMotifs(const vector<string>& motifs, int k){
    int score = 0;
    for(int j = 0; j < k; j++){
        vector<int> count(4, 0);
        for(int i = 0; i < motifs.size(); i++){
            char nucleotide = motifs[i][j];
            if(nucleotide == 'A')
            count[0]++;
            else if(nucleotide == 'C')
            count[1]++;
            else if(nucleotide == 'G')
            count[2]++;
            else if(nucleotide == 'T')
            count[3]++;
        }
        score += motifs.size() - *max_element(count.begin(), count.end());
    }
    return score;
}

//Randomly select kmers from the given sequence 
vector<string> randomlySelectedKmers(const vector<string>& dna, int k){
    vector<string> motifs;
    srand(time(0));
    for(const string& seq : dna){
        int randIndex = rand() % (seq.size() - k + 1);
        motifs.push_back(seq.substr(randIndex, k));
    }
    return motifs;
}

//Given profile -> We find most probable kmer
string findMostProbableKmer(const string& sequence, const vector<vector<double>>& profile, int k){
    double maxProb = -1.0;
    string bestKmer;
    for(int i = 0; i <= sequence.size() - k; i++){
        string kmer = sequence.substr(i, k);
        double prob = 1.0;
        for(int j = 0; j < k; j++){
            char nucleotide = kmer[j];
            if(nucleotide == 'A')
                prob *= profile[0][j];
            else if(nucleotide == 'C')
                prob *= profile[1][j];
            else if(nucleotide == 'G')
                prob *= profile[2][j];
            else if(nucleotide == 'T')
                prob *= profile[3][j];
            }
        if(prob > maxProb){
            maxProb = prob;
            bestKmer = kmer;
        }
    }
    return bestKmer;
}

//Iteratively improves motifs based on profile matrix
vector<string> randMotifSearch(const vector<string>& dna, int k, int t){
    vector<string> bestMotifs = randomlySelectedKmers(dna, k);
    int bestScore = scoreMotifs(bestMotifs, k);
    //While forever
    while(true){
        //Profile <- Profile(Motifs)
        vector<vector<double>> profile = createProfile(bestMotifs, k);
        vector<string> newMotifs;
        for(const string& seq : dna){
            newMotifs.push_back(findMostProbableKmer(seq, profile, k));
        }
        int currScore = scoreMotifs(newMotifs, k);
        if(currScore < bestScore){
            bestMotifs = newMotifs;
            bestScore = currScore;
        } else{
            return bestMotifs;
        }
    }
}

int main(int argc, char* argv[])
{
    if(argc != 2) {
        cerr << "Usage: " << argv[0] << "<input_file>" << endl;
        return 1;
    }

    string inputFileName = argv[1];
    ifstream inputFile(inputFileName);
    if (!inputFile.is_open()){
        cerr << "File Read Error" << endl;
        return 1;
    }

    int k, t;
    inputFile >> k >> t;
    vector<string> dna(t);
    for(int i = 0; i < t; i++){
        inputFile >> dna[i];
    }
    inputFile.close();

    vector<string> bestMotifs;
    int bestScore = INT_MAX;

    //Running algorithm 1500 times.
    for(int i = 0; i < 1500; i++){
        vector<string> motifs = randMotifSearch(dna, k, t);
        int currScore = scoreMotifs(motifs, k);
        if(currScore < bestScore){
            bestMotifs = motifs;
            bestScore = currScore;
        }
    }
    
    //Output file name extraction from input file.
    string qNum = "1";
    string testCaseNum = inputFileName.substr(inputFileName.find('_') + 1, inputFileName.find('.') - inputFileName.find('_' - 1));
    stringstream outputFileName;

    outputFileName << "sol_q" << qNum << "_t" << testCaseNum;

    //Write output to output file.
    ofstream outputFile(outputFileName.str());
    if(!outputFile){
        cerr << "Error: Unable to open output file " << outputFileName.str() << endl;
        return 1;
    }

    for(const string& motif : bestMotifs){
        outputFile << motif << endl;
        //Testing output
        //cout << motif << endl;
    }
    outputFile.close();

    return 0;

}