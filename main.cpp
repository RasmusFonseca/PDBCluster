/*
   PDBCluster - k-medoids and k-centers clustering of molecular structures.
   Copyright (C) 2014 Rasmus Fonseca

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Rasmus Fonseca - fonseca.rasmus@gmail.com
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <sys/time.h>

#include "Matrix.h"

using namespace std;

//From http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}

/** Read filenames (assumes they are valid pdb-files) and places atom-distance matrices in ret */
void fillMatrices(vector<string> &files, vector<string> &atomTypes, int firstRes, int lastRes, vector<Matrix*> &ret)
{
	//for( string fName: files ){
	for(unsigned int i=0;i<files.size(); i++){
		string fName = files[i];
		//cout<<fName<<endl;
		int atomCount = 0;
		ifstream fs(fName.c_str(), ifstream::in);
		string token;
		while (fs>>token) {
			if(token=="ATOM"){
				string name, tmp, resiStr;
				fs>>tmp>>name>>tmp>>tmp>>resiStr;
				int resi = atoi(resiStr.c_str());
				if(find(atomTypes.begin(), atomTypes.end(), name) != atomTypes.end() && resi>=firstRes && resi<=lastRes)
					atomCount++;
			}
		}
		fs.close();

		Matrix* coords = new Matrix(atomCount,3);
		atomCount=0;
		ifstream fs2(fName.c_str(), ifstream::in);
		string line;
		while (std::getline(fs2, line)){
			if(line.find("ATOM")==0){
				string name = line.substr(12,5);
				name = trim(name);
				string resiStr = line.substr(22,5);
				resiStr = trim(resiStr);
				int resi = atoi(resiStr.c_str());
				if(find(atomTypes.begin(), atomTypes.end(), name) != atomTypes.end() && resi>=firstRes && resi<=lastRes){
					string x = line.substr(30,8);
					string y = line.substr(38,8);
					string z = line.substr(46,8);
					coords->set(atomCount,0, atof(x.c_str()));
					coords->set(atomCount,1, atof(y.c_str()));
					coords->set(atomCount,2, atof(z.c_str()));
					atomCount++;
				}
			}
		}
		fs2.close();

		Matrix* d = new Matrix(atomCount, atomCount);
		for(int i=0;i<atomCount;i++){
			for(int j=i+1;j<atomCount;j++){
				double dx = coords->get(i,0)-coords->get(j,0);
				double dy = coords->get(i,1)-coords->get(j,1);
				double dz = coords->get(i,2)-coords->get(j,2);
				double dist = sqrt(dx*dx+dy*dy+dz*dz);
				d->set( i,j,dist );
				d->set( j,i,dist );
			}
		}
		ret.push_back(d);
		delete coords;
		//delete cTrans;
	}
}
/*
   ATOM      1 O5'    G A   1      59.712  40.180-111.625  1.00  0.00           O
   ATOM      2 C5'    G A   1      61.014  39.722-111.985  1.00  0.00           C
   ATOM      3 C4'    G A   1      61.399  38.449-111.242  1.00  0.00           C
   */

/** Return true iff vector contains val */
bool contains(vector<int>& vec, int val){
	for(unsigned int i=0;i<vec.size();i++){
		if(vec[i]==val) return true;
	}
	return false;
}

/** Associates the cost of a set of medoids with the average distance from any site to its nearest medoid */
double medoidCosts_distSum(vector<int> &medoids, Matrix &dMat){
	double ret = 0;
	for(unsigned int i=0;i<dMat.cols();i++){
		double min = 1000000000.0;

		for(unsigned int j=0;j<medoids.size();j++){
			int m = medoids[j];
			double d = dMat.get(m,i);
			if( d<min ) min = d;
		}
		ret+=min;
	}
	return ret/dMat.cols();

}

/** Associtates the cost of a set of medoids with the average largest distance between each medoid and the furthest member of its cluster. */
double medoidCosts_maxDistSum(vector<int> &medoids, Matrix &dMat){
	double maxDists[medoids.size()];
	for(unsigned int i=0;i<medoids.size(); i++) maxDists[i] = 0;

	for(unsigned int i=0;i<dMat.cols();i++){
		double min = 1000000000.0;
		int minMedoid = -1;

		for(unsigned int j=0;j<medoids.size();j++){
			int m = medoids[j];
			double d = dMat.get(m,i);
			if( d<min ) {
				min = d;
				minMedoid = j;
			}
		}
		if(min>maxDists[minMedoid]) maxDists[minMedoid] = min;
	}
	double ret = 0;
	for(unsigned int j=0;j<medoids.size();j++)
		ret+=maxDists[j];
	return ret/medoids.size();
}

/** Associates the cost of a set of medoids with the largest distance between any medoid and a member of its cluster. */
double medoidCosts_maxDist(vector<int> &medoids, Matrix &dMat){
	double maxDists[medoids.size()];
	for(unsigned int i=0;i<medoids.size(); i++) maxDists[i] = 0;

	for(unsigned int i=0;i<dMat.cols();i++){
		double min = 1000000000.0;
		int minMedoid = -1;

		for(unsigned int j=0;j<medoids.size();j++){
			int m = medoids[j];
			double d = dMat.get(m,i);
			if( d<min ) {
				min = d;
				minMedoid = j;
			}
		}
		if(min>maxDists[minMedoid]) maxDists[minMedoid] = min;
	}

	double ret = 0;
	for(unsigned int j=0;j<medoids.size();j++)
		if(maxDists[j]>ret) ret=maxDists[j];
	return ret;
}

/** Calls one of the cost-functions based on the clustering-type */
double medoidCosts(vector<int> &medoids, Matrix &dMat, int clusteringType){
	if(clusteringType==0) return medoidCosts_distSum(medoids, dMat);
	if(clusteringType==1) return medoidCosts_maxDistSum(medoids, dMat);
	if(clusteringType==2) return medoidCosts_maxDist(medoids, dMat);
	cerr<<"Unknown clustering type: "<<clusteringType<<endl;
	exit(-1);
}


double get_wall_time(){
	struct timeval time;
	if (gettimeofday(&time,NULL)){
		//  Handle error
		return 0;
	}
	return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){ return (double)clock() / CLOCKS_PER_SEC; }

int main(int argc, char** argv)
{
	//Print usage
	if( argc==1 ){
		cout<<"Usage: "<<endl<<argv[0]<<" -name <atomtype1> [-name <atomtype2> ..] -clusteringType <distSum|maxDistSum|maxDist> -firstRes <int> -lastRes <int> -clusters <integer> [list of pdb-files]"<<endl;
		cout<<"or: "<<endl<<argv[0]<<" -name <atomtype1> [-name <atomtype2> ..] -dumpDistances <filename>"<<endl;
		cout<<"Example of running k-centers for k=10:"<<endl;
		cout<<argv[0]<<" -name CA -name CB -clusteringType maxDist -clusters 10 sample1.pdb sample2.pdb sample3.pdb ..."<<endl;
		cout<<endl;
		cout<<"Standard value for -name is C4' and CA"<<endl;
		cout<<"Standard value for -clusteringType is distSum"<<endl;
		cout<<"Standard value for -clusters is 10"<<endl;
		cout<<"Standard value for -firstRes is -10000000"<<endl;
		cout<<"Standard value for -lastRes is 10000000"<<endl;
		cout<<endl;
		cout<<"\"-clusteringType distSum\" will perform standard k-medoids clustering where the cost associated with a set of medoids is the average distance from each member to its nearest medoid."<<endl;
		cout<<"\"-clusteringType maxDistSum\" will perform a variant of k-centers clustering where the cost associated with a set of medoids is the average distance from each medoid to its furthest associated member."<<endl;
		cout<<"\"-clusteringType maxDist\" will perform k-centers clustering where the cost associated with a set of medoids is the largest distance from any site to the nearest medoid."<<endl;
		cout<<"For all clustering types, the dRMSD metric is used to indicate the distance between two structures. If -dumpDistances is specified no clustering is performed, the structure dRMSD distance matrix is just written to the specified file."<<endl;
		exit(0);
	}

	//Parse arguments
	vector<string> atomTypes;
	vector<string> files;
	string distFileName;
	int clusters = 10;
	int clusteringType = 0;
	int firstRes = -10000000;
	int lastRes = 10000000;
	for(int i=1;i<argc;i++){
		string arg(argv[i]);
		if( arg=="-name" ){
			string nextArg(argv[++i]);
			atomTypes.push_back(nextArg);
		}else if( arg=="-clusteringType" ){
			string nextArg(argv[++i]);
			if(nextArg=="distSum") clusteringType = 0;
			if(nextArg=="maxDistSum") clusteringType = 1;
			if(nextArg=="maxDist") clusteringType = 2;
		}else if( arg=="-dumpDistances" ){
			string nextArg(argv[++i]);
			distFileName = nextArg;
		}else if( arg=="-clusters" ){
			clusters = atoi(argv[++i]);
		}else if( arg=="-firstRes" ){
			firstRes = atoi(argv[++i]);
		}else if( arg=="-lastRes" ){
			lastRes = atoi(argv[++i]);
		}else if( arg[0]=='-' ){
			cerr<<"Unknown argument: "<<arg<<endl;
			exit(-1);
		}else{
			files.push_back(arg);
		}
	}
	if(atomTypes.empty()){
		atomTypes.push_back("C4'");
		atomTypes.push_back("CA");
	}


	double start_cpu = get_cpu_time();
	double start_wall = get_wall_time();

	vector<Matrix*> strucMatrices;
	cout<<"Reading pdb-files"<<endl;
	fillMatrices(files, atomTypes, firstRes, lastRes, strucMatrices);

	cout<<"Filling matrix"<<endl;
	int s = strucMatrices.size();
	Matrix allDists(s,s);
	for(int i=0; i<s; i++){
		for(int j=i+1; j<s; j++){
			double d = strucMatrices[i]->getDifferenceSum(strucMatrices[j]);
			allDists.set(i,j,d);
			allDists.set(j,i,d);
		}
	}

	//If -dumpDistances was specified distFileName will be nonempty and we only print the distances and exit
	if(!distFileName.empty()){
		cout<<"Writing distances to "<<distFileName<<endl;
		//Check if all files are of the format "newpdb_.*.pdb"
		bool newpdbFiles = true;
		for(unsigned int i=0;i<files.size();i++){
			size_t pos1 = files[i].find("newpdb_");
			size_t pos2 = files[i].find(".pdb", pos1);
			if(pos1==string::npos || pos2==string::npos || pos2<pos1){
				newpdbFiles = false;
				break;
			}
		}

		
		ofstream distFile;
		distFile.open(distFileName);
		if(newpdbFiles){
			for(int i=0; i<s; i++){
				size_t pos1 = files[i].find("newpdb_");
				size_t pos2 = files[i].find(".pdb",pos1);
				string id_i = files[i].substr(pos1+7, pos2-pos1-7); 
				for(int j=i+1; j<s; j++){
					pos1 = files[j].find("newpdb_");
					pos2 = files[j].find(".pdb",pos1);
					string id_j = files[j].substr(pos1+7, pos2-pos1-7); 
					distFile<<id_i<<" "<<id_j<<" "<<allDists.get(i,j)<<endl;
				}
			}
		}else{
			for(int i=0; i<s; i++){
				for(int j=i+1; j<s; j++){
					distFile<<files[i]<<" "<<files[j]<<" "<<allDists.get(i,j)<<endl;
				}
			}
		}
		distFile.close();

		//ofstream distFile;
		//distFile.open(distFileName);
		//for(int i=0; i<s; i++){
		//	for(int j=i+1; j<s; j++){
		//		distFile<<i<<" "<<j<<" "<<allDists.get(i,j)<<endl;
		//	}
		//}
		//distFile.close();
		exit(0);
	}

	if(clusters>s){
		cerr<<"Cannot find "<<clusters<<" clusters in "<<s<<" structures"<<endl;
		exit(-1);
	}

	cout<<"Clustering"<<endl;
	vector<int> medoids;
	for(int i=0;i<clusters;i++){
		int val = rand()%s;
		while(contains(medoids, val)) val=rand()%s;
		medoids.push_back(val);
		//medoids.push_back(rand()%s);
		//medoids.push_back(i);
	}

	//Partitioning around medoids: http://en.wikipedia.org/wiki/K-medoids
	double oldCost = -1, cost = medoidCosts(medoids, allDists, clusteringType);
	while(cost!=oldCost){
		oldCost = cost;
		for(int i=0;i<clusters;i++){
			int medoid = medoids[i];
			for(int j=0;j<s;j++){
				if( contains(medoids,j) ) continue;
				medoids[i] = j;
				double newCost = medoidCosts(medoids, allDists, clusteringType );
				if(newCost<cost){
					cost=newCost;
					medoid = j;
					cout<<"Improved cost: "<<newCost<<endl;
				}
			}
			medoids[i] = medoid;
		}
	}
	sort(medoids.begin(), medoids.end());
	double end_cpu = get_cpu_time();
	double end_wall = get_wall_time();
	cout<<"CPU-time : "<<(end_cpu-start_cpu)<<" seconds"<<endl;
	cout<<"Wall-time: "<<(end_wall-start_wall)<<" seconds"<<endl;

	cout<<"Clusters:"<<endl;
	for(unsigned int i=0;i<files.size();i++){
		double min = 1000000000.0;
		int minMed = -1;

		for(unsigned int j=0;j<medoids.size();j++){
			int m = medoids[j];
			double d = allDists.get(m,i);
			if( d<min ) { min = d; minMed = j; }
		}
		cout<<files[i]<<" "<<minMed<<endl;

	}
	cout<<"Medoids:"<<endl;
	for( int m=0;m<clusters;m++)
		cout<<files[medoids[m]]<<endl;

	return 0;
}

