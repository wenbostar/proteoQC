#include <Rcpp.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath> 

//R Stuff

// [[Rcpp::plugins(cpp11)]]


extern "C" 
{
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
}


using namespace Rcpp;

void split(std::string& s, std::string& delim,std::vector< std::string >* ret)  
{  
    size_t last = 0;  
    size_t index=s.find_first_of(delim,last);  
    while (index!=std::string::npos)  
    {  
        ret->push_back(s.substr(last,index-last));  
        last=index+1;  
        index=s.find_first_of(delim,last);  
    }  
    if ((index-last)>0){  
        ret->push_back(s.substr(last,index-last));  
    }  
} 

double str2double(std::string &s) {
	double i;
	std::stringstream ss(s);
	ss >> i;
	return i;
}

RcppExport SEXP ChargeCount_Cpp(SEXP mgf_)
{
	std::string mgf = as<std::string> (mgf_);
	typedef std::map<std::string, int> hash_count;
	hash_count count;


	std::ifstream in(mgf.c_str());
	if(!in)
	{
		return wrap(NA_REAL);
	}
	std::string line;  
	getline(in, line);
    int total=0;
	int ms2flag=0;
	while(in)
	{
        size_t tfound = line.find("BEGIN IONS");
        if(tfound!=std::string::npos)
        {
            total++;
			ms2flag=1;
        }
		size_t endfound = line.find("END IONS");
		if(endfound!=std::string::npos){
			ms2flag=0;
		}
        size_t cfound = line.find("CHARGE=");
        if(cfound!=std::string::npos)
        {
            //std::string charge=line.substr(7,1);
			std::string charge=line.substr(7,std::string::npos);
			if(ms2flag==1){
				count[charge]++;
			}
            //std::cout<<l<<std::endl;
        }
		//std::cout<<line<<std::endl;
		getline(in, line);
	}
	in.close();
    
    Rcpp::List dict(count.size());

    
    int sum=0;
	for(hash_count::iterator it=count.begin();it != count.end();++it)
	{
        dict[it->first]=it->second;
        sum=sum+it->second;
	}
	if(total!=sum){
		dict["0"]=total-sum;
	}
	return(dict);

}


RcppExport SEXP LableRatio_Cpp(SEXP mgf_, SEXP iClass_, SEXP delta_)
{
	std::string mgf = as<std::string> (mgf_);
	double delta = as<double> (delta_);
	int iClass = as<int> (iClass_);
	typedef std::map<std::string, int> hash_count;
	hash_count count;
	hash_count nHash;
	nHash["1"]=0;
	nHash["2"]=0;
	nHash["3"]=0;
	nHash["4"]=0;
	nHash["5"]=0;
	nHash["6"]=0;
	nHash["7"]=0;
	nHash["8"]=0;
	nHash["9"]=0;


	std::ifstream in(mgf.c_str());
	if(!in)
	{
		return wrap(NA_REAL);
		
	}



	std::string line;  
	getline(in, line);
    int total=0;
	//int ms2flag=0;
	std::string sep=" ";

	double I113=113.1;
	double I114=114.1;
	double I115=115.1;
	double I116=116.1;
	double I117=117.1;
	double I118=118.1;
	double I119=119.1;
	double I121=121.1;

	double reporterInt[8]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
	std::vector<std::string> ret;
	while(in)
	{
        //size_t tfound = line.find("BEGIN IONS");
		std::string beginStr = line.substr(0,3);
		std::string firstNum = line.substr(0,1);
        //if(tfound!=std::string::npos)
		if(firstNum=="1"){
			//indicate that find the value
			//match to fragment ions
			ret.clear();
			split(line,sep,&ret);
			double mz = str2double(ret[0]);
			double intensity = str2double(ret[1]);
			if(mz <= 122.0 && iClass==1){
				////iTRAQ-8plex
				if(std::abs(mz-I113)<=delta){
					reporterInt[0]+=intensity;
				}else if(std::abs(mz-I114)<=delta){
					reporterInt[1]+=intensity;
				}else if(std::abs(mz-I115)<=delta){
					reporterInt[2]+=intensity;
				}else if(std::abs(mz-I116)<=delta){
					reporterInt[3]+=intensity;
				}else if(std::abs(mz-I117)<=delta){
					reporterInt[4]+=intensity;
				}else if(std::abs(mz-I118)<=delta){
					reporterInt[5]+=intensity;
				}else if(std::abs(mz-I119)<=delta){
					reporterInt[6]+=intensity;
				}else if(std::abs(mz-I121)<=delta){
					reporterInt[7]+=intensity;
				}
			}

		}else if(beginStr=="BEG"){
            total++;
			//ms2flag=1;
		}else if(beginStr=="END"){
			
		
		//size_t endfound = line.find("END IONS");
		//std::string endStr = line.substr(0,8);

		//if(endfound!=std::string::npos){
		//if(endStr=="END IONS"){
			//ms2flag=0;
			if(iClass==1){
				//iTRAQ-8plex
				std::string reporterMZ = "";
				if(reporterInt[0]>0.0){
					reporterMZ += "I113";
				}

				/*
				reporterMZ+=reporterInt[0]>0.0?"I113":"";
				reporterMZ+=reporterInt[1]>0.0?"I114":"";
				reporterMZ+=reporterInt[2]>0.0?"I115":"";
				reporterMZ+=reporterInt[3]>0.0?"I116":"";
				reporterMZ+=reporterInt[4]>0.0?"I117":"";
				reporterMZ+=reporterInt[5]>0.0?"I118":"";
				reporterMZ+=reporterInt[6]>0.0?"I119":"";
				reporterMZ+=reporterInt[7]>0.0?"I121":"";
				*/
				if(reporterInt[1]>0.0){
					reporterMZ += "I114";
				}
				
				if(reporterInt[2]>0.0){
					reporterMZ += "I115";
				}
				if(reporterInt[3]>0.0){
					reporterMZ += "I116";
				}
				if(reporterInt[4]>0.0){
					reporterMZ += "I117";
				}
				if(reporterInt[5]>0.0){
					reporterMZ += "I118";
				}
				if(reporterInt[6]>0.0){
					reporterMZ += "I119";
				}
				if(reporterInt[7]>0.0){
					reporterMZ += "I121";
				}

				if(reporterMZ!=""){
					count[reporterMZ]++;
				}else{
					count["none"]++;
				}

				//reset to zero
				//reporterInt ={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
				for(int i=0;i<8;i++){
					reporterInt[i]=0.0;
				}	
			}else{
			
			}

		}
		
		getline(in, line);
	}
	in.close();
    
	
    Rcpp::List dict(count.size());
	std::cout<<total<<std::endl;
    
    int sum=0;
	for(hash_count::iterator it=count.begin();it != count.end();++it)
	{
        dict[it->first]=it->second;
        sum=sum+it->second;
	}
	//if(total!=sum){
	//	dict["0"]=total-sum;
	//}
	return(dict);
	

}




/*
 * fixed the RECOMMENDED "registering-native-routines" of bioconductor'check.
*/
extern "C"
{
    void R_init_proteoQC(DllInfo *dll)
	{
		/*
		 * Register routines,allocate resources.
		 * Currently we call all of the functions whith .Call
		 */
		R_CallMethodDef callEntries[] = {
			{"ChargeCount_Cpp", (DL_FUNC) &ChargeCount_Cpp, 1},
			{"LableRatio_Cpp", (DL_FUNC) &LableRatio_Cpp, 3},
			{NULL,NULL,0}
		};
		R_registerRoutines(dll,NULL,callEntries,NULL,NULL);
	}
	void R_upload_proteoQC(DllInfo *dll)
	{
		/* Release resources. */
	}

}