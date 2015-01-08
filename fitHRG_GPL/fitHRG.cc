// ****************************************************************************************************
// *** COPYRIGHT NOTICE *******************************************************************************
// fitHRG - fits a hierarchical random graph (hrg) model to data
// Copyright (C) 2005-2009 Aaron Clauset
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// See http://www.gnu.org/licenses/gpl.txt for more details.
// 
// ****************************************************************************************************
// Author       : Aaron Clauset  ( aaronc@santafe.edu | http://www.santafe.edu/~aaronc/ )
// Collaborators: Cristopher Moore and Mark E.J. Newman
// Project      : Hierarchical Random Graphs
// Location     : University of New Mexico, Dept. of Computer Science AND Santa Fe Institute
// Created      : 26 October 2005
// Modified     : many, many times
//			   27 December 2007 (cleaned up for public consumption)
// 
// ****************************************************************************************************
// 
// This program fits the HRG model to an input graph by sampling the space of dendrogram models with
// probability proportional to the likelihood that they generate the input graph.
// 
// The program runs the MCMC forever, so when you feel it has converged on the equilibrium set of
// HRG models, you can terminate it at your convenience. The output .hrg file is updated periodically 
// to contain the HRG model with the highest likelihood. The corresponding .info file contains some
// summary information about this run of the program.
//
// ****************************************************************************************************
// 
//  See http://www.santafe.edu/~aaronc/randomgraphs/ for more information, updates to the code, etc.
//
// ****************************************************************************************************
// *** PROGRAM USAGE NOTES ****************************************************************************
// 
// The input to the algorithm must be a text file containing an edge list for the graph in question; 
// nodes are indexed by integers only, indices are separated by a tab, and edges are terminated by a 
// carriage return. Multi-edges and self-loops may appear, but will be stripped out automatically. 
// The graph can be comprised of disconnected components.
//
// For instance, here is a pair of triangles linked by a single edge:
//
// 1 	2
// 1		3
// 2		3
// 4		5
// 4		6
// 5		6
// 1		4
//
// If the input .pairs file is formatted incorrectly, the program will crash.
//
// ****************************************************************************************************

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include "stdlib.h"
#include "time.h"

#include "dendro.h"
#include "graph.h"
#include "rbtree.h"

using namespace std;

// ******** Function Prototypes ***************************************************************************

bool		markovChainMonteCarlo();
string    num2str(const unsigned int);
bool		parseCommandLine(int argc, char * argv[]);
bool		readInputFile();
void		recordHRG(const int, const int, const double);
void		recordNamesLUT();

// ******** Structures and Constants **********************************************************************

struct ioparameters {
	int			n;				// number vertices in input graph
	int			m;				// number of edges in input graph

	string		d_dir;			// working directory
	string		f_in;			// name of input file (either .pairs or .hrg)
	bool			flag_f;			// flag for if -f invoked
	string		f_dg;			// name of output hrg file
	string		f_dg_info;		// name of output information-on-hrg file
	string		f_stat;			// name of output statistics file
	string		f_pairs;			// name of output random graph file
	string		f_namesLUT;		// name of output names LUT file
	bool			flag_make;		// flag for if -make invoked
	string		s_scratch;		// filename sans extension
	string		s_tag;			// user defined filename tag
	string		start_time;		// time simulation was started
	int			timer;			// timer for reading input
	bool			flag_timer;		// flag for when timer fires
	bool			flag_compact;		// compact the Lxy file
};

// ******** Global Variables ******************************************************************************

ioparameters	ioparm;				// program parameters
rbtree		namesLUT;				// look-up table; translates input file node names to graph indices
dendro*		d;					// hrg data structure
unsigned int	t;					// number of time steps max = 2^32 ~ 4,000,000,000
double		bestL;				// best likelihood found so far
int			out_count;			// counts number of maximum found
unsigned int	period  = 10000;		// number of MCMC moves to do before writing stuff out; default: 10000
double*		Likeli;				// holds last k hrg likelihoods

// ******** Main Loop *************************************************************************************

int main(int argc, char * argv[]) {
	ioparm.n		= 0;					// 
	ioparm.timer   = 20;				// 
	ioparm.flag_make = ioparm.flag_f = false;
	ioparm.flag_compact = true;
	ioparm.s_tag   = "";				// 
	string input   = "";				// 
	t			= 1;					// 
	out_count		= 1;					// 
	time_t t1, t2; t1  = time(&t1);		// 
	int maxit=10000;

	if (parseCommandLine(argc, argv)) {
		d = new dendro;						// make the dendro-graph structure
		ioparm.start_time = asctime(localtime(&t1));
		if (ioparm.flag_f) {
			Likeli = new double [period];			// allocate space
			if (!readInputFile()) { cout << "Error: Malformed input file.\n"; return 0; }
			bestL = d->getLikelihood();			// store current likelihood
			//cout << "\nstep   \tLogL       \tbest LogL\tMC step\t\tdLogL\n";
			for (int i=0;i<maxit;i++) {
				if (!(markovChainMonteCarlo())) { return 0; }
				if (t >= 4294967294 or t < 0) { t = 1; }	// rollover step count
			}
		} else if (ioparm.flag_make) {
			if (!(d->importDendrogramStructure(ioparm.f_in))) { cout << "Error: Malformed input file.\n"; return 0; }
			d->makeRandomGraph();
			d->recordGraphStructure(ioparm.f_pairs);
		}
		return 1;
	} else { return 0; } 

}

// ******** Function Definitions **************************************************************************

bool markovChainMonteCarlo() {
	
	double* gptr;
	double  dL;
	bool    flag_taken;

	// Because moves in the dendrogram space are chosen (Monte Carlo) so that we sample dendrograms 
	// with probability proportional to their likelihood, a likelihood-proportional sampling of 
	// the dendrogram models would be equivalent to a uniform sampling of the walk itself. We would
	// still have to decide how often to sample the walk (at most once every n steps is recommended)
	// but for simplicity, the code here simply runs the MCMC itself. To actually compute something
	// over the set of sampled dendrogram models (in a Bayesian model averaging sense), you'll need
	// to code that yourself.
	
	// do 'period' MCMC moves before doing anything else
	for (unsigned int i=0; i<period; i++) {
		
		if (!(d->monteCarloMove(dL, flag_taken))) {	// make a MCMC move
			return false; }

		Likeli[i] = d->getLikelihood();			// get likelihood of this D given G
		if (Likeli[i] > bestL) {
			bestL = Likeli[i];					// store the current best likelihood
			recordHRG(t, out_count, bestL);		// write the hrg structure to file
			out_count++;						// increment count of record-breakings
		}									// 
		
		// Write some stuff to standard-out to describe the current state of things.
		/*		if ((ioparm.flag_compact and ((i+(t-period)) % (8192) == 1)) or !ioparm.flag_compact) {
			cout << "[" << t << "]\t" << Likeli[i];		 
			cout << "   \t(" << bestL << ")\t";			 
			if (flag_taken) { cout << "*\t"; } else { cout << " \t"; }
			cout << "\t" << dL << endl;			
			}
		*/
		if (t > 2147483640) { t = 1; } else { t++; }
	}										// 
	d->refreshLikelihood();						// corrects floating-point errors O(n)
	
	// record the last batch of likelihoods
	if (t-1==period) { ofstream ftest(ioparm.f_stat.c_str(), ios::trunc); ftest.close(); }
	ofstream fout(ioparm.f_stat.c_str(), ios::app);
	for (unsigned int i=0; i<period; i++) {
		if (ioparm.flag_compact and ((i+(t-period)) % (2048) == 1)) {
			fout << i+(t-period) << "\t" << Likeli[i] << "\n";
		} else if (!ioparm.flag_compact) {
			fout << i+(t-period) << "\t" << Likeli[i] << "\n";
		}
	}
	fout.close();
	
	return true;
}

// ********************************************************************************************************

string num2str(const unsigned int input) {
	// input must be a positive integer
	unsigned int temp = input;
	string str  = "";
	if (input == 0) { str = "0"; } else {
		while (temp != 0) {
			str  = char(int(temp % 10)+48) + str;
			temp = (unsigned int)temp/10;
		}
	}
	return str;
}

// ********************************************************************************************************

bool parseCommandLine(int argc, char * argv[]) {
	int argct = 1;
	string temp, ext;
	string::size_type pos;
	bool safeExit = false;

	if (argc==1) {
		cout << "\n  -- Hierarchical Random Graphs --\n";
		cout << "  by Aaron Clauset (copyright 2005-2008)\n\n";
		cout << "  fitHRG is a command line program that takes a simple graph file and runs\n";
		cout << "  a Markov chain Monte Carlo algorithm to sample dendrogram models with\n";
		cout << "  probability proportional to their likelihood of generating the observed\n";
		cout << "  data. Any such model can then be used as a null-model for subsequent\n";
		cout << "  statistical studies of the network's properties.\n\n";
		cout << "  -f <file>       Input .pairs graph file\n";
		cout << "  -make <file>    Build random graph from dendrogram file\n";
		cout << "  -t <string>     (optional) Label for this run\n";
		cout << "\n";
		cout << "  ./fitHRG -f graph.pairs\n";
		cout << "  ./fitHRG -f graph.pairs -t test\n";
		cout << "  ./fitHRG -make test-dendrogram.hrg\n";
		cout << "\n";
		return false;
		
	} else {
		while (argct < argc) {
			temp = argv[argct];
			
			if (temp == "-make" and !ioparm.flag_f) {
				ioparm.flag_make = true;				// -make is mutually exclusive with -f
				argct++;
				temp = argv[argct];
				ext = ".hrg";
				pos = temp.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .hrg format.\n"; return safeExit; }
				ioparm.f_in = temp;
				ext = "/";
				pos = string::npos;
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '/') { pos = i; } }
				if (pos != string::npos) {
					ioparm.d_dir = temp.substr(0, pos+1);
					temp = temp.substr(pos+1,temp.size()-pos-1);
				}
				// now grab the filename sans extension for building outputs files
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '.') { pos = i; } }
				ioparm.s_scratch = temp.substr(0,pos);
				ioparm.f_pairs   = ioparm.d_dir + ioparm.s_scratch + "-random.pairs";
				
				safeExit         = true;
				
			} else if (temp == "-t")       { argct++; ioparm.s_tag = argv[argct];
			} else if (temp == "-compact") { ioparm.flag_compact = true;
			} else if (temp == "-f" and !ioparm.flag_make) {
				ioparm.flag_f = true;				// -f is mutually exclusive with -make
				argct++;
				temp = argv[argct];
				ext = ".pairs";
				pos = temp.find(ext,0);
				if (pos == string::npos) { cout << " Error: Input file must claim to be .pairs format.\n"; return safeExit; }
				ioparm.f_in = temp;
				ext = "/";
				pos = string::npos;
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '/') { pos = i; } }
				if (pos != string::npos) {
					ioparm.d_dir = temp.substr(0, pos+1);
					temp = temp.substr(pos+1,temp.size()-pos-1);
				}
				// now grab the filename sans extension for building outputs files
				for (int i=0; i < temp.size(); i++) { if (temp[i] == '.') { pos = i; } }
				ioparm.s_scratch = temp.substr(0,pos);
//				ioparm.f_stat    = ioparm.d_dir + ioparm.s_scratch + "-L.xy";
				safeExit         = true;
				
			} else { cout << " Warning: ignored argument " << argct << " : " << temp << endl; }
			argct++;
		}
	}
	ioparm.f_namesLUT = ioparm.d_dir + ioparm.s_scratch + "-names.lut";
	if (ioparm.s_tag != "")    { ioparm.s_scratch += "_" + ioparm.s_tag; } 
	if (ioparm.flag_make)      { ioparm.f_pairs    = ioparm.d_dir + ioparm.s_scratch + "-random.pairs"; }
	if (ioparm.flag_f)         { ioparm.f_stat     = ioparm.d_dir + ioparm.s_scratch + "-L.xy";  }

	return safeExit;
}

// ********************************************************************************************************

bool readInputFile() {

	int n,m,s,f,a,b;    n = m = 0;
	int *degs;
	elementrb *item;
	time_t t1; t1 = time(&t1);
	time_t t2; t2 = time(&t2);

	// First, we scan through the input file to create a list of unique node names
	// (which we store in the namesLUT), and a count of the number of edges.
	//	cout << ">> input file scan ( " << ioparm.f_in << " )" << endl;
	//cout << "   edges: [0]"<<endl;
	ifstream fscan1(ioparm.f_in.c_str(), ios::in);
	while (fscan1 >> s >> f) {					// read friendship pair (s,f)
		if (s != f) {
			m++;								// count number of edges
			if (namesLUT.findItem(s) == NULL) { namesLUT.insertItem(s, n++); }
			if (namesLUT.findItem(f) == NULL) { namesLUT.insertItem(f, n++); }
		}
		
		if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
		  //	cout << "   edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;		// 
		}									// 
		t2=time(&t2);							// 
	}
	fscan1.close();
	//cout << "   edges: ["<<m<<"]"<<endl;
	d->g = new graph (n);						// make new graph with n vertices

	// Finally, we reparse the file and added edges to the graph
	m = 0;
	ioparm.flag_timer = false;					// reset timer
	
	//cout << ">> input file read ( " << ioparm.f_in << " )" << endl;
	//cout << "   edges: [0]"<<endl;
	ifstream fin(ioparm.f_in.c_str(), ios::in);
	while (fin >> s >> f) {
		m++;
		if (s != f) {
			item = namesLUT.findItem(s); a = item->value;
			item = namesLUT.findItem(f); b = item->value;
			if (!(d->g->doesLinkExist(a,b))) { if (!(d->g->addLink(a,b))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->g->getName(a) == "") { d->g->setName(a, num2str(s)); } }
			if (!(d->g->doesLinkExist(b,a))) { if (!(d->g->addLink(b,a))) { cout << "Error: (" << s << " " << f << ")" << endl; } else if (d->g->getName(b) == "") { d->g->setName(b, num2str(f)); } }
		}		
		if (t2-t1>ioparm.timer) {				// check timer; if necessarsy, display
		  //	cout << "   edges: ["<<m<<"]"<<endl;
			t1 = t2; ioparm.flag_timer = true;		// 
		}									// 
		t2=time(&t2);							// 
	}
	//	cout << ">> edges: ["<<m<<"]"<<endl;
	fin.close();
	ioparm.m = d->g->numLinks();					// store actual number of directional edges created
	ioparm.n = d->g->numNodes();					// store actual number of nodes used
	//cout << "vertices: ["<<ioparm.n<<"]"<<endl;
	
	recordNamesLUT();							// record names LUT to file for future reference
	//cout << ">> recorded names look-up table" << endl;
	
	d->buildDendrogram();

	return true;
}

// ********************************************************************************************************

void recordHRG(const int step, const int count, const double thisL) {
	
	time_t t1;
	
	// write hrg to file
	ioparm.f_dg = ioparm.d_dir + ioparm.s_scratch + "_best-dendro.hrg";
	d->recordDendrogramStructure(ioparm.f_dg);
	
	// write statistics about hrg to file
	ioparm.f_dg_info = ioparm.d_dir + ioparm.s_scratch + "_best-dendro.info";

	t1 = time(&t1); 
	ofstream fout(ioparm.f_dg_info.c_str(), ios::trunc);
	fout << "---HIERARCHICAL-RANDOM-GRAPH---\n";
	fout << "StartTime     : " << ioparm.start_time;
	fout << "InputFile     : " << ioparm.f_in    << "\n";
	fout << "Directory     : " << ioparm.d_dir   << "\n";
	fout << "---Basic-Graph-Information---\n";
	fout << "Nodes         : " << ioparm.n		<< "\n";
	fout << "Edges         : " << ioparm.m/2	<< "\n";
	fout << "---HRG-Information---\n";
	fout << "OutputTime    : " << asctime(localtime(&t1));
	fout << "NumStepsMCMC  : " << step			<< "\n";
	fout << "NumPrevBests  : " << count-1		<< "\n";
	fout << "LogLikelihood : " << thisL		<< "\n";
	fout << "HRG           : " << ioparm.s_scratch + "_best-dendro.hrg" << "\n";
	fout << "InfoFile      : " << ioparm.s_scratch + "_best-dendro.info"  << "\n";
	fout.close();
	
	return;
}

// ********************************************************************************************************

void recordNamesLUT() {
	rbtree reverseNamesLUT;
	keyValuePair *head, *prev;
	
	head = namesLUT.returnTreeAsList();
	while (head != NULL) {
		reverseNamesLUT.insertItem(head->y, head->x);
		prev = head;
		head = head->next;
		delete prev;
	}
	head = NULL; prev = NULL;
	
	elementrb *item;
	ofstream fout(ioparm.f_namesLUT.c_str(), ios::trunc);
	fout << "virtual\treal\n";
	for (int i=0; i<ioparm.n; i++) {
		item = reverseNamesLUT.findItem(i);
		fout << i << "\t" << item->value << "\n";
	}
	fout.close();
	
	return;
}

// ********************************************************************************************************
// ********************************************************************************************************
