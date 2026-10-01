#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <set>
#include <vector>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <map>
#include <string.h>
#include <algorithm>
#include <assert.h>
#include <ctime>

using namespace std;





// **************************
// **** Global Variables ****
// **************************

string progname;
string s1, s2, s1cf, s2cf, ref_struct, fit_struct;
string surf1 = "n";
string surf1cf = "n";
string s1ss = "n";
string s2ss = "n";
string s2w = "n";
int min_blosum_value = 0;
int num_iter = 3;
int num_iter_cf = 0;
int f_option = 1;
int surf_option = 0;
int ss_option = 0;
int n_option = 1;
int output_option = 1;

float dist_tolerance_init = 1.5;
float dist_tolerance_cf_init = 2.5;
float dist_tolerance = 1.5;
float dist_tolerance_cf = 2.5;
float d0 = 3.0;	

double u[3][3], t[3];
double u_best[3][3], t_best[3];

// ADDED (not upstream): which search loop is running, for writeCheckpoint().  The two clique
// search loops are textually identical apart from the data they work on, so the stage is set at
// the top of each rather than passed down.
const char *cp_stage = "search";
int cp_total = 0;

int ref_elements_in_product_graph[3000];
int fit_elements_in_product_graph[3000];

float shortestDist[3000];
double cost[3000][3000];
double score[3000][3000];

int spalte[3000];
int alignedRes[3000];
int alignedResBest[3000];

map<char, int> blosum62_map;
map<char, int> :: iterator blosum62_map_iter;

int aa_size = 22;
const char *tc_aa[] = { 
	"GLY", "ALA", "SER", "CYS", "VAL",
	"THR", "ILE", "PRO", "MET", "ASP",
	"ASN", "LEU", "LYS", "GLU", "GLN",
	"ARG", "HIS", "PHE", "TYR", "TRP",
	"CYX", "MSE" };

char oc_aa[] = { 
	'G', 'A', 'S', 'C', 'V',
	'T', 'I', 'P', 'M', 'D',
	'N', 'L', 'K', 'E', 'Q',
	'R', 'H', 'F', 'Y', 'W',
	'C', 'M' };

float blosum62[20][20] = { 
	{ 9, -1, -1, -3, 0, -3, -3, -3, -4, -3, -3, -3, -3, -1, -1, -1, -1, -2, -2, -2 },
	{ -1, 4, 1, -1, 1, 0, 1, 0, 0, 0, -1, -1, 0, -1, -2, -2, -2, -2, -2, -3 },
	{ -1, 1, 5, -1, 0, -2, 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, 0, -2, -2, -2 },
	{ -3, -1, -1, 7, -1, -2, -1, -1, -1, -1, -2, -2, -1, -2, -3, -3, -2, -4, -3, -4 },
	{ 0, 1, 0, -1, 4, 0, -1, -2, -1, -1, -2, -1, -1, -1, -1, -1, -2, -2, -2, -3 },
	{ -3, 0, -2, -2, 0, 6, -2, -1, -2, -2, -2, -2, -2, -3, -4, -4, 0, -3, -3, -2 },
	{ -3, 1, 0, -2, -2, 0, 6, 1, 0, 0, -1, 0, 0, -2, -3, -3, -3, -3, -2, -4 },
	{ -3, 0, -1, -1, -2, -1, 1, 6, 2, 0, -1, -2, -1, -3, -3, -4, -3, -3, -3, -4 },
	{ -4, 0, -1, -1, -1, -2, 0, 2, 5, 2, 0, 0, 1, -2, -3, -3, -3, -3, -2, -3 },
	{ -3, 0, -1, -1, -1, -2, 0, 0, 2, 5, 0, 1, 1, 0, -3, -2, -2, -3, -1, -2 },
	{ -3, -1, -2, -2, -2, -2, 1, -1, 0, 0, 8, 0, -1, -2, -3, -3, -2, -1, 2, -2 },
	{ -3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5, 2, -1, -3, -2, -3, -3, -2, -3 },
	{ -3, 0, -1, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5, -1, -3, -2, -3, -3, -2, -3 },
	{ -1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5, 1, 2, -2, 0, -1, -1 },
	{ -1, -2, -1, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4, 2, 1, 0, -1, -3 },
	{ -1, -2, -1, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4, 3, 0, -1, -2 },
	{ -1, -2, 0, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4, -1, -1, -3 },
	{ -2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6, 3, 1 },
	{ -2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7, 2 },
	{ -2, -3, -2, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2,	11 }
	};

// ************************





// *****************
// **** Classes ****
// *****************

class Atom{
	public:
		Atom();
	    ~Atom();
		Atom( Atom & );
		
	public:		
		char res_name;
		int res_Seq;
		char cf[2];
		
		double R[3];   
		double SR[3];
		
		double R_sup[3];   
		double SR_sup[3];
		
		char ss;   // A: alpha, B: beta, C: coil
};

Atom :: Atom()
{}

Atom :: ~Atom()
{}

Atom :: Atom( Atom & )
{}



class Pair{
	public:
		Pair();
	    ~Pair();
		Pair( Pair & );
		
	public:		
		int ref_index;
	    int fit_index;
};

Pair :: Pair()
{}

Pair :: ~Pair()
{}

Pair :: Pair( Pair & )
{}



class Edge{
	public:
		Edge();
	    ~Edge();
		Edge( Edge & );
		
	public:		
		int pair_i;
	    int pair_j;
};

Edge :: Edge()
{}

Edge :: ~Edge()
{}

Edge :: Edge( Edge & )
{}




class Fragment{
	public:
		Fragment();
	    ~Fragment();
		Fragment( Fragment & );
		
	public:		
		int index[3];		
};

Fragment :: Fragment()
{}

Fragment :: ~Fragment()
{}

Fragment :: Fragment( Fragment & )
{}









// *********************************
// ********* Parse Command *********
// *********************************

void help( string msg ){
  cout << msg
  << "\n\n"
  << "options [required]:                                                                                  " << "\n"
  << "    -s1      structure 1 (PDB format)                                                                " << "\n"
  << "    -s1cf    chemical feature file for structure 1                                                   " << "\n"
  << "    -s2      structure 2 (PDB format)                                                                " << "\n"
  << "    -s2cf    chemical feature file for structure 2                                                   " << "\n"
  << "\n"
  << "options [advanced]:                                                                                  " << "\n"
  << "    -b       BLOSUM62 cutoff value for filtering (default=0)                                         " << "\n"    
  << "    -iter    the number of iteration for CA-based maximum clique search (1-3, default=3)             " << "\n"
  << "                < 1: don't use  CA-based maximum clique search                                       " << "\n"  
  << "    -itercf  the number of iteration for chemical feature-based maximum clique search (1-3, default=0)" << "\n"
  << "                < 1: don't use chemical feature-based maximum clique search                          " << "\n"  
  << "    -f       the number of conserved residue(s) in a fragment pair (1-3, default=1)                  " << "\n"
  << "                >3 : don't use fragment superposition                                                " << "\n"
  << "\n"
  << "    -surf    use surface structure for alignment                                                     " << "\n"
  << "             (This option is only used when a whole protein structure is used for s1)                " << "\n"
  << "               0: don't use this option (default)                                                    " << "\n"
  << "               1: use this option                                                                    " << "\n"
  << "    -surf1   surface residues for structure 1 (n: don't use this option, default=n)                  " << "\n"
  << "    -surf1cf chemical feature for the surface residues structure (n: don't use this option, default=n)" << "\n"
  << "\n"
  << "    -s2w     additional structure to be transferred (PDB format) (n: don't use this option, default=n)" << "\n"
  << "\n"
  << "    -ss      check secondary structure conservation                                                  " << "\n"
  << "               0: don't use this option (default)                                                    " << "\n"
  << "               1: use this option                                                                    " << "\n"
  << "    -s1ss    secondary structure file for structure 1 (n: don't use this option, default=n)          " << "\n"
  << "    -s2ss    secondary structure file for structure 2 (n: don't use this option, default=n)          " << "\n"
  << "\n"
  << "    -n       normalization option                                                                    " << "\n"
  << "               1: normalize using structure with smaller number of chemical feature points (N_min) (default)" << "\n"
  << "               2: normalize using structure with larger number of chemical feature points (N_max)    " << "\n"
  << "               3: normalize using the average of N_min and N_max                                     " << "\n"
  << "               4: normalize using s1 structure                                                       " << "\n"
  << "               5: normalize using s2 structure                                                       " << "\n"
  << "\n"
  << "    -o       output option                                                                           " << "\n"
  << "               0: score only                                                                         " << "\n"
  << "               1: score and superposed s2 structure with matrix (default)                            " << "\n"
  << "               2: score with matrix                                                                  " << "\n"
  << "\n"
  << "    -help" << "\n"
  << "\n"
  << "outputs:                                                                                             " << "\n"
  << "     GA-score           : G-LoSA alignment score                                                     " << "\n"
  << "     ali_struct.pdb     : the PDB coordinates of s2 structure aligned onto s1 structure              " << "\n"
  << "     matrix.txt         : translational and rotational matrix to align s2 structure onto s1 structure" << "\n"
  << "     ali_struct_with.pdb: the PDB coordinates of s2w structure transferred by the matrix             " << "\n";

  exit(1);
}


void parseCommandLine( int argc, char *argv[] ){
  string a;
  progname = *argv;

  if( argc == 1 ) help( "" );

  for( ++argv ; --argc ; ++argv ){
    a = *argv;

    if( a == "-s1" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s1 = *++argv;
    } 	
    else if( a == "-s1cf" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s1cf = *++argv;
    }
    else if( a == "-s1ss" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s1ss = *++argv;
    }
	else if( a == "-s2" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s2 = *++argv;
    }
    else if( a == "-s2cf" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s2cf = *++argv;
    }
    else if( a == "-surf" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      surf_option = (int)atof( *++argv );
    }
    else if( a == "-surf1" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      surf1 = *++argv;
    }
    else if( a == "-surf1cf" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      surf1cf = *++argv;
    }
    else if( a == "-s2ss" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s2ss = *++argv;
    }
    else if( a == "-s2w" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      s2w = *++argv;
    }
    else if( a == "-b" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      min_blosum_value = (int)atof( *++argv );
    } 
    else if( a == "-iter" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      num_iter = (int)atof( *++argv );
    }  
    else if( a == "-itercf" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      num_iter_cf = (int)atof( *++argv );
    }  
    else if( a == "-f" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      f_option = (int)atof( *++argv );
    } 
    else if( a == "-ss" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      ss_option = (int)atof( *++argv );
    } 
    else if( a == "-n" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      n_option = (int)atof( *++argv );
    } 
    else if( a == "-o" ){
      if( --argc == 0 ) help( "ERROR: " + a + " requires value." );
      output_option = (int)atof( *++argv );
    }  
	else if (a=="-help") help( "" );
    else help( "ERROR: Bad option: " + a );
  }
}

// ************************






// ************************************************************************
// *********** Class and Functions for Maximum Clique Detection ***********
// ************************************************************************
// ===============================================================
// There are a class and functions to identify maximum clique based on an improved
// branch and bound algorithm
// Reference:
// Janez Konc and Dusanka Janezic. An improved branch and bound algorithm for the 
// maximum clique problem. MATCH Commun. Math. Comput. Chem., 2007, 58, 569-590.
//
// Modified by H. S. Lee
// ==============================================================

class Maxclique {
	const bool* const* e;
	int pk, level;
	const float Tlimit;
  
	class Vertices {
		class Vertex {
  	    	int i, d;
    	
    		public:
      			void set_i( const int ii ) { i = ii; }
      			int get_i() const { return i; }
      			void set_degree( int dd ) { d = dd; }
      			int get_degree() const { return d; }
      	};
    
    	Vertex *v;
    	int sz;
    	static bool desc_degree( const Vertex vi, const Vertex vj) { return (vi.get_degree() > vj.get_degree()); }
  		
  		public:
    		Vertices( int size ) : sz(0) { v = new Vertex[size]; }
    		~Vertices () {}
    		void dispose() { if (v) delete [] v; }
    		void sort() { std::sort( v, v+sz, desc_degree ); }
    		void init_colors();
    		void set_degrees( Maxclique& );
    		int size() const { return sz; }
    		void push( const int ii ) { v[sz++].set_i(ii); };
   			void pop() { sz--; };
    		Vertex& at( const int ii ) const { return v[ii]; };
    		Vertex& end() const { return v[sz - 1]; };
	};
  
	class ColorClass {
    	int *i;
    	int sz;
  
  		public:
    		ColorClass() : sz(0), i(0) {}
    		ColorClass( const int sz ) : sz(sz), i(0) { init(sz); }
    		~ColorClass() { if (i) delete [] i; }
    		void init( const int sz ) { i = new int[sz]; rewind(); }
    		void push( const int ii ) { i[sz++] = ii; };
    		void pop() { sz--; };
   			void rewind() { sz = 0; };
    		int size() const { return sz; }
    		int& at(const int ii) const { return i[ii]; }
    		ColorClass& operator=( const ColorClass& dh ) {
	      		for (int j = 0; j < dh.sz; j++) i[j] = dh.i[j];      		
    	  		sz = dh.sz;      		
      			return *this;
			}
	};
	
	Vertices V;
	ColorClass *C, QMAX, Q;
	
	class StepCount {
    	int i1, i2;
  
  		public:
    		StepCount() : i1(0), i2(0) {}
    		void set_i1( const int ii )  { i1 = ii; }
    		int get_i1() const { return i1; }
   			void set_i2( const int ii )  { i2 = ii; }
   			int get_i2() const { return i2; }
    		void inc_i1()  { i1++; }
	};
  
	StepCount *S;
  	bool connection(const int i, const int j ) const { return e[i][j]; }
  	bool cut1( const int, const ColorClass& );
 	void cut2( const Vertices&, Vertices& );
  	void color_sort( Vertices& );
  	void expand( Vertices );
  	void expand_dyn( Vertices );
  	void _mcq( int*&, int&, bool );
  	void degree_sort( Vertices &R ) { R.set_degrees(*this); R.sort(); }
	
	public:
		Maxclique( const bool* const*, const int, const float=0.025 );
  		int steps() const { return pk; }
  		void mcq( int* &maxclique, int &sz ) { _mcq( maxclique, sz, false ); }
  		void mcqdyn( int* &maxclique, int &sz ) { _mcq(maxclique, sz, true); }
  		~Maxclique() {
    		if( C ) delete [] C;
    		if( S ) delete [] S;
   			V.dispose();
  	};
};



Maxclique::Maxclique ( const bool* const* conn, const int sz, const float tt) : pk(0), level(1), Tlimit(tt), V(sz), Q(sz), QMAX(sz) {
	assert(conn!=0 && sz>0);
  	
  	for (int i=0; i < sz; i++) V.push(i);
  	
  	e = conn;
  	C = new ColorClass[sz + 1];
  	
  	for(int i=0; i < sz + 1; i++) C[i].init(sz + 1);
  	
  	S = new StepCount[sz + 1];
}



void Maxclique::_mcq( int* &maxclique, int &sz, bool dyn ) { 
	V.set_degrees(*this);
  	V.sort();
  	V.init_colors();
  
  	if( dyn ) {
    	for ( int i=0; i < V.size() + 1; i++ ) {
      		S[i].set_i1(0);
      		S[i].set_i2(0);
    	}
   
		expand_dyn(V);
	}
  	else
    	expand(V);
  	
  	maxclique = new int[QMAX.size()]; 
  	
  	for( int i=0; i<QMAX.size(); i++ ) { 
    	maxclique[i] = QMAX.at(i);
  	}
  	
  	sz = QMAX.size();
}



void Maxclique::Vertices::init_colors() { 
	const int max_degree = v[0].get_degree();
  
  	for( int i = 0; i < max_degree; i++ )
   		v[i].set_degree(i + 1);
  
  	for( int i = max_degree; i < sz; i++ )
    	v[i].set_degree(max_degree + 1);
}



void Maxclique::Vertices::set_degrees( Maxclique &m ) { 
  for( int i=0; i < sz; i++ ) {
    int d = 0;
    
    for( int j=0; j < sz; j++ )
    	if( m.connection(v[i].get_i(), v[j].get_i()) ) d++;
    
    	v[i].set_degree(d);
	}
}



bool Maxclique::cut1( const int pi, const ColorClass &A ) {
	for( int i = 0; i < A.size(); i++ )
    	if( connection(pi, A.at(i)) ) return true;
  
  	return false;
}



void Maxclique::cut2( const Vertices &A, Vertices &B ) {
	for( int i = 0; i < A.size() - 1; i++ ) {
    	if( connection(A.end().get_i(), A.at(i).get_i()) ) 
    		B.push(A.at(i).get_i());
	}
}



void Maxclique::color_sort( Vertices &R ) {
	int j = 0;
  	int maxno = 1;
  	int min_k = QMAX.size() - Q.size() + 1;
  	C[1].rewind();
  	C[2].rewind();
 	int k = 1;
 	
 	for( int i=0; i < R.size(); i++ ) {
    	int pi = R.at(i).get_i();
   		k = 1;
    	
    	while( cut1(pi, C[k]) )
      		k++;
    
    	if( k > maxno ) {
      		maxno = k;
      		C[maxno + 1].rewind();
    	}
    
    	C[k].push(pi);
    
    	if( k < min_k ) {
      		R.at(j++).set_i(pi);
    	}
  	}
  
  	if( j > 0 ) R.at(j-1).set_degree(0);
  	if( min_k <= 0 ) min_k = 1;
 	
 	for( k = min_k; k <= maxno; k++ ) {
    	for( int i = 0; i < C[k].size(); i++ ) {
      		R.at(j).set_i(C[k].at(i));
      		R.at(j++).set_degree(k);
    	}
    }
}



void Maxclique::expand( Vertices R ) {
	while( R.size() ) {
    	if( Q.size() + R.end().get_degree() > QMAX.size() ) {
      		Q.push(R.end().get_i());
      		Vertices Rp(R.size());
      		cut2(R, Rp);
      		
      		if( Rp.size() ) {
        		color_sort(Rp);
				pk++;
        		expand(Rp);
      		}
      		else if( Q.size() > QMAX.size() ) { 
				QMAX = Q;
      		}    
      		
      		Rp.dispose();
      		Q.pop();
    	}
    	else {
      		return;
    	}
    	
    	R.pop();
  	}
}



void Maxclique::expand_dyn( Vertices R ) {
	S[level].set_i1(S[level].get_i1() + S[level - 1].get_i1() - S[level].get_i2());
  	S[level].set_i2(S[level - 1].get_i1());
  
  	while( R.size() ) {
    	if( Q.size() + R.end().get_degree() > QMAX.size() ) {
      		Q.push(R.end().get_i());
     		Vertices Rp(R.size());
      		cut2(R, Rp);
      
      		if( Rp.size() ) {
        		if( (float)S[level].get_i1()/++pk < Tlimit ) {
          			degree_sort(Rp);
        		}
        		
        		color_sort(Rp);
				S[level].inc_i1();
				level++;
				expand_dyn(Rp);
				level--;
      		}
      		else if( Q.size() > QMAX.size() ) { 
				QMAX = Q;
      		}    
     
     		Rp.dispose();
      		Q.pop();
    	}
    	else{
      		return;
    	}
    
    	R.pop();
  	}
}
// *************************************************




char getResidueName( const char* res ){
	for( int i = 0 ; i < aa_size ; i++ ){
		if( strcmp( tc_aa[i], res ) == 0 ){
			return oc_aa[i];
		}
	}

	return 'X';
}




void setBLOSUN62Map(){
	//C  S	 T	 P 	 A	 G	 N	 D	 E 	 Q	 H	 R	 K	 M	 I	 L	 V	 F	 Y	 W
	blosum62_map.insert( pair< char, int > ( 'C', 0 ) );
	blosum62_map.insert( pair< char, int > ( 'S', 1 ) );
	blosum62_map.insert( pair< char, int > ( 'T', 2 ) );
	blosum62_map.insert( pair< char, int > ( 'P', 3 ) );
	blosum62_map.insert( pair< char, int > ( 'A', 4 ) );
	blosum62_map.insert( pair< char, int > ( 'G', 5 ) );
	blosum62_map.insert( pair< char, int > ( 'N', 6 ) );
	blosum62_map.insert( pair< char, int > ( 'D', 7 ) );
	blosum62_map.insert( pair< char, int > ( 'E', 8 ) );
	blosum62_map.insert( pair< char, int > ( 'Q', 9 ) );
	blosum62_map.insert( pair< char, int > ( 'H', 10 ) );
	blosum62_map.insert( pair< char, int > ( 'R', 11 ) );
	blosum62_map.insert( pair< char, int > ( 'K', 12 ) );
	blosum62_map.insert( pair< char, int > ( 'M', 13 ) );
	blosum62_map.insert( pair< char, int > ( 'I', 14 ) );
	blosum62_map.insert( pair< char, int > ( 'L', 15 ) );
	blosum62_map.insert( pair< char, int > ( 'V', 16 ) );
	blosum62_map.insert( pair< char, int > ( 'F', 17 ) );
	blosum62_map.insert( pair< char, int > ( 'Y', 18 ) );
	blosum62_map.insert( pair< char, int > ( 'W', 19 ) );
}



int getValueBLOSUN62Map( char res ){	
	blosum62_map_iter = blosum62_map.find( res );

	if( blosum62_map_iter == blosum62_map.end() )
		return -1;
	
    return blosum62_map_iter->second;
}




int getNumberOfResidues( string filename ){
	char buffer[100];
	int num_residue = 0;

	FILE *pInputFile = fopen( filename.c_str(), "r" );

	while( fgets( buffer, 100, pInputFile ) != NULL ){
		if( buffer[0] == 'T' || buffer[0] == 'E' )
			break;

		if( buffer[0] == 'A' && buffer[1] == 'T' && buffer[13] == 'C' && buffer[14] == 'A' ){
			num_residue++;			
		}
	}

	fclose( pInputFile );
	
	return num_residue;
}





void readCAAtomWithSideChainCentroidFromPDBFile( string filename, vector< Atom * > &vStruct ){
	char buffer[100];

	Atom *atom; 
	char res[4];
	char resSeq[5];
	char X_char[9], Y_char[9], Z_char[9];

	FILE *pInputFile = fopen( filename.c_str(), "r" );

	int num_CA_atoms = 0;
	int num_sc_atoms = 0;
	double sum_sc_X = 0;
	double sum_sc_Y = 0;
	double sum_sc_Z = 0;

	while( fgets( buffer, 100, pInputFile ) != NULL ){		
		if( buffer[0] == 'T' || buffer[0] == 'E' ){
			if( num_sc_atoms == 0 ){
				atom->SR[0] = atom->R[0];
				atom->SR[1] = atom->R[1];
				atom->SR[2] = atom->R[2];	
			}			
			else{
				atom->SR[0] = sum_sc_X/num_sc_atoms;
				atom->SR[1] = sum_sc_Y/num_sc_atoms;
				atom->SR[2] = sum_sc_Z/num_sc_atoms;							
			}	

			vStruct.push_back( atom );	

			break;
		}

		if( buffer[0] == 'A' && buffer[1] == 'T' && buffer[13] == 'C' && buffer[14] == 'A' ){
			num_CA_atoms++;

			if( num_CA_atoms > 1 ){
				if( num_sc_atoms == 0 ){
					atom->SR[0] = atom->R[0];
					atom->SR[1] = atom->R[1];
					atom->SR[2] = atom->R[2];
				}				
				else{
					atom->SR[0] = sum_sc_X/num_sc_atoms;
					atom->SR[1] = sum_sc_Y/num_sc_atoms;
					atom->SR[2] = sum_sc_Z/num_sc_atoms;										
				}

				vStruct.push_back( atom );	

				num_sc_atoms = 0;
				sum_sc_X = 0;
				sum_sc_Y = 0;
				sum_sc_Z = 0;	
			}

			// Residue name
			for( int i = 0 ; i < 3 ; i++ ){
				res[i] = buffer[17+i];
			}

			res[3] = '\0';		
			
			
			// Residue sequence
			for( int i = 0 ; i < 4 ; i++ )
				resSeq[i] = buffer[22+i]; 
		
			resSeq[4] = '\0';		


			// Coordinates
			for( int i = 0 ; i < 8 ; i++ ){
				X_char[i] = buffer[30+i];
				Y_char[i] = buffer[38+i];
				Z_char[i] = buffer[46+i];
			}

			X_char[8] = '\0';
			Y_char[8] = '\0';
			Z_char[8] = '\0';

			atom = new Atom();
			atom->res_name = getResidueName( res );
			atom->res_Seq = (int)atof( resSeq );
			atom->R[0] = (double)atof( X_char );
			atom->R[1] = (double)atof( Y_char );
			atom->R[2] = (double)atof( Z_char );			
		}


		if( buffer[0] == 'A' && buffer[1] == 'T' && !( buffer[13] == 'N' && buffer[14] == ' ' ) && !( buffer[13] == 'C' && buffer[14] == 'A' ) 
			&& !( buffer[13] == 'C' && buffer[14] == ' ' ) && !( buffer[13] == 'O' && buffer[14] == ' ' ) ){
			num_sc_atoms++;

			// Coordinates
			for( int i = 0 ; i < 8 ; i++ ){
				X_char[i] = buffer[30+i];
				Y_char[i] = buffer[38+i];
				Z_char[i] = buffer[46+i];
			}

			X_char[8] = '\0';
			Y_char[8] = '\0';
			Z_char[8] = '\0';

			sum_sc_X += (double)atof( X_char );
			sum_sc_Y += (double)atof( Y_char );
			sum_sc_Z += (double)atof( Z_char );	
		}
	}

	fclose( pInputFile );
}





void readChemicalFeatures( string filename, vector< Atom * > &vStruct ){
	char buffer[100];

	Atom *atom; 
	char res[4];
	char resSeq[5];
	char cf[2];
	char X_char[9], Y_char[9], Z_char[9];

	FILE *pInputFile = fopen( filename.c_str(), "r" );

	while( fgets( buffer, 100, pInputFile ) != NULL ){
		if( buffer[0] == 'T' || buffer[0] == 'E' )
			break;

		if( buffer[0] == 'A' && buffer[1] == 'T' ){
			// Residue name
			for( int i = 0 ; i < 3 ; i++ ){
				res[i] = buffer[17+i];
			}

			res[3] = '\0';
			
			
			// Residue sequence
			for( int i = 0 ; i < 4 ; i++ )
				resSeq[i] = buffer[22+i]; 
		
			resSeq[4] = '\0';	
			
			
			// Chemical feature
			cf[0] = buffer[13]; 
			cf[1] = buffer[14]; 


			// Coordinates
			for( int i = 0 ; i < 8 ; i++ ){
				X_char[i] = buffer[30+i];
				Y_char[i] = buffer[38+i];
				Z_char[i] = buffer[46+i];
			}

			X_char[8] = '\0';
			Y_char[8] = '\0';
			Z_char[8] = '\0';

			atom = new Atom();
			atom->res_name = getResidueName( res );
			atom->res_Seq = (int)atof( resSeq );
			atom->cf[0] = cf[0];
			atom->cf[1] = cf[1];
			atom->R[0] = (double)atof( X_char );
			atom->R[1] = (double)atof( Y_char );
			atom->R[2] = (double)atof( Z_char );
			vStruct.push_back( atom );
		}
	}

	fclose( pInputFile );
}





void setSecondaryStructure( string filename, vector< Atom * > &vStruct ){
	char buffer[100];
	int i = -1;

	FILE *pInputFile = fopen( filename.c_str(), "r" );

	while( fgets( buffer, 100, pInputFile ) != NULL ){
		i++;
		
		if( buffer[0] == 'T' || buffer[0] == 'E' )
			break;

		vStruct[i]->ss = buffer[13];
	}

	fclose( pInputFile );
}





void readCAAtomsFromPDBFile( string filename, vector< Atom * > &vStruct ){
	char buffer[100];

	Atom *atom; 
	char res[4];
	char resSeq[5];
	char X_char[9], Y_char[9], Z_char[9];

	FILE *pInputFile = fopen( filename.c_str(), "r" );

	while( fgets( buffer, 100, pInputFile ) != NULL ){
		if( buffer[0] == 'T' || buffer[0] == 'E' )
			break;

		if( buffer[0] == 'A' && buffer[1] == 'T' && buffer[13] == 'C' && buffer[14] == 'A' ){
			// Residue name
			for( int i = 0 ; i < 3 ; i++ ){
				res[i] = buffer[17+i];
			}

			res[3] = '\0';
			
			
			// Residue sequence
			for( int i = 0 ; i < 4 ; i++ )
				resSeq[i] = buffer[22+i]; 
		
			resSeq[4] = '\0';	


			// Coordinates
			for( int i = 0 ; i < 8 ; i++ ){
				X_char[i] = buffer[30+i];
				Y_char[i] = buffer[38+i];
				Z_char[i] = buffer[46+i];
			}

			X_char[8] = '\0';
			Y_char[8] = '\0';
			Z_char[8] = '\0';

			atom = new Atom();
			atom->res_name = getResidueName( res );
			atom->res_Seq = (int)atof( resSeq );
			atom->R[0] = (double)atof( X_char );
			atom->R[1] = (double)atof( Y_char );
			atom->R[2] = (double)atof( Z_char );
			vStruct.push_back( atom );
		}
	}

	fclose( pInputFile );
}





void GeneratePairs( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, vector< Pair * > &vPair ){
	int num_element_ref = vRefStruct.size();
	int num_element_fit = vFitStruct.size();

	Pair *p; 
	
	if( ss_option == 0 ){
		for( int i = 0 ; i < num_element_ref ; i++ ){
			for( int j = 0 ; j < num_element_fit ; j++ ){
				if( blosum62[getValueBLOSUN62Map( vRefStruct[i]->res_name )][getValueBLOSUN62Map( vFitStruct[j]->res_name )] >= min_blosum_value ){
					p = new Pair();
			    	p->ref_index = i;
					p->fit_index = j;

					vPair.push_back( p );
				}
			}
		}
	}
	else{
		for( int i = 0 ; i < num_element_ref ; i++ ){
			for( int j = 0 ; j < num_element_fit ; j++ ){
				if( ( blosum62[getValueBLOSUN62Map( vRefStruct[i]->res_name )][getValueBLOSUN62Map( vFitStruct[j]->res_name )] >= min_blosum_value ) && ( vRefStruct[i]->ss == vFitStruct[j]->ss ) ){
					p = new Pair();
			    	p->ref_index = i;
					p->fit_index = j;

					vPair.push_back( p );
				}
			}
		}
	}
}




void GeneratePairsForCF( vector< Atom * > &vRefStructCF, vector< Atom * > &vFitStructCF, vector< Pair * > &vPair ){
	int num_element_ref = vRefStructCF.size();
	int num_element_fit = vFitStructCF.size();

	Pair *p; 
	
	for( int i = 0 ; i < num_element_ref ; i++ ){
		for( int j = 0 ; j < num_element_fit ; j++ ){
			if( ( vRefStructCF[i]->cf[0] == vFitStructCF[j]->cf[0] ) && ( vRefStructCF[i]->cf[1] == vFitStructCF[j]->cf[1] ) ){
				p = new Pair();
			    p->ref_index = i;
				p->fit_index = j;

				vPair.push_back( p );
			}
		}
	}
}




void GenerateEdges( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, vector< Pair * > &vPair, vector< Edge * > &vEdge, float dist_tol_cutoff ){
	Edge *e; 
	int num_pairs = vPair.size();
	double dist_1, dist_2, d_dist;	

	for( int i = 0 ; i < num_pairs ; i++ ){
		for( int j = i+1 ; j < num_pairs ; j++ ){
			dist_1 = sqrt( pow( vRefStruct[vPair[i]->ref_index]->R[0] - vRefStruct[vPair[j]->ref_index]->R[0], 2 ) 
				+ pow( vRefStruct[vPair[i]->ref_index]->R[1] - vRefStruct[vPair[j]->ref_index]->R[1], 2 ) 
				+ pow( vRefStruct[vPair[i]->ref_index]->R[2] - vRefStruct[vPair[j]->ref_index]->R[2], 2 ) );

			dist_2 = sqrt( pow( vFitStruct[vPair[i]->fit_index]->R[0] - vFitStruct[vPair[j]->fit_index]->R[0], 2 ) 
				+ pow( vFitStruct[vPair[i]->fit_index]->R[1] - vFitStruct[vPair[j]->fit_index]->R[1], 2 ) 
				+ pow( vFitStruct[vPair[i]->fit_index]->R[2] - vFitStruct[vPair[j]->fit_index]->R[2], 2 ) );

			d_dist = dist_1 - dist_2;

			if( d_dist < dist_tol_cutoff && d_dist > -dist_tol_cutoff ){     
				e = new Edge();
				e->pair_i = i;
				e->pair_j = j;

				vEdge.push_back( e );
			}
		}
	}
}




void resetArrayForElementsInProductGraph(){
	for( int i = 0 ; i < 3000 ; i++ ){
		ref_elements_in_product_graph[i] = 0;
		fit_elements_in_product_graph[i] = 0;
	}
}




void setElementsInProductGraph( vector< Pair * > &vPair, vector< Edge * > &vEdge ){
	int num_edges = vEdge.size();
	int pair_i, pair_j;
	int ref_index_i, fit_index_i;
	int ref_index_j, fit_index_j;
	
	for( int i = 0 ; i < num_edges ; i++ ){
		pair_i = vEdge[i]->pair_i;
		pair_j = vEdge[i]->pair_j;
		
		ref_index_i = vPair[pair_i]->ref_index;
		fit_index_i = vPair[pair_i]->fit_index;
		
		ref_index_j = vPair[pair_j]->ref_index;
		fit_index_j = vPair[pair_j]->fit_index;
		
		if( ( ( ref_index_i > ref_index_j ) && ( fit_index_i > fit_index_j ) ) || ( ( ref_index_i < ref_index_j ) && ( fit_index_i < fit_index_j ) ) )
		
		ref_elements_in_product_graph[ref_index_i] = 1;
		ref_elements_in_product_graph[ref_index_j] = 1;
		
		fit_elements_in_product_graph[fit_index_i] = 1;
		fit_elements_in_product_graph[fit_index_j] = 1;
	}
}




void writePairs( vector< Pair * > &vPair ){
	int num_pairs = vPair.size();
	
	FILE *pOutputFile_pair = fopen( "pairs.rst", "w" );

	fprintf( pOutputFile_pair, "%s", "pair_index\tref_index\tfit_index\n" );

	for( int i = 0 ; i < num_pairs ; i++ ){
		fprintf( pOutputFile_pair, "%d", i );
		fprintf( pOutputFile_pair, "%s", "\t" );
		fprintf( pOutputFile_pair, "%d", vPair[i]->ref_index );
		fprintf( pOutputFile_pair, "%s", "\t" );
		fprintf( pOutputFile_pair, "%d", vPair[i]->fit_index );
		fprintf( pOutputFile_pair, "\n" );
	}

	fprintf( pOutputFile_pair, "END\n" );

	fclose( pOutputFile_pair );
}
	
	
	
	
void writeProductGraph( int num_pairs, vector< Edge * > &vEdge ){
	FILE *pOutputFile = fopen( "product_graph.rst", "w" );	

	fprintf( pOutputFile, "%s", "p edge " );
	fprintf( pOutputFile, "%d", num_pairs );
	fprintf( pOutputFile, "%s", " P=0.5 SEED=12415\n" );
	
	int num_edges = vEdge.size();

	for( int i = 0 ; i < num_edges ; i++ ){
		fprintf( pOutputFile, "%s", "e " );

		if( vEdge[i]->pair_i < vEdge[i]->pair_j ){
			fprintf( pOutputFile, "%d", vEdge[i]->pair_i );
			fprintf( pOutputFile, "%s", " " );
			fprintf( pOutputFile, "%d", vEdge[i]->pair_j );
			fprintf( pOutputFile, "\n" );
		}
		else{
			fprintf( pOutputFile, "%d", vEdge[i]->pair_j );
			fprintf( pOutputFile, "%s", " " );
			fprintf( pOutputFile, "%d", vEdge[i]->pair_i );
			fprintf( pOutputFile, "\n" );
		}
	}

	fclose( pOutputFile );
}
	

	
	
int read_dimacs( string name, bool** &conn, int &size ){
	ifstream f (name.c_str());
    string buffer;
    assert( f.is_open() );
    set<int> v;
    multimap<int,int> e;
    int num_line = 0;
    
    while( !getline( f, buffer ).eof() ){
    	if( buffer[0] == 'e' ){
    		num_line++;
      		int vi, vj;
      		sscanf( buffer.c_str(), "%*c %d %d", &vi, &vj );
      		v.insert(vi);
      		v.insert(vj);
      		e.insert( make_pair( vi, vj ) );
    	}
  	}
  	
  	if( num_line == 0 )
  		return 0;

  	size = *v.rbegin() + 1; 
  	conn = new bool*[size];
 
  	for( int i = 0 ; i < size ; i++ ){
    	conn[i] = new bool[size];
    	memset(conn[i], 0, size * sizeof(bool));
  	}
 
  	for( multimap<int,int>::iterator it = e.begin() ; it != e.end() ; it++ ){
    	conn[it->first][it->second] = true;
    	conn[it->second][it->first] = true;
  	}
  
  	//cout << "|E| = " << e.size() << "  |V| = " << v.size() << " p = " << (double) e.size() / (v.size() * (v.size() - 1) / 2) << endl;
  
  	f.close();
  	
  	return 1;
}




int identifyMaxClique( vector< int > &vPairs ){
	bool **conn;
    int size;		
    
  	int read_dimacs_ok = read_dimacs( "product_graph.rst", conn, size );
  	
  	if( read_dimacs_ok == 0 )
  		return 0;
  	
  	clock_t start1, start2;
  	int *qmax;
  	int qsize;

	start1 = time(NULL);
    start2 = clock();
    
    Maxclique md( conn, size, 0.025 );
    md.mcqdyn( qmax, qsize ); 

       	    	
    // =================================	    	
    // *** set aligned residue pairs ***
    // =================================
       
    for( int i = 0; i < qsize; i++ ){
    	vPairs.push_back( qmax[i] );
    }
    	
    delete [] qmax;
    
    for (int i = 0 ; i < size ; i++ )
      delete [] conn[i];
      
    delete [] conn;
    
        	
    if( qsize < 3 ){
    	return 0;
    }
    
    
    return 1;
}



	
void superposeFitStructure( string filename, string output_name ){
	char buffer[100];

	double r[3];
	double r_prime[3];

	char X_char[9], Y_char[9], Z_char[9];

	FILE *pInputFile = fopen( filename.c_str(), "r" );
	FILE *pOutputFile = fopen( output_name.c_str(), "w" );

	while( 1 ){
		fgets( buffer, 100, pInputFile );

		if( buffer[0] == 'T' || buffer[0] == 'E' )
			break;

		if( ( buffer[0] == 'A' && buffer[1] == 'T' ) || ( buffer[0] == 'H' && buffer[1] == 'E' ) ){
			for( int i = 0 ; i < 30 ; i++ )
				fprintf( pOutputFile, "%c", buffer[i] );

			for( int i = 0 ; i < 8 ; i++ ){
				X_char[i] = buffer[30+i];
				Y_char[i] = buffer[38+i];
				Z_char[i] = buffer[46+i];
			}

			X_char[8] = '\0';
			Y_char[8] = '\0';
			Z_char[8] = '\0';

			r[0] = (double)atof( X_char );
			r[1] = (double)atof( Y_char );
			r[2] = (double)atof( Z_char );

			r_prime[0] = t[0] + u[0][0]*r[0] + u[0][1]*r[1] + u[0][2]*r[2];
			r_prime[1] = t[1] + u[1][0]*r[0] + u[1][1]*r[1] + u[1][2]*r[2];
			r_prime[2] = t[2] + u[2][0]*r[0] + u[2][1]*r[1] + u[2][2]*r[2];		
			
			fprintf( pOutputFile, "%8.3f", r_prime[0] );
			fprintf( pOutputFile, "%8.3f", r_prime[1] );
			fprintf( pOutputFile, "%8.3f", r_prime[2] );
			
			for( int i = 54 ; i < 78 ; i++ )
				fprintf( pOutputFile, "%c", buffer[i] );

			fprintf( pOutputFile, "\n" );
		}		
	}

	fprintf( pOutputFile, "TER" );

	fclose( pOutputFile );
	fclose( pInputFile );	
}





void setSuperposedCoordinatesForFitStructure( vector< Atom * > &vFitStruct ){
	double r[3];
	double r_prime[3];
	int num_residues = vFitStruct.size();
	
	for( int i = 0 ; i < num_residues ; i++ ){		
		// CA atom
		r[0] = vFitStruct[i]->R[0];
		r[1] = vFitStruct[i]->R[1];
		r[2] = vFitStruct[i]->R[2];
			
		r_prime[0] = t[0] + u[0][0]*r[0] + u[0][1]*r[1] + u[0][2]*r[2];
		r_prime[1] = t[1] + u[1][0]*r[0] + u[1][1]*r[1] + u[1][2]*r[2];
		r_prime[2] = t[2] + u[2][0]*r[0] + u[2][1]*r[1] + u[2][2]*r[2];	
			
		vFitStruct[i]->R_sup[0] = r_prime[0];
		vFitStruct[i]->R_sup[1] = r_prime[1];
		vFitStruct[i]->R_sup[2] = r_prime[2];	
		
		// side chain centroid
		r[0] = vFitStruct[i]->SR[0];
		r[1] = vFitStruct[i]->SR[1];
		r[2] = vFitStruct[i]->SR[2];
			
		r_prime[0] = t[0] + u[0][0]*r[0] + u[0][1]*r[1] + u[0][2]*r[2];
		r_prime[1] = t[1] + u[1][0]*r[0] + u[1][1]*r[1] + u[1][2]*r[2];
		r_prime[2] = t[2] + u[2][0]*r[0] + u[2][1]*r[1] + u[2][2]*r[2];	
			
		vFitStruct[i]->SR_sup[0] = r_prime[0];
		vFitStruct[i]->SR_sup[1] = r_prime[1];
		vFitStruct[i]->SR_sup[2] = r_prime[2];			
	}
}





void setRotationMatrix( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	int n = vRefStruct.size();    // The sizes of vRefStruct and vFitStruct must be identical
	int mode = 1;
	int ier;

	int gj, gk, gl, gm1, gm;
	int ip[9] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
	int ip1201[4] = { 1, 2, 0, 1 };

	double r[3][3], xc[3], yc[3], wc;
	double a[3][3], b[3][3], e[3], rr[6], ss[6];
	double e0, d, spur, det, cof, h, g;
	double cth, sth, sqrth, p, sigma;
	double sqrt3, tol;

	sqrt3 = 1.73205080756888;
    tol = 0.01;

	wc = 0.0;
    double rms = 0.0;
    e0 = 0.0;

    for( int i = 0 ; i < 3 ; i++ ){
		xc[i] = 0.0;
        yc[i] = 0.0;
        t[i] = 0.0;

		 for( int j = 0 ; j < 3 ; j++ ){
			 r[i][j] = 0.0;
             u[i][j] = 0.0;
             a[i][j] = 0.0;
			 
			 if( i == j ){
               u[i][j] = 1.0;
               a[i][j] = 1.0;
			 }
		 }
	}

	ier = -1;

	if( n <  1 ) 
		return;
    
	ier = -2;

	for( int m = 0 ; m < n ; m++ ){
		wc = wc + 1.0;

		for( int i = 0 ; i < 3 ; i++ ){
			xc[i] = xc[i] + vFitStruct[m]->R[i];
            yc[i] = yc[i] + vRefStruct[m]->R[i];
		}
	}

	if( wc <= 0 ) return;

	for( int i = 0 ; i < 3 ; i++ ){
		xc[i] = xc[i] / wc;
        yc[i] = yc[i] / wc;
	}

	for( int m = 0 ; m < n ; m++ ){
		for( int i = 0 ; i < 3 ; i++ ){
			e0 = e0 + ( pow( vFitStruct[m]->R[i] - xc[i], 2 ) + pow( vRefStruct[m]->R[i] - yc[i], 2 ) );
            d = vRefStruct[m]->R[i] - yc[i];

			for( int j = 0 ; j < 3 ; j++ ){
				r[i][j] = r[i][j] + d * ( vFitStruct[m]->R[j] - xc[j] );
			}
		}
	}

	det = r[0][0] * ( r[1][1] * r[2][2] - r[1][2] * r[2][1] ) 
		- r[0][1] * ( r[1][0] * r[2][2] - r[1][2] * r[2][0] )
		+ r[0][2] * ( r[1][0] * r[2][1] - r[1][1] * r[2][0] );

	sigma = det;

	int index = -1;

	for( int j = 0 ; j < 3 ; j++ ){
		for( int i = 0 ; i <= j ; i++ ){
			index++;
            rr[index] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
		}
	}

	 spur = ( rr[0] + rr[2] + rr[5] ) / 3.0;

	 cof = ( rr[2]*rr[5] - rr[4]*rr[4] + rr[0]*rr[5] - rr[3]*rr[3] + rr[0]*rr[2] - rr[1]*rr[1] ) / 3.0;

	 det = det * det;

	 for( int i = 0 ; i < 3 ; i++ ){
		 e[i] = spur;
	 }	 	 

	 if( spur <= 0 ){
		 for( int i = 0 ; i < 3 ; i++ ){
			 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
		 }

		 for( int i = 0 ; i < 3 ; i++ ){
			 if( e[i] < 0 ) 
				 e[i] = 0.0;

			 e[i] = sqrt( e[i] );
		 }

		 ier = 0;

		 if( e[1] <= ( e[0] * 0.00005 ) ) 
			 ier = -1;

		 d = e[2];

		 if( sigma < 0.0 ){
			 d = -d;
			  
			 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
				 ier = -1;
		 }

		 d = ( d + e[1] ) + e[0];

		 rms = ( e0 - d ) - d;

		 if( rms < 0.0 ) 
			 rms = 0.0;

		 return;
	 }

	 d = spur * spur;
     h = d - cof;
     g = ( spur * cof - det ) / 2.0 - spur * h;

	 if( h <= 0 ){ 
		 if( mode == 0 ){
			 for( int i = 0 ; i < 3 ; i++ ){
				 if( e[i] < 0 ) 
					 e[i] = 0.0;

				 e[i] = sqrt( e[i] );
			 }

			 ier = 0;

             if( e[1] <= ( e[0] * 0.00005 ) ) 
				 ier = -1;

			 d = e[2];

			 if( sigma < 0.0 ){
				 d = -d;
				  
				 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
					 ier = -1;
			 }

			 d = ( d + e[1] ) + e[0];

			 rms = ( e0 - d ) - d;

			 if( rms < 0.0 ) 
				 rms = 0.0;

			 return;
		 }
		 else{
			 for( int l = 0 ; l < 2 ; l++ ){
				 d = 0.0;

				 for( int i = 0 ; i < 3 ; i++ ){
					 b[i][l] = r[i][0] * a[0][l] + r[i][1] * a[1][l] + r[i][2] * a[2][l];
					 d = d + pow( b[i][l], 2 );
				 }

				 if( d > 0 ) 
					 d = 1.0 / sqrt( d );

				 for( int i = 0 ; i < 3 ; i++ ){
					 b[i][l] = b[i][l] * d;
				 }
			 }

             d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
             p = 0.0;

			 for( int i = 0 ; i < 3 ; i++ ){
				 b[i][1] = b[i][1] - d * b[i][0];
			     p = p + pow( b[i][1], 2 );
			 }

			 if( p <= tol ){
				  p = 1.0;

				  for( int i = 0 ; i < 3 ; i++ ){
					  if( p < abs( b[i][0] ) ){}
					  else{
						  p = abs( b[i][0] );
						  gj = i;
					  }
				  }

				  gk = ip1201[gj];
				  gl = ip1201[gj+1];
				  p = sqrt( pow( b[gk][0], 2 ) + pow( b[gl][0], 2 ) );

				  if( p <= tol ){
					 for( int i = 0 ; i < 3 ; i++ ){
						 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
					 }

					 for( int i = 0 ; i < 3 ; i++ ){
						 if( e[i] < 0 ) 
							 e[i] = 0.0;

						 e[i] = sqrt( e[i] );
					 }

					 ier = 0;

					 if( e[1] <= ( e[0] * 0.00005 ) ) 
						 ier = -1;

					 d = e[2];

					 if( sigma < 0.0 ){
						 d = -d;
						  
						 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
							 ier = -1;
					 }

					 d = ( d + e[1] ) + e[0];

					 rms = ( e0 - d ) - d;

					 if( rms < 0.0 ) 
						 rms = 0.0;

					 return;
				  }

				  b[gj][1] = 0.0;
                  b[gk][1] = -b[gl][0] / p;
                  b[gl][1] = b[gk][0] / p;
			 }
			 else{
				 p = 1.0 / sqrt( p );

				for( int i = 0 ; i < 3 ; i++ )
					b[i][1] = b[i][1] * p;
			 }

			 b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
             b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
             b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];

			 for( int i = 0 ; i < 3 ; i++ ){
				 for( int j = 0 ; j < 3 ; j++ ){
					 u[i][j] = b[i][0] * a[j][0] + b[i][1] * a[j][1] + b[i][2] * a[j][2];
				 }
			 }

			 for( int i = 0 ; i < 3 ; i++ ){
				 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
			 }

			 for( int i = 0 ; i < 3 ; i++ ){
				 if( e[i] < 0 ) 
					 e[i] = 0.0;

				 e[i] = sqrt( e[i] );
			 }

			 ier = 0;

             if( e[1] <= ( e[0] * 0.00005 ) ) 
				 ier = -1;

			 d = e[2];

			 if( sigma < 0.0 ){
				 d = -d;
				  
				 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
					 ier = -1;
			 }

			 d = ( d + e[1] ) + e[0];

			 rms = ( e0 - d ) - d;

			 if( rms < 0.0 ) 
				 rms = 0.0;

			 return;
		 }
	 }
	 
     sqrth = sqrt( h );
	 d = h*h*h - g*g;

	 if( d < 0 ) d = 0.0;

	 d = atan2( sqrt(d), -g ) / 3.0;
     cth = sqrth * cos(d);
     sth = sqrth * sqrt3 * sin(d);
     e[0] = ( spur + cth ) + cth;
     e[1] = ( spur - cth ) + sth;
     e[2] = ( spur - cth ) - sth;

	 if( mode == 0 ){
		 for( int i = 0 ; i < 3 ; i++ ){
			 if( e[i] < 0 ) 
				 e[i] = 0.0;

			 e[i] = sqrt( e[i] );
		 }

		 ier = 0;

		 if( e[1] <= ( e[0] * 0.00005 ) ) 
			 ier = -1;

		 d = e[2];

		 if( sigma < 0.0 ){
			 d = -d;
			  
			 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
				 ier = -1;
		 }

		 d = ( d + e[1] ) + e[0];

		 rms = ( e0 - d ) - d;

		 if( rms < 0.0 ) 
			 rms = 0.0;

		 return;
	 }

     for( int l = 0 ; l < 3 ; l += 2 ){
		 d = e[l];
         ss[0] = ( d - rr[2] ) * ( d - rr[5] ) - rr[4] * rr[4];
         ss[1] = ( d - rr[5] ) * rr[1] + rr[3] * rr[4];
         ss[2] = ( d - rr[0] ) * ( d - rr[5] ) - rr[3] * rr[3];
         ss[3] = ( d - rr[2] ) * rr[3] + rr[1] * rr[4];
         ss[4] = ( d - rr[0] ) * rr[4] + rr[1] * rr[3];
         ss[5] = ( d - rr[0] ) * ( d - rr[2] ) - rr[1] * rr[1];

		 if( abs( ss[0] ) >= abs( ss[2] ) ){
			 gj = 0;

			 if( abs( ss[0] ) < abs( ss[5] ) )
				 gj = 2;
		 }
		 else if( abs( ss[2] ) >= abs( ss[5] ) )
			 gj = 1;
		 else
			 gj = 2;

		 d = 0.0;
         gj = 3 * gj;

		 for( int i = 0 ; i < 3 ; i++ ){
			 gk = ip[i+gj];
             a[i][l] = ss[gk];
             d = d + ss[gk] * ss[gk];
		 }

		 if( d > 0 )
			 d = 1.0 / sqrt(d);

		 for( int i = 0 ; i < 3 ; i++ ){
			 a[i][l] = a[i][l] * d;
		 }
	 }

	 d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];

	 if( ( e[0] - e[1] ) > ( e[1] - e[2] ) ){ 
         gm1 = 2;
         gm = 0;
	 }
	 else{
         gm1 = 0;
         gm = 2;
	 }

	 p = 0.0;

     for( int i = 0 ; i < 3 ; i++ ){
		 a[i][gm1] = a[i][gm1] - d * a[i][gm];
         p = p + pow( a[i][gm1], 2 );
	 }

	 if( p <= tol ){
		  p = 1.0;

		  for( int i = 0 ; i < 3 ; i++ ){
			  if( p < abs( a[i][gm] ) ){}
			  else{
				  p = abs( a[i][gm] );
			      gj = i;
			  }
		  }

		  gk = ip1201[gj];
          gl = ip1201[gj+1];
		  p = sqrt( pow( a[gk][gm], 2 ) + pow( a[gl][gm], 2 ) );

		  if( p <= tol ){
			 for( int i = 0 ; i < 3 ; i++ ){
				 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
			 }

			 for( int i = 0 ; i < 3 ; i++ ){
				 if( e[i] < 0 ) 
					 e[i] = 0.0;

				 e[i] = sqrt( e[i] );
			 }

			 ier = 0;

             if( e[1] <= ( e[0] * 0.00005 ) ) 
				 ier = -1;

			 d = e[2];

			 if( sigma < 0.0 ){
				 d = -d;
				  
				 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
					 ier = -1;
			 }

			 d = ( d + e[1] ) + e[0];

			 rms = ( e0 - d ) - d;

			 if( rms < 0.0 ) 
				 rms = 0.0;

			 return;
		  }

		  a[gj][gm1] = 0.0;
          a[gk][gm1] = -a[gl][gm] / p;
          a[gl][gm1] =  a[gk][gm] / p;
	 }
	 else{
		 p = 1.0 / sqrt(p);

		 for( int i = 0 ; i < 3 ; i++ ){
			 a[i][gm1] = a[i][gm1] * p;
		 }
	 }

	 a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
     a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
     a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];

	 for( int l = 0 ; l < 2 ; l++ ){
		 d = 0.0;

		 for( int i = 0 ; i < 3 ; i++ ){
			 b[i][l] = r[i][0] * a[0][l] + r[i][1] * a[1][l] + r[i][2] * a[2][l];
			 d = d + pow( b[i][l], 2 );
		 }

		 if( d > 0 ) 
			 d = 1.0 / sqrt( d );

		 for( int i = 0 ; i < 3 ; i++ ){
			 b[i][l] = b[i][l] * d;
		 }
	 }

	 d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
	 p = 0.0;

	 for( int i = 0 ; i < 3 ; i++ ){
		 b[i][1] = b[i][1] - d * b[i][0];
		 p = p + pow( b[i][1], 2 );
	 }

	 if( p <= tol ){
		  p = 1.0;

		  for( int i = 0 ; i < 3 ; i++ ){
			  if( p < abs( b[i][0] ) ){}
			  else{
				  p = abs( b[i][0] );
				  gj = i;
			  }
		  }

		  gk = ip1201[gj];
		  gl = ip1201[gj+1];
		  p = sqrt( pow( b[gk][0], 2 ) + pow( b[gl][0], 2 ) );

		  if( p <= tol ){
			 for( int i = 0 ; i < 3 ; i++ ){
				 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
			 }

			 for( int i = 0 ; i < 3 ; i++ ){
				 if( e[i] < 0 ) 
					 e[i] = 0.0;

				 e[i] = sqrt( e[i] );
			 }

			 ier = 0;

             if( e[1] <= ( e[0] * 0.00005 ) ) 
				 ier = -1;

			 d = e[2];

			 if( sigma < 0.0 ){
				 d = -d;
				  
				 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
					 ier = -1;
			 }

			 d = ( d + e[1] ) + e[0];

			 rms = ( e0 - d ) - d;

			 if( rms < 0.0 ) 
				 rms = 0.0;

			 return;
		  }

		  b[gj][1] = 0.0;
		  b[gk][1] = -b[gl][0] / p;
		  b[gl][1] = b[gk][0] / p;
	 }
	 else{
		 p = 1.0 / sqrt( p );

		 for( int i = 0 ; i < 3 ; i++ )
			b[i][1] = b[i][1] * p;
	 }

	 b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
	 b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
	 b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];

	 for( int i = 0 ; i < 3 ; i++ ){
		 for( int j = 0 ; j < 3 ; j++ ){
			 u[i][j] = b[i][0] * a[j][0] + b[i][1] * a[j][1] + b[i][2] * a[j][2];
		 }
	 }

	 for( int i = 0 ; i < 3 ; i++ ){
		 t[i] = yc[i] - u[i][0] * xc[0] - u[i][1] * xc[1] - u[i][2] * xc[2];
	 }

	 for( int i = 0 ; i < 3 ; i++ ){
		 if( e[i] < 0 ) 
			 e[i] = 0.0;

		 e[i] = sqrt( e[i] );
	 }

	 ier = 0;

	 if( e[1] <= ( e[0] * 0.00005 ) ) 
		 ier = -1;

	 d = e[2];

	 if( sigma < 0.0 ){
		 d = -d;
		  
		 if( ( e[1] - e[2] ) <= ( e[0] * 0.00005 ) ) 
			 ier = -1;
	 }

	 d = ( d + e[1] ) + e[0];

	 rms = ( e0 - d ) - d;

	 if( rms < 0.0 ) 
		 rms = 0.0;

	 return;
}  





void setMaxCliqueAlignedNodes( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, vector< Atom * > &vRefAligned, vector< Atom * > &vFitAligned, vector< int > &vPairs ){
	char buffer[100];
	char index[5];
	
	FILE *pInputFile_pairs = fopen( "pairs.rst", "r" );
	int pair_index, res_ref, res_fit;
	int start_pos;

	fgets( buffer, 100, pInputFile_pairs );
	//cout << buffer;

	while( 1 ){
		fgets( buffer, 100, pInputFile_pairs );

		if( buffer[0] == 'E' )
			break;

		int pos = 0;

		for( int i = 0 ; i < 100 ; i++ ){
			if( buffer[i] == '\t' ){
				for( int j = 0 ; j < i ; j++ )
					index[j] = buffer[j];													
				
				index[i] = '\0';
				pair_index = (int)atof( index );
				start_pos = i+1;

				break;
			}		
		}

		for( int i = 0 ; i < vPairs.size() ; i++ ){
			if( pair_index == vPairs[i] ){
				for( int j = start_pos ; j < 100 ; j++ ){
					if( buffer[j] == '\t' ){
						int pos = 0;

						for( int k = start_pos ; k < j ; k++ ){
							index[pos] = buffer[k];		
							pos++;
						}
						
						index[pos] = '\0';
						res_ref = (int)atof( index );
						start_pos = j+1;

						break;
					}
				}		
				
				for( int j = start_pos ; j < 100 ; j++ ){
					if( buffer[j] == '\n' ){
						int pos = 0;

						for( int k = start_pos ; k < j ; k++ ){
							index[pos] = buffer[k];		
							pos++;
						}
						
						index[pos] = '\0';
						res_fit = (int)atof( index );

						break;
					}
				}				

				//cout << pair_index << "\t" << res_ref << "\t" << res_fit << "\n";

				vRefAligned.push_back( vRefStruct[res_ref] );
				vFitAligned.push_back( vFitStruct[res_fit] );
			}
		}
	}

	fclose( pInputFile_pairs );
}





void writeRotationMatrix( string filename ){
	FILE *pOutputFile = fopen( filename.c_str(), "w" );
	
	fprintf( pOutputFile, "*** Rotation matrix to superpose a structure to reference structure ***\n" );
	fprintf( pOutputFile, "i         t(i)           u(i,1)        u(i,2)         u(i,3)\n" );
	
	fprintf( pOutputFile, "1   " );
	fprintf( pOutputFile, "%15.10f", t[0] );
	fprintf( pOutputFile, "%15.10f", u[0][0] );
	fprintf( pOutputFile, "%15.10f", u[0][1] );
	fprintf( pOutputFile, "%15.10f", u[0][2] );
	fprintf( pOutputFile, "\n" );
	
	fprintf( pOutputFile, "2   " );
	fprintf( pOutputFile, "%15.10f", t[1] );
	fprintf( pOutputFile, "%15.10f", u[1][0] );
	fprintf( pOutputFile, "%15.10f", u[1][1] );
	fprintf( pOutputFile, "%15.10f", u[1][2] );
	fprintf( pOutputFile, "\n" );
	
	fprintf( pOutputFile, "3   " );	
	fprintf( pOutputFile, "%15.10f", t[2] );
	fprintf( pOutputFile, "%15.10f", u[2][0] );
	fprintf( pOutputFile, "%15.10f", u[2][1] );
	fprintf( pOutputFile, "%15.10f", u[2][2] );
	fprintf( pOutputFile, "\n" );
	
	fprintf( pOutputFile, "\n" );	
	fprintf( pOutputFile, "X = t[1] + u[1][1]*x + u[1][2]*y + u[1][3]*z\n" );	
	fprintf( pOutputFile, "Y = t[2] + u[2][1]*x + u[2][2]*y + u[2][3]*z\n" );	
	fprintf( pOutputFile, "Z = t[3] + u[3][1]*x + u[3][2]*y + u[3][3]*z\n" );	
	
	fclose( pOutputFile );
}





void setBsetRotationMatrix(){
	for( int i = 0 ; i < 3 ; i++ ){
		 t_best[i] = t[i];
		 
		 u_best[i][0] = u[i][0];
		 u_best[i][1] = u[i][1];
		 u_best[i][2] = u[i][2];
	}
}




// ===============================================================
// ADDED (not upstream): publish the best alignment found so far.
//
// Everything this program writes is produced after the last search iteration finishes, and the
// search can run for hours, so a run in progress is otherwise completely opaque -- it prints
// nothing and touches no file.  This is called once per iteration of each search loop and
// writes the current and best GA-score together with u_best/t_best, which is a complete and
// usable transformation.  The web front end reads it to report progress and to draw the
// alignment as it stands.
//
// It only writes.  No state is modified and nothing is read back, so the resulting GA-score is
// bit-for-bit what the unpatched program produces; infrastructure/smoke_test.py checks that
// against a known value on every build.
// ===============================================================
void writeCheckpoint( const char *stage, int iter, int total, double score, double best,
                      int N_target ){
	FILE *pOutputFile = fopen( "checkpoint.tmp", "w" );
	
	if( pOutputFile == NULL )
		return;
	
	// The scores tracked inside the search loops are unnormalised sums; the GA-score that gets
	// printed at the end is best_score/N_target.  The same division is applied here so the
	// reported number is on the scale callers expect rather than an internal one.
	double norm = ( N_target > 0 ) ? (double)N_target : 1.0;
	
	fprintf( pOutputFile, "stage\t%s\n", stage );
	fprintf( pOutputFile, "iteration\t%d\t%d\n", iter, total );
	fprintf( pOutputFile, "score\t%15.10f\n", score / norm );
	fprintf( pOutputFile, "best_score\t%15.10f\n", best / norm );
	
	for( int i = 0 ; i < 3 ; i++ ){
		fprintf( pOutputFile, "matrix\t%15.10f%15.10f%15.10f%15.10f\n",
		         t_best[i], u_best[i][0], u_best[i][1], u_best[i][2] );
	}
	
	fclose( pOutputFile );
	
	// Renamed into place rather than written in place: the reader polls on a timer and would
	// otherwise sometimes get a half-written file.
	rename( "checkpoint.tmp", "checkpoint.rst" );
	
	
	// The same numbers appended rather than overwritten, which gives the convergence history
	// for free.  This costs one open/append/close per search iteration -- a dozen or so over an
	// entire run, against a search that can take hours -- so it is not worth sampling or
	// buffering.  clock() rather than time() because it starts at program start on its own and
	// has sub-second resolution; the search is single-threaded and CPU-bound, so CPU seconds
	// and wall seconds are the same number here.
	static int cp_logged = 0;
	
	FILE *pLogFile = fopen( "convergence.tsv", "a" );
	
	if( pLogFile == NULL )
		return;
	
	if( cp_logged == 0 ){
		fprintf( pLogFile, "seconds\tstage\titeration\ttotal\tscore\tbest_score\n" );
		cp_logged = 1;
	}
	
	fprintf( pLogFile, "%.2f\t%s\t%d\t%d\t%.10f\t%.10f\n",
	         (double)clock() / CLOCKS_PER_SEC, stage, iter, total, score / norm, best / norm );
	
	fclose( pLogFile );
}





void CopyBestRotationMatrix(){
	for( int i = 0 ; i < 3 ; i++ ){
		 t[i] = t_best[i];
		 
		 u[i][0] = u_best[i][0];
		 u[i][1] = u_best[i][1];
		 u[i][2] = u_best[i][2];
	}
}





// ===============================================================
// This is a function to accurately identify aligned residue pairs
// by solving Linear Sum Assignment Problem (LSAP)
// Reference:
// Derigs,U. (1985) The shortest augmenting path method for solving assignment problems - Motivation and computational experience. In: Monma,C.L. (ed.)
// Algorithms and Software for Optimization. Baltzer, Basel, pp. 57–102.
// 
// Modified by H. S. Lee
// ==============================================================

void assignAlignedResiduePairs( vector< Atom * > vStruct1, vector< Atom * > vStruct2, int align_mode ){
	int n_1 = vStruct1.size();	  // vStruct1.size() >= vStruct2.size();
	int n_2 = vStruct2.size();
	int n;
	
	if( n_1 >= n_2 )
		n = n_1;
	else
		n = n_2;
	
	float d_2 = 0;
	float d0_2 = d0 * d0;
	
	for( int i = 1 ; i <= n_1 ; i++ ){
		for( int j = 1 ; j <= n_2 ; j++ ){
			cost[i][j] = 2;
		}
	}
	
	if( align_mode == 1 ){
		for( int i = 0 ; i < n_1 ; i++ ){
			for( int j = 0 ; j < n_2 ; j++ ){
				d_2 = pow( vStruct1[i]->R[0] - vStruct2[j]->R_sup[0], 2 ) + pow( vStruct1[i]->R[1] - vStruct2[j]->R_sup[1], 2 ) + pow( vStruct1[i]->R[2] - vStruct2[j]->R_sup[2], 2 );
            	score[i+1][j+1] = 1 / ( 1 + d_2 / d0_2 );            	
        	}
		}
	}	
	else if( align_mode == 2 ){
		for( int i = 0 ; i < n_1 ; i++ ){
			for( int j = 0 ; j < n_2 ; j++ ){
				d_2 = pow( vStruct1[i]->R_sup[0] - vStruct2[j]->R[0], 2 ) + pow( vStruct1[i]->R_sup[1] - vStruct2[j]->R[1], 2 ) + pow( vStruct1[i]->R_sup[2] - vStruct2[j]->R[2], 2 );
            	score[i+1][j+1] = 1 / ( 1 + d_2 / d0_2 );            	
        	}
		}
	}
	
		
	for( int i = 1 ; i <= n_1 ; i++ ){
		for( int j = 1 ; j <= n_2 ; j++ ){
            cost[i][j] = 2 - score[i][j];
        }
	}
	

	// **********************************************
	//   solve Linear Sum Assignment Problem (LSAP)
	// **********************************************
	
	double eps = 1E-10;
	double sup = 1E+9;
	double cc, d, ui, vgl, vj, z;
	double dminus[3000];
	double dplus[3000];
	double ys[3000];
	double yt[3000];
	int i, j, u, w, ind, index, j0;
	int zeile[3000];
	int vor[3000];
	int label[3000];
	
	
	// construct an initial partial assignment
	for( i = 1 ; i <= n ; i++ ){
		zeile[i] = 0;
		spalte[i] = 0;
		vor[i] = 0;
		ys[i] = 0.0E+00;
		yt[i] = 0.0E+00;
	}
	
	for( i = 1 ; i <= n_2 ; i++ ){
		for( j = 1 ; j <= n_1 ; j++ ){
			cc = cost[j][i];
			
			if( j != 1 ){
				if( ( cc - ui ) >= eps ){
				}
				else{
					ui = cc;
					j0 = j;
				}
			}
			else{
				ui = cc;
				j0 = j;
			}
		}	
			
		ys[i] = ui;
					
		if( zeile[j0] == 0 ){
			zeile[j0] = i;
			spalte[i] = j0;
		}
	}
	
	for( j = 1 ; j <= n_1 ; j++ ){
		yt[j] = 0;
		
		if( zeile[j] == 0 ){
			yt[j] = sup; 
		}
	}
	
	
	for( i = 1 ; i <= n_2 ; i++ ){
		ui = ys[i];
		
		for( j = 1 ; j <= n_1 ; j++ ){
			vj = yt[j];
			
			if( eps < vj ){
				cc = cost[j][i] - ui;
				
				if( cc + eps < vj ){				
					yt[j] = cc;
					vor[j] = i;
				}
			}
		}
	}
	
	
	for( j = 1 ; j <= n_1 ; j++ ){
		i = vor[j];
		
		if( i != 0 ){
			if( spalte[i] == 0 ){
				spalte[i] = j;
                zeile[j] = i;
			}
		}		
	}
	
	
	for( i = 1 ; i <= n_2 ; i++ ){
		if( spalte[i] == 0 ){
			ui = ys[i];
			
			for( j = 1 ; j <= n_1 ; j++ ){
				if ( zeile[j] == 0 ){
					cc = cost[j][i];
              		
              		if ( cc - ui - yt[j] <= -eps ){
              			spalte[i] = j;
               			zeile[j] = i;
              		} 
				}
			}		
		}
	}
	
	
	for( u = 1 ; u <= n_2 ; u++ ){
		if( spalte[u] == 0 ){
			for( i = 1 ; i <= n_1 ; i++ ){
				vor[i] = u;
          		label[i] = 0;
          		dplus[i] = sup;
          		dminus[i] = cost[i][u] - ys[u] - yt[i];
			}
			
			dplus[u] = 0.0E+00;
			
			while( 1 ){
				d = sup;
				index = 0;
			
				for( i = 1 ; i <= n_1 ; i++ ){
					if( label[i] == 0 ){
						if( dminus[i] + eps < d ){
							d = dminus[i];
              				index = i;
						}
					}
				}
			
				if( index == 0 ){
					//cout << "LSAPR - Fatal error! No unlabeled node with DMINUS" << "\n";
            		return;
				}
			
				if( zeile[index] <= 0 ){
					while( 1 ){
						w = vor[index];
						zeile[index] = w;
        				ind = spalte[w];
        				spalte[w] = index;
        			
        				if( w != u ){
          					index = ind;
        				}
        				else{
        					break;
        				}
					}
				
					for( i = 1 ; i <= n_2 ; i++ ){
          				if( dplus[i] != sup ){
           					ys[i] = ys[i] + d - dplus[i];
          				}

          				if( dminus[i] + eps < d ){
            				yt[i] = yt[i] + dminus[i] - d;
          				}
        			}
        			
      				break;
				}
				else{
					label[index] = 1;
        			w = zeile[index];
        			dplus[w] = d;

        			for( i = 1 ; i <= n_1 ; i++ ){
          				if( label[i] == 0 ){
            				vgl = d + cost[i][w] - ys[w] - yt[i];
            			
            				if( vgl + eps < dminus[i] ){
              					dminus[i] = vgl;
              					vor[i] = w;
              				}
            			}
          			}
				}							
			}
		}
	}	
}





void executeAssignAlignedResiduePairs( vector< Atom * > vRefStruct, vector< Atom * > vFitStruct, int align_mode ){
	if( align_mode == 1 )
		assignAlignedResiduePairs( vRefStruct, vFitStruct, align_mode );	
	else if( align_mode == 2 )
		assignAlignedResiduePairs( vFitStruct, vRefStruct, align_mode );	
}		
				
			
				

void setAlignedResiduePairs( int align_mode, int n_ref, int n_fit ){
	if( align_mode == 1 ){
		for( int i = 0 ; i < n_ref ; i++ ){
			alignedRes[i] = -1;
		}
	
		for( int i = 1 ; i <= n_fit ; i++ ){
			if( spalte[i] != 0 ){
				alignedRes[spalte[i]-1] = i-1;
			}
		}
	}	
	else if( align_mode == 2 ){
		for( int i = 0 ; i < n_ref ; i++ ){
			alignedRes[i] = -1;
		}
		
		for( int i = 1 ; i <= n_ref ; i++ ){
			if( spalte[i] != 0 ){
				alignedRes[i-1] = spalte[i]-1;
			}
		}	
	}
}




void filterAlignedResiduePairs( vector< Atom * > vRefStruct, vector< Atom * > vFitStruct, float cutoff ){
	float dist;
	
	for( int i = 0 ; i < vRefStruct.size() ; i++ ){
		if( alignedRes[i] != -1 ){
			dist = sqrt( pow( vRefStruct[i]->R[0] - vFitStruct[alignedRes[i]]->R_sup[0], 2 ) + pow( vRefStruct[i]->R[1] - vFitStruct[alignedRes[i]]->R_sup[1], 2 ) + pow( vRefStruct[i]->R[2] - vFitStruct[alignedRes[i]]->R_sup[2], 2 ) );
		
		    if( dist > cutoff ){
		    	alignedRes[i] = -1;
		    }
		}
	}
}



float getCoverage( int n_ref ){
	int num_ali = 0;
	
	for( int i = 0 ; i < n_ref ; i++ ){
		if( alignedRes[i] != -1 )
			num_ali++;			
	}
	
	return num_ali/(float)n_ref;
}




void setAlignedResiduePairsWithBestScore( int n_ref ){
	for( int i = 0 ; i < n_ref ; i++ ){
		alignedResBest[i] = alignedRes[i];		
	}		
}




void showAlignedResiduePairs( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	cout << "\n" << "# s1 - s2" << "\n";
	
	for( int i = 0 ; i < vRefStruct.size() ; i++ ){
		if( alignedResBest[i] != -1 )
			cout << vRefStruct[i]->res_Seq << " - " << vFitStruct[alignedResBest[i]]->res_Seq << "\n";
	}	
}




void showAlignedChemicalFeaturePairs( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	cout << "\n" << "# s1cf - s2cf" << "\n";
	
	for( int i = 0 ; i < vRefStruct.size() ; i++ ){
		if( alignedResBest[i] != -1 )
			cout << vRefStruct[i]->res_Seq << " (" << vRefStruct[i]->cf[0] << vRefStruct[i]->cf[1] << ") - " << vFitStruct[alignedResBest[i]]->res_Seq << " (" << vFitStruct[alignedResBest[i]]->cf[0] << vFitStruct[alignedResBest[i]]->cf[1] << ")" << "\n";
	}	
}




void setFragment( vector< Atom * > &vStruct, vector< Fragment * > &vFragment ){
	Fragment *frag;
	
	int n = vStruct.size();
	
	for( int i = 0 ; i < n-2 ; i++ ){
		frag = new Fragment();
		frag->index[0] = i;
		frag->index[1] = i+1;	
		frag->index[2] = i+2;	
					
		vFragment.push_back( frag );	
	}
}




void setFragmentUsingElementsInProductGraph ( vector< Atom * > &vStruct, vector< Fragment * > &vFragment, int mode ){
	Fragment *frag;
	
	int n = vStruct.size();
	
	for( int i = 0 ; i < n-2 ; i++ ){
		if( mode == 1 ){
			if( ( ref_elements_in_product_graph[i] == 1 ) && ( ref_elements_in_product_graph[i+1] == 1 ) && ( ref_elements_in_product_graph[i+2] == 1 ) ){
				frag = new Fragment();
				frag->index[0] = i;
				frag->index[1] = i+1;	
				frag->index[2] = i+2;	
					
				vFragment.push_back( frag );
			}
		}
		else{
			if( ( fit_elements_in_product_graph[i] == 1 ) && ( fit_elements_in_product_graph[i+1] == 1 ) && ( fit_elements_in_product_graph[i+2] == 1 ) ){
				frag = new Fragment();
				frag->index[0] = i;
				frag->index[1] = i+1;	
				frag->index[2] = i+2;	
					
				vFragment.push_back( frag );
			}		
		}
	}
}




void setNonOverlappingFragment( vector< Atom * > &vStruct, vector< Fragment * > &vFragment ){
	Fragment *frag;
	
	int n = vStruct.size();
	
	for( int i = 0 ; i < n-2 ; i += 3 ){
		frag = new Fragment();
		frag->index[0] = i;
		frag->index[1] = i+1;	
		frag->index[2] = i+2;	
					
		vFragment.push_back( frag );	
	}
}




int checkFragmentPairSimilarity( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, Fragment * ref_frag, Fragment * fit_frag ){
	int num_similar_residues = 0;
	
	if( ss_option == 0 ){
		for( int i = 0 ; i < 3 ; i++ ){
			if( blosum62[getValueBLOSUN62Map( vRefStruct[ref_frag->index[i]]->res_name )][getValueBLOSUN62Map( vFitStruct[fit_frag->index[i]]->res_name )] >= min_blosum_value )	
				num_similar_residues++;
		}
	}
	else{
		for( int i = 0 ; i < 3 ; i++ ){
			if( ( blosum62[getValueBLOSUN62Map( vRefStruct[ref_frag->index[i]]->res_name )][getValueBLOSUN62Map( vFitStruct[fit_frag->index[i]]->res_name )] >= min_blosum_value ) && ( vRefStruct[ref_frag->index[i]]->ss == vFitStruct[fit_frag->index[i]]->ss ) )	
				num_similar_residues++;
		}
	}
	
	if( num_similar_residues >= f_option )
		return 1;
	else
		return 0;
}




void clearFragmentVector( vector< Fragment * > &vFragment ){
	for( int i = 0 ; i < vFragment.size() ; i++ ){
		delete vFragment[i];
	}

	vFragment.clear();			
	vFragment.resize( 0 ); 
}





void setAlignedAtomPairsUsingFragments( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, Fragment * ref_frag, Fragment * fit_frag, vector< Atom * > &vRefAligned, vector< Atom * > &vFitAligned ){
	for( int i = 0 ; i < 3 ; i++ ){
		vRefAligned.push_back( vRefStruct[ref_frag->index[i]] );		
		vFitAligned.push_back( vFitStruct[fit_frag->index[i]] );		
	}			
}





void setAlignedAtomPairsUsingAlignedResArray( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, vector< Atom * > &vRefAligned, vector< Atom * > &vFitAligned ){
	Atom * atom;
	int n_ref = vRefStruct.size();
	
	for( int i = 0 ; i < n_ref ; i++ ){
		if( alignedRes[i] != -1 ){
			atom = new Atom();				
			atom->R[0] = vRefStruct[i]->R[0];
			atom->R[1] = vRefStruct[i]->R[1];
			atom->R[2] = vRefStruct[i]->R[2];				
			vRefAligned.push_back( atom );
				
			atom = new Atom();				
			atom->R[0] = vFitStruct[alignedRes[i]]->R[0];
			atom->R[1] = vFitStruct[alignedRes[i]]->R[1];
			atom->R[2] = vFitStruct[alignedRes[i]]->R[2];				
			vFitAligned.push_back( atom );
		}
	}	
}






double getGLoSAScore( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	int n_ref = vRefStruct.size();	 
	int n_fit = vFitStruct.size();
	
	double d0 = 2.5;
 	double score = 0; 
 	double d_R_2, d_SR_2;
 	char ref_cf[2];
 	char fit_cf[2];
 		
	double d0_2 = d0 * d0;
	double q = 0;
	
	for( int i = 0 ; i < n_ref ; i++ ){
		if( alignedRes[i] != -1 ){
			d_R_2 = pow( vRefStruct[i]->R[0] - vFitStruct[alignedRes[i]]->R_sup[0], 2 ) 
					+ pow( vRefStruct[i]->R[1] - vFitStruct[alignedRes[i]]->R_sup[1], 2 ) 
					+ pow( vRefStruct[i]->R[2] - vFitStruct[alignedRes[i]]->R_sup[2], 2 );
					
			ref_cf[0] = vRefStruct[i]->cf[0];
			ref_cf[1] = vRefStruct[i]->cf[1];
			
			fit_cf[0] = vFitStruct[alignedRes[i]]->cf[0];
			fit_cf[1] = vFitStruct[alignedRes[i]]->cf[1];
			
			
			// chemical feature similarity
			if( ( ref_cf[0] == fit_cf[0] ) && ( ref_cf[1] == fit_cf[1] ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'D' ) && ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) )
				q = 1;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'D' ) && ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'A' ) && ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) )
				q = 1;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'A' ) && ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'D' ) && ( fit_cf[0] == 'P' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'D' ) && ( ref_cf[0] == 'P' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'A' ) && ( fit_cf[0] == 'N' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'A' ) && ( ref_cf[0] == 'N' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) && ( fit_cf[0] == 'P' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) && ( ref_cf[0] == 'P' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) && ( fit_cf[0] == 'N' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) && ( ref_cf[0] == 'N' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'A' ) && ( ref_cf[1] == 'R' ) && ( fit_cf[0] == 'A' ) && ( fit_cf[1] == 'L' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'A' ) && ( fit_cf[1] == 'R' ) && ( ref_cf[0] == 'A' ) && ( ref_cf[1] == 'L' ) )
				q = 0.8;
			else
				q = 0.5;			
			
					
			score += q / ( 1 + ( d_R_2 / d0_2 ) );   
		}
	}	
	
	return score;
}






double getNormalizedGLoSAScore( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct, int N_target ){
	int n_ref = vRefStruct.size();	 
	int n_fit = vFitStruct.size();
	
	double d0;
	double p = 1.0/2.0;
	
	// *** universal purpose ***
	d0 = 0.27 * pow( N_target - 6.0, p ) + 0.98;
		
 	double score = 0; 
 	double d_R_2, d_SR_2;
 	char ref_cf[2];
 	char fit_cf[2];
 		
	double d0_2 = d0 * d0;
	double q = 0;
	
	for( int i = 0 ; i < n_ref ; i++ ){
		if( alignedRes[i] != -1 ){
			d_R_2 = pow( vRefStruct[i]->R[0] - vFitStruct[alignedRes[i]]->R_sup[0], 2 ) 
					+ pow( vRefStruct[i]->R[1] - vFitStruct[alignedRes[i]]->R_sup[1], 2 ) 
					+ pow( vRefStruct[i]->R[2] - vFitStruct[alignedRes[i]]->R_sup[2], 2 );
					
			ref_cf[0] = vRefStruct[i]->cf[0];
			ref_cf[1] = vRefStruct[i]->cf[1];
			
			fit_cf[0] = vFitStruct[alignedRes[i]]->cf[0];
			fit_cf[1] = vFitStruct[alignedRes[i]]->cf[1];
			
			
			// chemical feature similarity
			if( ( ref_cf[0] == fit_cf[0] ) && ( ref_cf[1] == fit_cf[1] ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'D' ) && ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) )
				q = 1;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'D' ) && ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'A' ) && ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) )
				q = 1;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'A' ) && ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) )
				q = 1;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'D' ) && ( fit_cf[0] == 'P' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'D' ) && ( ref_cf[0] == 'P' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'H' ) && ( ref_cf[1] == 'A' ) && ( fit_cf[0] == 'N' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'H' ) && ( fit_cf[1] == 'A' ) && ( ref_cf[0] == 'N' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) && ( fit_cf[0] == 'P' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) && ( ref_cf[0] == 'P' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'O' ) && ( ref_cf[1] == 'H' ) && ( fit_cf[0] == 'N' ) && ( fit_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'O' ) && ( fit_cf[1] == 'H' ) && ( ref_cf[0] == 'N' ) && ( ref_cf[1] == 'C' ) )
				q = 0.8;
			else if( ( ref_cf[0] == 'A' ) && ( ref_cf[1] == 'R' ) && ( fit_cf[0] == 'A' ) && ( fit_cf[1] == 'L' ) )
				q = 0.8;
			else if( ( fit_cf[0] == 'A' ) && ( fit_cf[1] == 'R' ) && ( ref_cf[0] == 'A' ) && ( ref_cf[1] == 'L' ) )
				q = 0.8;
			else
				q = 0.5;			
			
					
			score += q / ( 1 + ( d_R_2 / d0_2 ) );   
		}
	}	
	
	return score;
}






// *** the scoring function of previous version of G-LoSA ***
double getScore( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	int n_ref = vRefStruct.size();	

	int N_align = 0;
	double rmsd = 0;
	double temp = 0;

	for( int i = 0 ; i < vRefStruct.size() ; i++ ){
		if( alignedRes[i] != -1 ){
			N_align++;

			rmsd += pow( vRefStruct[i]->R[0] - vFitStruct[alignedRes[i]]->R_sup[0], 2 ) 
				+ pow( vRefStruct[i]->R[1] - vFitStruct[alignedRes[i]]->R_sup[1], 2 ) 
				+ pow( vRefStruct[i]->R[2] - vFitStruct[alignedRes[i]]->R_sup[2], 2 );

			rmsd += pow( vRefStruct[i]->SR[0] - vFitStruct[alignedRes[i]]->SR_sup[0], 2 ) 
				+ pow( vRefStruct[i]->SR[1] - vFitStruct[alignedRes[i]]->SR_sup[1], 2 ) 
				+ pow( vRefStruct[i]->SR[2] - vFitStruct[alignedRes[i]]->SR_sup[2], 2 );
		}
	}
	
	rmsd = sqrt( rmsd / ( N_align * 2 ) );
	
	if( N_align == 0 ){
		return 0;
	}
	else{
		if( rmsd <= 0.5 )   
			temp = 0.5;
		else
			temp = rmsd;
			
		return ( N_align * N_align ) / temp;
	}
}





double getAvgDistance( vector< Atom * > &vRefStruct, vector< Atom * > &vFitStruct ){
	int n_ref = vRefStruct.size();	
	int N_align = 0; 
	
	double d = 0;
	double avg_d = 0;
	
	for( int i = 0 ; i < n_ref ; i++ ){
		if( alignedRes[i] != -1 ){
			N_align++;
			
			d += sqrt( pow( vRefStruct[i]->R[0] - vFitStruct[alignedRes[i]]->R_sup[0], 2 ) 
					+ pow( vRefStruct[i]->R[1] - vFitStruct[alignedRes[i]]->R_sup[1], 2 ) 
					+ pow( vRefStruct[i]->R[2] - vFitStruct[alignedRes[i]]->R_sup[2], 2 ) );
		}
	}	
	
	avg_d = d / (double)N_align;		

	return avg_d;
}





void clearEgdeVector( vector< Edge * > &vEdge ){
	for( int i = 0 ; i < vEdge.size() ; i++ )	    
		delete vEdge[i];
				
	vEdge.clear();			
	vEdge.resize( 0 ); 
}




void clearPairsVector( vector< int > &vPairs ){
	vPairs.clear();			
	vPairs.resize( 0 ); 	
}



void clearPairVector( vector< Pair * > &vPair ){
	for( int i = 0 ; i < vPair.size() ; i++ )	    
		delete vPair[i];

	vPair.clear();			
	vPair.resize( 0 ); 
}


void clearAlignedResidueVector( vector< Atom * > &vAligned ){
	vAligned.clear();			
	vAligned.resize( 0 ); 
}


void clearAtomVector( vector< Atom * > &vStruct ){
	for( int i = 0 ; i < vStruct.size() ; i++ )	    
		delete vStruct[i];

	vStruct.clear();			
	vStruct.resize( 0 ); 	
}

// ************************








// *********************************
// ********** G-LoSA MAIN **********
// *********************************

int main( int argc,char *argv[] ){
	//clock_t begin = clock();
	
	// *** Parsing Arguments & Setup Initial Parameters ***
	parseCommandLine( argc, argv ); 
	setBLOSUN62Map();


	// ============================
	//    Generate product graph 
	// ============================	
	
	double score = 0;
	double best_score = 0;
	double score_temp = 0;
	double best_score_temp = 0;
	
	double best_dist_tolerance = 1.0;
	double avg_dist = 0;
	double best_avg_dist = 0;
	
	vector< Atom * > vRefStruct; 
	vector< Atom * > vRefStructCF; 
	vector< Atom * > vFitStruct;		
	vector< Atom * > vFitStructCF;		
	vector< Atom * > vRefSurface; 
	vector< Atom * > vRefSurfaceCF;
	
	vector< Pair * > vPair;	
	vector< Edge * > vEdge;
	
	vector< int > vPairs;
	vector< Atom * > vRefAligned;
	vector< Atom * > vFitAligned;
	
	vector< Fragment * > vRefFragment;
	vector< Fragment * > vFitFragment;
	
    int n_ref, n_fit, n_ref_cf, n_fit_cf;
    int align_mode = 0;
    int align_mode_cf = 0;	
    int N_target, N_target_min, N_target_max;
    double score_ref, score_min;

        
    
	
	
	// =============================
	// *** read input structures ***
	// =============================
	
	readCAAtomWithSideChainCentroidFromPDBFile( s1, vRefStruct );
	readChemicalFeatures( s1cf, vRefStructCF );
	
	readCAAtomWithSideChainCentroidFromPDBFile( s2, vFitStruct );
	readChemicalFeatures( s2cf, vFitStructCF );
	
	if( ( surf1 != "n" ) && ( surf1cf != "n" ) && ( surf_option == 1 ) ){
		readCAAtomWithSideChainCentroidFromPDBFile( surf1, vRefSurface );
		readChemicalFeatures( surf1cf, vRefSurfaceCF );
	}
	
	if( ( s1ss != "n" ) && ( s2ss != "n" ) && ( ss_option = 1 ) ){
		setSecondaryStructure( s1ss, vRefStruct );
		setSecondaryStructure( s2ss, vFitStruct );
	}
	
	n_ref = vRefStruct.size();
	n_fit = vFitStruct.size();
	
	n_ref_cf = vRefStructCF.size();
	n_fit_cf = vFitStructCF.size();
	
	
	if( n_ref >= n_fit )
		align_mode = 1;
	else
		align_mode = 2;
	
	if( n_ref_cf >= n_fit_cf ){
		align_mode_cf = 1;
		N_target_min = n_fit_cf;
		N_target_max = n_ref_cf;
	}
	else{
		align_mode_cf = 2;
		N_target_min = n_ref_cf;
		N_target_max = n_fit_cf;
	}
	
	
	if( n_option == 1 )
		N_target = N_target_min;
	else if( n_option == 2 )
		N_target = N_target_max;	
	else if( n_option == 3 )
		N_target = ( N_target_min + N_target_max ) / 2;
	else if( n_option == 4 )
		N_target = n_ref_cf;
	else if( n_option == 5 )
		N_target = n_fit_cf;
	
		
   
			
			
			
	
	
	// ==========================================================
	// *** Iterative maximum clique-based structure alignment ***
	// ==========================================================
	
	// *** Iterative MCS using Ca atom coordinates ***
	if( surf_option == 0 )
		GeneratePairs( vRefStruct, vFitStruct, vPair );
	else
		GeneratePairs( vRefSurface, vFitStruct, vPair );
		
	writePairs( vPair );	
	
	if( num_iter > 3 )
		num_iter = 3;
	
	cp_stage = "CA clique search";	cp_total = num_iter;	// ADDED
	
	for( int i = 0 ; i < num_iter ; i++ ){
	    dist_tolerance = dist_tolerance_init + i * 2.0;
	    
	    if( surf_option == 0 )
			GenerateEdges( vRefStruct, vFitStruct, vPair, vEdge, dist_tolerance );
		else
			GenerateEdges( vRefSurface, vFitStruct, vPair, vEdge, dist_tolerance );
			
		writeProductGraph( vPair.size(), vEdge );
		
		int max_clique_result = identifyMaxClique( vPairs );
	    
		if( max_clique_result == 0 )
			score = 0;
		else{
			if( surf_option == 0 )		
				setMaxCliqueAlignedNodes( vRefStruct, vFitStruct, vRefAligned, vFitAligned, vPairs );
			else
				setMaxCliqueAlignedNodes( vRefSurface, vFitStruct, vRefAligned, vFitAligned, vPairs );
			
			setRotationMatrix( vRefAligned, vFitAligned );	
			setSuperposedCoordinatesForFitStructure( vFitStruct );
			
			executeAssignAlignedResiduePairs( vRefStruct, vFitStruct, align_mode );
			setAlignedResiduePairs( align_mode, n_ref, n_fit );	
			filterAlignedResiduePairs( vRefStruct, vFitStruct, 8.0 );
			clearAlignedResidueVector( vRefAligned );
			clearAlignedResidueVector( vFitAligned );		
			
			setAlignedAtomPairsUsingAlignedResArray( vRefStruct, vFitStruct, vRefAligned, vFitAligned );
			setRotationMatrix( vRefAligned, vFitAligned );	
			
			
			setSuperposedCoordinatesForFitStructure( vFitStructCF );
			executeAssignAlignedResiduePairs( vRefStructCF, vFitStructCF, align_mode_cf );
			setAlignedResiduePairs( align_mode_cf, n_ref_cf, n_fit_cf );		
			//avg_dist = getAvgDistance( vRefStructCF, vFitStructCF );
			//score = getGLoSAScore( vRefStructCF, vFitStructCF );
			score = getNormalizedGLoSAScore( vRefStructCF, vFitStructCF, N_target );
		}
		
		if( score > best_score ){
			best_dist_tolerance = dist_tolerance;			
			best_score = score;
			
			setBsetRotationMatrix();	
			setAlignedResiduePairsWithBestScore( n_ref_cf );
			
			//best_avg_dist = avg_dist;
		}	
		
		writeCheckpoint( cp_stage, i + 1, cp_total, score, best_score, N_target );	// ADDED
		
		clearEgdeVector( vEdge );
		clearPairsVector( vPairs );		
		clearAlignedResidueVector( vRefAligned );
		clearAlignedResidueVector( vFitAligned );		
	}
		
	clearPairVector( vPair );
	
	
	
	
	
	
	
	// *** Iterative MCS using CFP coordinates *** 
	if( surf_option == 0 )
		GeneratePairsForCF( vRefStructCF, vFitStructCF, vPair );
	else
		GeneratePairsForCF( vRefSurfaceCF, vFitStructCF, vPair );
	
	writePairs( vPair );
	
	if( num_iter_cf > 3 )
		num_iter_cf = 3;
		
	cp_stage = "feature clique search";	cp_total = num_iter_cf;	// ADDED
	
	for( int i = 0 ; i < num_iter_cf ; i++ ){
	    dist_tolerance_cf = dist_tolerance_cf_init + i * 0.5;
	    
	    if( surf_option == 0 )
			GenerateEdges( vRefStructCF, vFitStructCF, vPair, vEdge, dist_tolerance_cf );
		else
			GenerateEdges( vRefSurfaceCF, vFitStructCF, vPair, vEdge, dist_tolerance_cf );
		
		writeProductGraph( vPair.size(), vEdge );
		
		int max_clique_result = identifyMaxClique( vPairs );
	    
		if( max_clique_result == 0 ){
			score = 0;
		}
		else{		
			if( surf_option == 0 )		
				setMaxCliqueAlignedNodes( vRefStructCF, vFitStructCF, vRefAligned, vFitAligned, vPairs );
			else
				setMaxCliqueAlignedNodes( vRefSurfaceCF, vFitStructCF, vRefAligned, vFitAligned, vPairs );
		
			setRotationMatrix( vRefAligned, vFitAligned );	
			setSuperposedCoordinatesForFitStructure( vFitStructCF );
			
			executeAssignAlignedResiduePairs( vRefStructCF, vFitStructCF, align_mode_cf );
			setAlignedResiduePairs( align_mode_cf, n_ref_cf, n_fit_cf );	
			filterAlignedResiduePairs( vRefStructCF, vFitStructCF, 8.0 );
			clearAlignedResidueVector( vRefAligned );
			clearAlignedResidueVector( vFitAligned );		
			
			setAlignedAtomPairsUsingAlignedResArray( vRefStructCF, vFitStructCF, vRefAligned, vFitAligned );
			setRotationMatrix( vRefAligned, vFitAligned );	
			
			
			setSuperposedCoordinatesForFitStructure( vFitStructCF );
			executeAssignAlignedResiduePairs( vRefStructCF, vFitStructCF, align_mode_cf );
			setAlignedResiduePairs( align_mode_cf, n_ref_cf, n_fit_cf );		
			//avg_dist = getAvgDistance( vRefStructCF, vFitStructCF );
			//score = getGLoSAScore( vRefStructCF, vFitStructCF );
			score = getNormalizedGLoSAScore( vRefStructCF, vFitStructCF, N_target );
		}
		
		if( score > best_score ){
			best_dist_tolerance = dist_tolerance;			
			best_score = score;
			
			setBsetRotationMatrix();	
			setAlignedResiduePairsWithBestScore( n_ref_cf );
			
			//best_avg_dist = avg_dist;
		}	
		
		writeCheckpoint( cp_stage, i + 1, cp_total, score, best_score, N_target );	// ADDED
		
		clearEgdeVector( vEdge );
		clearPairsVector( vPairs );		
		clearAlignedResidueVector( vRefAligned );
		clearAlignedResidueVector( vFitAligned );		
	}
		
	clearPairVector( vPair );
	
	
	
	
	
	
	// =======================================================
	// ***************** Fragment alignments *****************
	// =======================================================
	
	if( ( f_option >= 0 ) && ( f_option <= 3 ) ){
		if( n_ref > n_fit ){
			if( surf_option == 0 )
				setNonOverlappingFragment( vRefStruct, vRefFragment );
			else
				setNonOverlappingFragment( vRefSurface, vRefFragment );
				
			setFragment( vFitStruct, vFitFragment );
		}
		else{
			if( surf_option == 0 )
				setFragment( vRefStruct, vRefFragment );
			else
				setFragment( vRefSurface, vRefFragment );
				
			setNonOverlappingFragment( vFitStruct, vFitFragment );
		}
		
		int doFS = 0;
		
		for( int i = 0 ; i < vRefFragment.size() ; i++ ){
			for( int j = 0 ; j < vFitFragment.size() ; j++ ){
				doFS = 0;
				
				if( surf_option == 0 ){
					if( checkFragmentPairSimilarity( vRefStruct, vFitStruct, vRefFragment[i], vFitFragment[j] ) == 1 ){
						doFS = 1;
						setAlignedAtomPairsUsingFragments( vRefStruct, vFitStruct, vRefFragment[i], vFitFragment[j], vRefAligned, vFitAligned );
					}
				}
				else{
					if( checkFragmentPairSimilarity( vRefSurface, vFitStruct, vRefFragment[i], vFitFragment[j] ) == 1 ){
						doFS = 1;
						setAlignedAtomPairsUsingFragments( vRefSurface, vFitStruct, vRefFragment[i], vFitFragment[j], vRefAligned, vFitAligned );
					}
				}
				
				if( doFS == 1 ){				
					setRotationMatrix( vRefAligned, vFitAligned );	
					setSuperposedCoordinatesForFitStructure( vFitStruct );
			
					executeAssignAlignedResiduePairs( vRefStruct, vFitStruct, align_mode );
					setAlignedResiduePairs( align_mode, n_ref, n_fit );			
					filterAlignedResiduePairs( vRefStruct, vFitStruct, 8.0 );
					clearAlignedResidueVector( vRefAligned );
					clearAlignedResidueVector( vFitAligned );
			
					setAlignedAtomPairsUsingAlignedResArray( vRefStruct, vFitStruct, vRefAligned, vFitAligned );
					setRotationMatrix( vRefAligned, vFitAligned );	
					
					setSuperposedCoordinatesForFitStructure( vFitStructCF );
					executeAssignAlignedResiduePairs( vRefStructCF, vFitStructCF, align_mode_cf );
					setAlignedResiduePairs( align_mode_cf, n_ref_cf, n_fit_cf );	
					//avg_dist = getAvgDistance( vRefStructCF, vFitStructCF );	
					//score = getGLoSAScore( vRefStructCF, vFitStructCF );
					score = getNormalizedGLoSAScore( vRefStructCF, vFitStructCF, N_target );
					
					clearAlignedResidueVector( vRefAligned );
					clearAlignedResidueVector( vFitAligned );	
			
					if( score > best_score ){
						best_score = score;
			
						setBsetRotationMatrix();	
						setAlignedResiduePairsWithBestScore( n_ref_cf );
						
						//best_avg_dist = avg_dist;
					}	
				}
			}
			
			// ADDED: once per reference fragment rather than per pair, which would be
			// thousands of writes on a large structure.
			writeCheckpoint( "fragment superposition", i + 1, vRefFragment.size(),
			                 score, best_score, N_target );
		}
	    
		clearFragmentVector( vRefFragment );
		clearFragmentVector( vFitFragment );	
	}
	
	
	
	
	
	// ==============
	// *** Output ***
	// ==============
	
	cout << "\n";
	cout << "Structure 1 (s1): " << s1 << "\n";
	cout << "Structure 2 (s2): " << s2 << "\n";
	cout << "\n";
	
	cout << "GA-score: " << ( 1 / (double)N_target ) * best_score << "\n";
	cout << "\n";
	
	
	// *** only for random set
	//cout << ( 1 / (double)N_target ) * best_score << "\t" << best_avg_dist << "\n";
	
	
	if( output_option == 1 ){
		CopyBestRotationMatrix(); 
		writeRotationMatrix( "matrix.txt" );  
		superposeFitStructure( s2, "ali_struct.pdb" );	
		
		cout << "*** Output Files ***" << "\n";
		cout << "1. ali_struct.pdb     : the PDB coordinates of s2 aligned onto s1" << "\n";
		cout << "2. matrix.txt         : translational and rotational matrix to align s2 onto s1." << "\n";
			 
		if( s2w != "n" ){
			superposeFitStructure( s2w, "ali_struct_with.pdb" );
			cout << "3. ali_struct_with.pdb: the PDB coordinates of s2w structure transferred by the matrix" << "\n";	
		}
	}
	else if( output_option == 2 ){
		CopyBestRotationMatrix(); 
		writeRotationMatrix( "matrix.txt" );  
	}
	

	
	


	// ====== Delete remaining pointers ======
	clearAtomVector( vRefStruct );	
    clearAtomVector( vFitStruct );
    clearAtomVector( vRefStructCF );	
    clearAtomVector( vFitStructCF );
    clearAtomVector( vRefSurface );	
    clearAtomVector( vRefSurfaceCF );
    // =======================================
		
	
    
    //clock_t end = clock();
    //double elapsed_secs = double( end - begin) / CLOCKS_PER_SEC;
    //cout << elapsed_secs << "\n";
  
    
	return 0;
}
