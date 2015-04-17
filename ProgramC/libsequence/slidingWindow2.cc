#include <Sequence/Fasta.hpp>
#include <Sequence/Alignment.hpp>
#include <Sequence/PolySites.hpp>
#include <Sequence/PolySNP.hpp>
#include <Sequence/PolyTableSlice.hpp>
#include <iostream>


//these are made explicit for example purposes
using std::vector;
using std::cout;
using std::endl;

//run a non-overlapping 100bp window over a SNP data set

int main(int argc, char **argv)
{
  const char *infilename = argv[1];

  vector<Sequence::Fasta> data;
  
  Sequence::Alignment::GetData(data,infilename);

  if ( Sequence::Alignment::IsAlignment(data) &&
       Sequence::Alignment::validForPolyAnalysis(data.begin(),data.end()) )
    {
      const unsigned alignmentLength = data[0].length();

      Sequence::PolySites SNPtable(data);

      Sequence::PolySNP analyzeRegion(&SNPtable);

      cout << "Tajima's D for the whole dataset is: "
	   << analyzeRegion.TajimasD()
	   << endl;

      Sequence::PolyTableSlice<Sequence::PolySites> windows(SNPtable.sbegin(),
							    SNPtable.send(),
							    100, //window length (bp)
							    100, //step size (bp)
							    alignmentLength);
      
      for(unsigned i = 0 ; i < windows.size() ; ++i)
	{
	  Sequence::PolySites window(windows[i]); //use copy constructor
	  Sequence::PolySNP analyzeWindow(&window);
	  cout << "Tajima's D for window "
	       << i
	       << " is: "
	       << analyzeWindow.TajimasD() 
	       << endl;
	}
    }
}
