// 1-D histogram drawing options
#include "TInterpreter.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TMath.h"
   
/* customizable code to read data from file and plot in graph*/
   
void graphs6_1_2()
{ 
   const char* input = "Ergebnisse6_1_2_Gar_MAJ.txt";
   //const Int_t xmax = 9700;
   const Int_t xmax = 20;
   //const Int_t mymax = xmax;
   //const Int_t mymin = 478866;
   const TString titlex = "Majority based size of length vectors wrt. Garside length";
   const TString titley = "Occurence";
	
   TCanvas *c1 = new TCanvas("c1","c1");
   c1->SetGrid();
   
   //TH1I *graph = new TH1I("graph", "", (mymax-mymin), mymin, mymax);
   TH1I *graph = new TH1I("graph", "", xmax+1, 0, xmax+1);

   //Fill with data
   Int_t x,y,ymin,ymax,xmin,mxmax,bin;
   Int_t nlines = 0;
   ymin = -1;
   ymax = 0;
   xmin = -1;
   mxmax = 0;
   
   ifstream in;
   in.open(input);

   while (1) {
		in >>"size" >> x >> ":" >> y;
		if (!in.good()){
			break;
		}
		if (ymin == -1){
			ymin = y;
			xmin = x;
		}
		if (y<ymin){
			ymin = y;
		}
		if (y>ymax){
			ymax = y;
		}
		if (x<xmin){
			xmin = x;
		}
		if (x>mxmax){
			mxmax = x;
		}
		graph->Fill(x,y);
   }
   
   printf("xmin=%i, xmax=%i", xmin, mxmax);
   printf("xmin=%i, xmax=%i", xmin, mxmax);
   
   graph->LabelsOption("X","a");
		
   graph->SetStats(0);
   graph->SetFillColor(38);
   
   for (int i=1;i<=xmax+1;i++)
   {
	   printf("\nBin %d: %d",i-1,graph->GetBinContent(i));
   }
   
   graph->Draw();
   
   //Axis stuff
   TAxis *xaxis = graph->GetXaxis();
   TAxis *yaxis = graph->GetYaxis();
   xaxis->SetTitle(titlex);
   yaxis->SetTitle(titley);
   xaxis->CenterTitle();
   yaxis->CenterTitle();
   
   in.close();
}
