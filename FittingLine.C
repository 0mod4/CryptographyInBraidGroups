#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

/*Line Fit*/

Double_t fitFunction(Double_t *x, Double_t *par) {
  return par[0] + 1*par[1]*x[0] ;
}

void FittingLine() {
	
   TCanvas *c1 = new TCanvas("c1","Fitting",10,10,700,500);
   c1->SetFillColor(0);
   c1->SetFrameFillColor(0);
   c1->SetGrid();
 
   const int nBins = 17;
   Double_t data1[nBins] 	= {0.066048, 0.071443, 0.074411, 0.081578, 0.091777, 0.100589, 0.102928, 0.121907, 0.137952, 0.135277, 0.142871, 0.156199, 0.14026, 0.190504, 0.193613, 0.167253, 0.171193};
   
   TH1F *hist1 = new TH1F("hist1","",nBins,2,18);
   hist1->SetMarkerStyle(20);
   hist1->SetMarkerSize(0.4);
   hist1->SetMarkerColor(4);
   hist1->SetStats(0);

   for(int i=0; i < nBins;  i++)
   {
		if (data1[i] != 0)
			hist1->SetBinContent(i+1,data1[i]);
   }

   //Fitting Function
   TF1 *fit = new TF1("fit",fitFunction,2,18,3);
   fit->SetNpx(500);
   fit->SetLineWidth(4);
   
   // Fitting
   TFitResultPtr fitres;
   Double_t par[17][3];
   Double_t err[17][3];
     
   fit->SetLineColor(4);
   hist1->Fit("fit","V+","p");
   fit->GetParameters(par[1]);
   err[1][0] = fit->GetParError(0);
   err[1][1] = fit->GetParError(1);
   err[1][2] = fit->GetParError(2);
   
printf("N=%d;%f;%f;%f\n", 1, par[1][0], par[1][1], par[1][2]);
}
