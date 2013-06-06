#include <iostream>
#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TSlider.h>

#include "PadWithSlider.hpp"

using namespace std;

PadWithSlider::PadWithSlider()
{
  counter = 0;
}


void PadWithSlider::Show()
{
  TRandom3* r = new TRandom3();
  h = new TH1F("h","h",500,-10,10);
  for (int p=0; p<1000; p++) {
    r->SetSeed();
    h->Fill(r->Uniform(-5,5));
  }

  c1 = new TCanvas("c1");
  pad = new TPad("pad","lego pad",0.1,0.1,0.98,0.98);
  pad->SetFillColor(33);
  pad->Draw();
  pad->cd();
  h->SetFillColor(46);
  h->Draw();
  c1->cd();
  xslider = new TSlider("xslider","x",0.1,0.02,0.98,0.08);
  xslider->SetRange(0,0.1);
  cout <<"slider created" << endl;
  xslider->SetObject(this);
  cout <<"object set" << endl;
}                         




// Overriding the ExecuteEvent function of TPad.
void PadWithSlider::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  cout << "Here! " << counter << endl;
  float xmin = -10 + 10*xslider->GetMinimum();
  float xmax = h->GetXaxis()->GetXmax();
  h->SetAxisRange(xmin,xmax);
  pad->cd();
  h->Draw();
  c1->Update();
  counter++;
}
