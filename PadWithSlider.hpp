class PadWithSlider: public TPad{
public:
  PadWithSlider();
  void Show();
  
  
private:
  void ExecuteEvent(Int_t event, Int_t px, Int_t py);
  TCanvas* c1;
  TPad *pad;
  TH1F* h;
  TSlider* xslider;
  int counter;

  ClassDef(PadWithSlider,1);

};
