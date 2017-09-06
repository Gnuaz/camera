#include <time.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButton.h>
#include <TGSlider.h>
#include <TGLabel.h>
#include <TMath.h>
#include "TGButton.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLayout.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGTextEntry.h"
#include "TGTripleSlider.h"
using namespace std;



//runge-kutta method or euler method can be choosed x( )
//approximation of gravity                           (x)
//some interface of drag and drop.                   (x)
//2D pendulum                                        ( )
//approximation 2 no spin                            (x)
void Check(Int_t event,Int_t x,Int_t y,TObject *selected);
void InitialCondition();
TText gl1 (-.5,.25,"1) Click where you want");
TText gl2 (-.5,.15,"   to set the pendulum");
TText gl3 (-.5,0,"2) And Drag to give");
TText gl4 (-.5,-.1,"  some momentum");
TText xaxis (-.1,.9, "Y");
TText yaxis (.9,.1,"X");
TEllipse *el;
TArrow *ar1, *ar2, *ar3;
TArrow *axisx, *axisy;
double v,r,l,g,s,xpos,ypos,xvel,yvel,zpos,n=2;
double xpo, ypo,xpoa,ypoa;
int t;

  Int_t gid,lid, spinid, tid,rid ;
int axis = 1;
int broke = 0;

struct presets{
  int ID;
  float lat;
  float len;
  float gra;
  float spi;
  float tim;
  float rad;
};


  


class testframe : public TGMainFrame{


 private:

  //declaration for gui
  
  TTimer *timer = new TTimer("pendulum()",20);
  TGTextButton *start, *reset, *exit, *preset, *save, *clear;
  TRootEmbeddedCanvas *fCanvas;
  TGLayoutHints *fLcan;
  TGCheckButton *enable, *approx1, *approx2, *approx3;
  TGHSlider *gravity, *latitude, *spintime, *time, *length,*viewpoint, *radius,*viewangle;
  TGNumberEntry *ngravity, *nlatitude, *nspintime, *ntime, *nlength, *nviewpoint,*nradius,*nviewangle;
  TGLayoutHints *fly, *fly1, *fly2, *flys, *flyne,*flyl;
  TGVButtonGroup *prop, *initial;
  TGLabel *label0,*label1,*label2,*label3,*label4,*label5,*label6,*label7;
  TGLabel *preset0, *preset1, *preset2;
  
  TGHorizontalFrame *fcframe, *fhframe0, *fhframe1, *fhframe2,*gravity1,*gravity2, *latitude1, *latitude2, *spintime1, *spintime2, *time1, *time2, *length1, *length2,*fhview,*fhpbutton, *fhpreset1,*fhpreset2, *viewpoint1, *viewpoint2,*viewangle2,*viewangle1,*radius1,*radius2,*trace;
  TGVerticalFrame *fvframe1, *fvframe2,*fcheck;
  TGHorizontal3DLine *hsep0,*hsep1,*hsep2,*hsep3,*hsep4;
  TGVertical3DLine *separator;
  TGComboBox *fCombo, *fPreset;

  //for visualizing
  TH2D         *foucault2d   = new TH2D("f1","f1",250,-50,50,250,-50,50);
  TGeoManager  *geom         = new TGeoManager("pendulum","Foucault`s Pendulum");
  TGeoMaterial *matsph       = new TGeoMaterial("sphere"      ,55.845,26,7.87);
  TGeoMedium   *Sph          = new TGeoMedium("sphere" , 1, matsph);
  TGeoMaterial *mattube      = new TGeoMaterial("tube"      , .935, 0., 10000.);
  TGeoMedium   *Tube         = new TGeoMedium("tube" , 1, mattube);
  TGeoMaterial *matcone      = new TGeoMaterial("cone"      ,  .938,1.,10000.);
  TGeoMedium   *Cone         = new TGeoMedium("cone" , 1, matcone);
  TGeoMaterial *matEmpty     = new TGeoMaterial("Empty"      , 0, 0, 0);
  TGeoMedium   *Empty        = new TGeoMedium("cone" , 1, matcone);

  
  Double_t worldx = 20., worldy = 20., worldz = 20.;
  TGeoVolume   *top          = geom->MakeBox("WORLD", Empty, worldx, worldy, worldz);
  TGeoVolume   *pendulumbox  = geom->MakeBox("PENDULUM", Empty, 4,4,200);
  TGeoVolume   *box          = geom->MakeBox("box", Tube, .5,.5,3);

  Double_t PendulumRadius = 3;
  Double_t TubeRadius     = 0.1;
  Double_t ConeRadius     = 0.5;
  Double_t ConeHeight     = 1;
  //Double_t TubeHeight;
  
  TGeoVolume * sphere  = geom->MakeSphere("sphere",Sph,0,PendulumRadius);
  TGeoVolume * tube;
  TGeoVolume * cone    = geom->MakeCone("cone",Cone,4,0,0,ConeRadius,0);
  
  
  Bool_t sstart,isstart,senable,sspin,saprox1,saprox2;
  int vid=2;

  //Values for calculating.
  Double_t tl,tg,tr,ts,tt,trr;
  Double_t txpos,typos, txvel,tyvel,tzpos;
  Double_t ttheta, tphi,teta=0;

  //some presets
  presets set[20];

  //declaration for 2d pendulum
  TLine *line;
  TEllipse *pend, *pend1, *pend2;
  TArrow *garrow, *carrow, *varrow, *garrow1, *carrow1, *varrow1;
  TText *legend1;
  TText *legend2;
  TText *legend3;





  
 public:


  
  TCanvas *c1;
  testframe();
  virtual ~testframe();
  void lsdo  (Int_t lpos = 0);
  void rsdo  (Int_t rpos = 0);
  void gsdo  (Int_t gpos = 0);
  void ssdo  (Int_t spos = 0);
  void tsdo  (Int_t tpos = 0);
  void vsdo  (Int_t vpos = 0);
  void asdo  (Int_t apos = 0);
  void rasdo (Int_t rapos = 0);
  void rado  (Long_t raval );
  void ado   (Long_t aval  );
  void vdo   (Long_t vval  );
  void gdo   (Long_t gval  );
  void ldo   (Long_t lval  );
  void tdo   (Long_t tval  );
  void sdo   (Long_t sval  );
  void rdo   (Long_t rval  );
  void fstart();
  void freset();
  void fpreset();
  void fclear();
  // void fsave();
  void EnableCondition();
  void setspin();
  void rotate();
  void fPresetdo(Int_t id);
  void graphicdrawing();
  //function for calculating
  Double_t Getxpos()  {return txpos;};
  Double_t Getypos()  {return typos;};
  Double_t Getxvel()  {return txvel;};
  Double_t Getyvel()  {return tyvel;};
  Double_t Getzpos()  {return tzpos;};
  Double_t Geteta()   {return teta;};
  Bool_t   Getsspin() {return sspin;};
  Bool_t   Getaprox2(){return saprox2;};
  Double_t Getrad()   {return (nradius   ->GetNumberEntry()->GetNumber());};
  Double_t Getl()     {return (nlatitude ->GetNumberEntry()->GetNumber());};
  Double_t Getr()     {return (nlength   ->GetNumberEntry()->GetNumber());};
  Double_t Getg()     {return (ngravity  ->GetNumberEntry()->GetNumber());};
  Double_t Gett()     {return (ntime     ->GetNumberEntry()->GetNumber());};
  Double_t Gets()     {return (nspintime ->GetNumberEntry()->GetNumber());};
  Double_t Getv()     {return (nviewpoint->GetNumberEntry()->GetNumber());};
  
  const TGCheckButton *Getapprox1(){return approx1;};  
  const TGCheckButton *Getspin()   {return enable;};
  TTimer              *GetTimer()  {return timer;};
  
  void Setxpos  (Double_t xpval) {txpos=xpval;};
  void Setypos  (Double_t ypval) {typos=ypval;};
  void Setxvel  (Double_t xval)  {txvel=xval;};
  void Setyvel  (Double_t yval)  {tyvel=yval;};
  void Setzpos  (Double_t zval)  {tzpos=zval;};
  void Settheta (Double_t tval)  {ttheta=tval;};
  void Setphi   (Double_t pval)  {tphi=pval;};
  void Seteta   (Double_t eval)  {teta=eval;};

  //function for tracing
  void graphdrawing();
  void geomsetting();
  void nospin();
  //apply();
  void InitialCondition();

  //void Check(Int_t event,Int_t x,Int_t y,TObject *selected);
  //void SetCondition();
  
  
  // void pendulum_broken();
};


testframe::testframe() : TGMainFrame(gClient->GetRoot(), 100,100,kHorizontalFrame)
{
  DontCallClose();
  SetCleanup(kDeepCleanup);
  geom->SetTopVolume(top);
  tube = geom->MakeBox("tube", Tube,0.1,0.1 , 65/2);
  pendulumbox-> AddNode(sphere,0 , new TGeoTranslation(0,0,0));
  pendulumbox-> AddNode(tube,1 , new TGeoTranslation(0,0,65));
  pendulumbox-> AddNode(cone,2 , new TGeoTranslation(0,0,-3));
  
  //foucault2d->SetPoint(0,20,20,1);
  //foucault2d->SetPoint(1,-20,-20,1);
  foucault2d->SetNameTitle("Tracing","Foucault`s Pendulum");
  foucault2d->SetBinContent(0,0,100);
  TPaveStats *ptstats = new TPaveStats(0.1,0.2,0.1,0.2,"brNDC");
  //TPaletteAxis *palette = new TPaletteAxis(0.111,0.213,0.123,0.2323,foucault2d);
  gStyle->SetOptStat("n");
  foucault2d->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(foucault2d);
  txpos=0;
  typos=-20;
  txvel=3;
  tyvel=0;
  
  //Create canvas and frame for canvas and add to the main frame, and 15 pixel margins all around 1 = 3 pixel.
  fcframe = new TGHorizontalFrame(this,705,705);
  fCanvas = new TRootEmbeddedCanvas("Canvas",fcframe,800,800);
  fLcan = new TGLayoutHints(kLHintsCenterX | kLHintsCenterY, 5,5,5,5);
  AddFrame(fcframe, fLcan);
  fcframe->AddFrame(fCanvas, fLcan);

  
  //setting large vertical frame and inner horizontal basic frames and seperator and layouthints.

  fly = new TGLayoutHints(kLHintsRight | kLHintsExpandY,0,5,0,0);
  
  fvframe1 = new TGVerticalFrame(this,250,700,kFixedWidth);
  fcheck = new TGVerticalFrame(fvframe1,250,700,kFixedWidth);
  AddFrame(fvframe1, fly);
  
  separator = new TGVertical3DLine(this);
  AddFrame(separator, fly);
  

  
  //and positioning every frames.
  fly1= new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,0,2);
   fly2= new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,0,5);
  fhframe0 = new TGHorizontalFrame(fvframe1);
  fhframe1 = new TGHorizontalFrame(fvframe1);
  fhframe2 = new TGHorizontalFrame(fvframe1);
  
  gravity1   = new TGHorizontalFrame(fvframe1);
  gravity2   = new TGHorizontalFrame(fvframe1);
  latitude1  = new TGHorizontalFrame(fvframe1);
  latitude2  = new TGHorizontalFrame(fvframe1);
  radius1    = new TGHorizontalFrame(fvframe1);
  radius2    = new TGHorizontalFrame(fvframe1);
  spintime1  = new TGHorizontalFrame(fvframe1);
  spintime2  = new TGHorizontalFrame(fvframe1);
  time1      = new TGHorizontalFrame(fvframe1);
  time2      = new TGHorizontalFrame(fvframe1);
  length1    = new TGHorizontalFrame(fvframe1);
  length2    = new TGHorizontalFrame(fvframe1);
  viewpoint1 = new TGHorizontalFrame(fvframe1);
  viewpoint2 = new TGHorizontalFrame(fvframe1);
  viewangle1 = new TGHorizontalFrame(fvframe1);
  viewangle2 = new TGHorizontalFrame(fvframe1);
  fhview     = new TGHorizontalFrame(fvframe1);
  fhpbutton  = new TGHorizontalFrame(fvframe1);
  fhpreset1  = new TGHorizontalFrame(fvframe1);
  //fhpreset2  = new TGHorizontalFrame(fvframe1);
  trace      = new TGHorizontalFrame(fvframe1);
  
  hsep0      = new TGHorizontal3DLine(fvframe1);
  hsep1      = new TGHorizontal3DLine(fvframe1);
  hsep2      = new TGHorizontal3DLine(fvframe1);

  
  fvframe1->AddFrame(fhframe0,fly1);

  fvframe1->AddFrame(fhframe2,new TGLayoutHints(kLHintsBottom| kLHintsExpandX,0,0,0,5));
  fvframe1->AddFrame(fhframe1,new TGLayoutHints(kLHintsBottom| kLHintsExpandX,0,0,0,5));
  fvframe1->AddFrame(fcheck,new TGLayoutHints(kLHintsBottom| kLHintsExpandX,0,0,0,5));

  fvframe1->AddFrame(fhpreset1 ,new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,10,5));
  // fvframe1->AddFrame(fhpreset2 ,fly2);
  fvframe1->AddFrame(hsep2     ,new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,5,10));
  fvframe1->AddFrame(latitude1 ,fly1);
  fvframe1->AddFrame(latitude2 ,fly2);
  fvframe1->AddFrame(length1   ,fly1);
  fvframe1->AddFrame(length2   ,fly2);
  fvframe1->AddFrame(radius1   ,fly1);
  fvframe1->AddFrame(radius2   ,fly2);
  fvframe1->AddFrame(gravity1  ,fly1);
  fvframe1->AddFrame(gravity2  ,fly2);
  fvframe1->AddFrame(spintime1 ,fly1);
  fvframe1->AddFrame(spintime2 ,fly2);
  fvframe1->AddFrame(time1     ,fly1);
  fvframe1->AddFrame(time2     ,fly2);
  fvframe1->AddFrame(fhpbutton ,fly2);
  fvframe1->AddFrame(hsep0     ,new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,5,10));
  fvframe1->AddFrame(fhview    ,fly2);
  fvframe1->AddFrame(trace     ,fly2);
  fvframe1->AddFrame(viewpoint1,fly1);
  fvframe1->AddFrame(viewpoint2,fly2);
  fvframe1->AddFrame(viewangle1,fly1);
  fvframe1->AddFrame(viewangle2,fly2);
  fvframe1->AddFrame(hsep1     ,new TGLayoutHints(kLHintsTop | kLHintsExpandX,0,0,5,10));
  
  

  
  //setting buttons.
  start     = new TGTextButton(fhframe2   , "&Start");
  reset     = new TGTextButton(fhframe2   , "&Reset");
  exit      = new TGTextButton(fhframe2   , "&Exit","gApplication->Terminate(0)");
  preset    = new TGTextButton(fhpbutton  , "Reset &Properties");
  //save      = new TGTextButton(fhpreset2  , "Sa&ve");
  clear     = new TGTextButton(trace   , "&Clear Trace");
  
  sstart    = kFALSE;
  isstart=kFALSE;
  
  reset     ->SetToolTipText("Reset the pendulum");
  preset    ->SetToolTipText("Reset all properties");
  exit      ->SetToolTipText("Quit this Window");
  //save      ->SetToolTipText("Save all properties of now");
  clear     ->SetToolTipText("Clear tracing of pendulum");

  //fhpreset2  ->AddFrame(save    , new TGLayoutHints(kLHintsRight  ,2,2,5,0));
  fhpbutton  ->AddFrame(preset  , new TGLayoutHints(kLHintsTop    | kLHintsExpandX ,0,0,10,5));
  fhframe2   ->AddFrame(start   , new TGLayoutHints(kLHintsTop    | kLHintsExpandX ,2,2,2,2));
  fhframe2   ->AddFrame(reset   , new TGLayoutHints(kLHintsTop    | kLHintsExpandX ,2,2,2,2));
  fhframe2   ->AddFrame(exit    , new TGLayoutHints(kLHintsTop    | kLHintsExpandX ,2,2,2,2));
  trace   ->AddFrame(clear    , new TGLayoutHints(kLHintsTop    | kLHintsExpandX ,2,2,2,2));
  //save       ->Resize(100,20);
 
  
  //checkbutton
  enable = new TGCheckButton(fcheck, new TGHotString("Spinning coordinate/\nStatic coordinate"));
  enable->SetToolTipText("State of the coordinate");
  fcheck->AddFrame(enable, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,10));
  //enable->SetOn(kFALSE);
  //enable->Connect("Clicked()","testframe",this,"setspin()");
  sspin=kFALSE;

  approx1 = new TGCheckButton(fcheck, new TGHotString("Approximation #1 : \n   Spinning of planet effect only Colioli`s effect"));
  approx2 = new TGCheckButton(fcheck, new TGHotString("Approximation #2 : \n   There is no spin of planet"));
  fcheck->AddFrame(approx1, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5));

  fcheck->AddFrame(approx2, new TGLayoutHints(kLHintsTop|kLHintsLeft,0,0,5));
  approx2->Connect("Clicked()","testframe",this,"nospin()");
  saprox1=kFALSE;
  saprox2=kFALSE;

  
  //slider and number entry and label.
  
   //set LayoutHints
   flyne = new TGLayoutHints(kLHintsTop|kLHintsRight,0,10);
   flys = new TGLayoutHints(kLHintsTop|kLHintsCenterX|kLHintsExpandX);
   flyl = new TGLayoutHints(kLHintsBottom|kLHintsLeft,8);


  //latitude
  label0 = new TGLabel(latitude1, "Latitude(Angle)");
  latitude1->AddFrame (label0,flyl);
  nlatitude = new TGNumberEntry(latitude1, 37.5,5,999,TGNumberFormat::kNESRealOne,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,
			       -180,180);
  latitude1->AddFrame(nlatitude,flyne);
 
 
  latitude = new TGHSlider(latitude2,100,kSlider1|kScaleBoth);
  latitude->SetScale(15);
  latitude->SetRange(-180,180);
  latitude->SetPosition(75);
  lid = latitude->WidgetId();
  latitude->Connect("PositionChanged(Int_t)", "testframe", this, "lsdo(Int_t)");
  latitude2->AddFrame(latitude, flys);
  //nlatitude->Connect("ValueChanged(Long_t)","testframe",this,"ldo(Long_t)");
  nlatitude->Connect("ValueSet(Long_t)","testframe",this,"ldo(Long_t)");
  //(nlatitude->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "ldo(Long_t)");
  tl= (nlatitude ->GetNumberEntry()->GetNumber());
  
  //length
  label1  = new TGLabel(length1, "Length(m)");
  length1 ->AddFrame (label1,flyl);
  nlength = new TGNumberEntry(length1, 65,5,999,TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
			       20, 99999);
  length1 ->AddFrame(nlength,flyne);
 
  length  =new TGHSlider(length2,100,kSlider1|kScaleBoth);
  rid     =length->WidgetId();
  length  ->Connect("PositionChanged(Int_t)", "testframe", this, "rsdo(Int_t)");
  length  ->SetScale(20);
  length  ->SetRange(20,200);
  length  ->SetPosition(65);
  length2 ->AddFrame(length, flys);
  //(nlength->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "rdo(Long_t)");
  //nlength ->Connect("ValueChanged(Long_t)","testframe",this,"rdo(Long_t)");
  nlength ->Connect("ValueSet(Long_t)","testframe",this,"rdo(Long_t)");
  tr =(nlength   ->GetNumberEntry()->GetNumber());

  //radius
  label7 = new TGLabel(radius1, "Radius of the planet(Km)");
  radius1->AddFrame (label7,flyl);
  nradius = new TGNumberEntry(radius1, 6400,5,999,TGNumberFormat::kNESInteger,
				                       TGNumberFormat::kNEANonNegative,
				TGNumberFormat::kNELLimitMinMax,0,99999);
  radius1->AddFrame(nradius,flyne);
  trr = (nradius   ->GetNumberEntry()->GetNumber());
  radius = new TGHSlider(radius2,10000,kSlider1|kScaleBoth);
  //(nradius->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "rado(Long_t)");
  radius->Connect("PositionChanged(Int_t)", "testframe", this, "rasdo(Int_t)");
  // nradius->Connect("ValueChanged(Long_t)","testframe",this,"sdo(Long_t)");
  nradius->Connect("ValueSet(Long_t)","testframe",this,"rado(Long_t)");
  
  radius->SetRange(1,99999);
  radius->SetPosition(6400);
  radius2->AddFrame(radius, flys);
  
  //gravity
  label2   =  new TGLabel(gravity1, "Gravity(m/s^2)");
  gravity1 -> AddFrame (label2,flyl);
  ngravity =  new TGNumberEntry(gravity1, 9.8,5,999,TGNumberFormat::kNESRealOne,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
			       1, 99999);
  
  gravity1 -> AddFrame(ngravity,flyne);
  
  gravity  =  new TGHSlider(gravity2,100, kSlider1 | kScaleBoth );
  gid      =  gravity->WidgetId();
  gravity  -> Connect("PositionChanged(Int_t)", "testframe",
                     this, "gsdo(Int_t)");
  (ngravity->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "gdo(Long_t)");
  //ngravity -> Connect("ValueChanged(Long_t)","testframe",this,"gdo(Long_t)");
  ngravity -> Connect("ValueSet(Long_t)","testframe",this,"gdo(Long_t)");
  gravity  -> SetRange(1,200);
  gravity  -> SetPosition(10);
  gravity2 -> AddFrame(gravity, flys);
  tg= (ngravity  ->GetNumberEntry()->GetNumber());
  
  //spintime
  label3 = new TGLabel(spintime1, "Spinning time of earth(s)");
  spintime1->AddFrame (label3,flyl);
  nspintime = new TGNumberEntry(spintime1, 86400,7,999,TGNumberFormat::kNESInteger,
				                       TGNumberFormat::kNEANonNegative,
				                       TGNumberFormat::kNELLimitMinMax,0, 9999999);
  spintime1->AddFrame(nspintime,flyne);
 
  spintime = new TGHSlider(spintime2,1000,kSlider1|kScaleBoth);
  // (nspintime->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "sdo(Long_t)");
  spintime->Connect("PositionChanged(Int_t)", "testframe", this, "ssdo(Int_t)");
  // nspintime->Connect("ValueChanged(Long_t)","testframe",this,"sdo(Long_t)");
  nspintime->Connect("ValueSet(Long_t)","testframe",this,"sdo(Long_t)");
  spinid = spintime->WidgetId();
  spintime->SetRange(100,99999);
  spintime->SetPosition(86400);
  spintime2->AddFrame(spintime, flys);
  ts =(nspintime ->GetNumberEntry()->GetNumber());
  
  //time
   label4 = new TGLabel(time1, "Time Multiplier(*n)");
  time1->AddFrame (label4,flyl);
  ntime = new TGNumberEntry(time1, 1,5,999,TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
			       0, 99999);
  time1->AddFrame(ntime,flyne);
 
  time = new TGHSlider(time2,100,kSlider1|kScaleBoth);
  time->Connect("PositionChanged(Int_t)", "testframe", this, "tsdo(Int_t)");
  //(ntime->GetNumberEntry())->Connect("ReturnPressed()","testframe", this, "tdo(Long_t)");
  //ntime->Connect("ValueChanged(Long_t)","testframe",this,"tdo(Long_t)");
  ntime->Connect("ValueSet(Long_t)","testframe",this,"tdo(Long_t)");
  tid = time->WidgetId();
  time->SetRange(1,100);
  time->SetPosition(1);
  time2->AddFrame(time, flys);
  tt=(ntime     ->GetNumberEntry()->GetNumber());
  
  //viewer selecter sector
  label5 = new TGLabel(fhview, "Viewer                        : ");
  fCombo = new TGComboBox(fhview, -1,
                  kHorizontalFrame | kSunkenFrame | kDoubleBorder,
                  GetWhitePixel());
 
  fCombo  ->AddEntry("3D Model",1);
  fCombo  ->AddEntry("Tracing",2);
  //fCombo ->AddEntry("2D Graphics",3);
  fCombo  ->Select(2);
  fCombo  ->Resize(100,20);
  fhview  ->AddFrame(label5 , new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsCenterY));
  fhview  ->AddFrame(fCombo , new TGLayoutHints(kLHintsRight|kLHintsExpandX));

  //viewpoint setting
  
  label6  = new TGLabel(viewpoint1, "Viewpoint(m) (Only for 3Dmodel)");
  viewpoint1 ->AddFrame (label6,new TGLayoutHints(kLHintsBottom|kLHintsLeft));
  nviewpoint = new TGNumberEntry(viewpoint1, 20,5,999,TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
			       10, 2000);
  viewpoint1 ->AddFrame(nviewpoint,new TGLayoutHints(kLHintsBottom|kLHintsRight));
 
  viewpoint  =new TGHSlider(viewpoint2,100,kSlider1|kScaleBoth);
  viewpoint  ->Connect("PositionChanged(Int_t)", "testframe", this, "vsdo(Int_t)");
  viewpoint  ->SetScale(20);
  viewpoint  ->SetRange(10,200);
  viewpoint  ->SetPosition(20);
  viewpoint2 ->AddFrame(viewpoint,flys);
  // nviewpoint ->Connect("ValueChanged(Long_t)","testframe",this,"vdo(Long_t)");
  // nviewpoint ->Connect("ValueSet(Long_t)","testframe",this,"vdo(Long_t)");
  
//viewangle setting
  label7     = new TGLabel(viewangle1, "Viewangle(Angle)");
  viewangle1 ->AddFrame (label6,new TGLayoutHints(kLHintsBottom|kLHintsLeft));
  nviewangle = new TGNumberEntry(viewangle1,37.5,5,999,TGNumberFormat::kNESRealOne,TGNumberFormat::kNEAAnyNumber,TGNumberFormat::kNELLimitMinMax,
			       -180,180);
  viewangle1 ->AddFrame(nviewangle,new TGLayoutHints(kLHintsBottom|kLHintsRight));
 
  viewangle  =new TGHSlider(viewangle2,100,kSlider1|kScaleBoth);
  viewangle  ->Connect("PositionChanged(Int_t)", "testframe", this, "asdo(Int_t)");
  viewangle  ->SetScale(20);
  viewangle  ->SetRange(-180,180);
  viewangle  ->SetPosition(75);
  viewangle2 ->AddFrame(viewangle,flys);
  
  //preset selecter sector
  
  preset0    = new TGLabel   (fhpreset1  , "Preset                        : ");
  fPreset    = new TGComboBox(fhpreset1  , -2,
                  kHorizontalFrame | kSunkenFrame | kDoubleBorder,
                  GetWhitePixel());
  fPreset    ->AddEntry("Seoul, Earth"      ,1);
  fPreset    ->AddEntry("North Pole, Earth" ,2);
  fPreset    ->AddEntry("South Pole, Earth" ,3);
  fPreset    ->AddEntry("North Pole, Moon"  ,4);
  fPreset    ->AddEntry("South Pole, Moon"  ,5);
  fPreset    ->AddEntry("North Pole, Mars"  ,6);
  fPreset    ->AddEntry("South Pole, Mars"  ,7);
  fPreset    ->Select(1);
  fPreset    ->Resize(100,20);
  fhpreset1  ->AddFrame(preset0 , new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsCenterY));
  fhpreset1  ->AddFrame(fPreset , new TGLayoutHints(kLHintsRight|kLHintsExpandX));
  fPreset    ->Connect("Selected(Int_t)","testframe",this,"fPresetdo(Int_t)");
  
  //seoul
  set[0].ID      = 1;
  set[0].lat     = 37.5;
  set[0].len     = 65;
  set[0].gra     = 9.8;
  set[0].spi     = 86400;
  set[0].tim     = 1;
  set[0].rad     = 6400;
  //northpole
  set[1].ID      = 2;
  set[1].lat     = 90;
  set[1].len     = 65;
  set[1].gra     = 9.8;
  set[1].spi     = 86400;
  set[1].tim     = 1;
  set[1].rad     = 6400;
  //southpole
  set[2].ID      = 3;
  set[2].lat     = -90;
  set[2].len     = 65;
  set[2].gra     = 9.8;
  set[2].spi     = 86400;
  set[2].tim     = 1;
  set[2].rad     = 6400;
  //n moon
  set[3].ID      = 4;
  set[3].lat     = 90;
  set[3].len     = 65;
  set[3].gra     = 1.6;
  set[3].spi     = 2360534;
  set[3].tim     = 1;
  set[3].rad     = 1735;
  //s moon
  set[4].ID      = 5;
  set[4].lat     = -90;
  set[4].len     = 65;
  set[4].gra     = 1.6;
  set[4].spi     = 2360534;
  set[4].tim     = 1;
  set[4].rad     = 1735;
  //n mars
  set[5].ID      = 6;
  set[5].lat     = 90;
  set[5].len     = 65;
  set[5].gra     = 3.7;
  set[5].spi     = 88642;
  set[5].tim     = 1;
  set[5].rad     = 3390;
  //s mars
  set[6].ID      = 7;
  set[6].lat     = -90;
  set[6].len     = 65;
  set[6].gra     = 3.7;
  set[6].spi     = 88642;
  set[6].tim     = 1;
  set[6].rad     = 3390;
  
  
  
  //button Signal setting
  start   ->Connect("Clicked()","testframe",this,"fstart()");
  reset   ->Connect("Clicked()","testframe",this,"freset()");
  preset  ->Connect("Clicked()","testframe",this,"fpreset()");
  clear   ->Connect("Clicked()","testframe",this,"fclear()");

  //fstart();
  SetWindowName("Foucault's Pendulum");
  MapSubwindows();
  Resize();
  MapRaised();
}

testframe::~testframe()
{
   // Clean up all widgets, frames and layouthints that were used
   Cleanup();
}


  //basical control function of sliders and numberentries
void testframe::lsdo(Int_t lpos){
  double_t npos = lpos;
  nlatitude->Clear();
  nlatitude->SetNumber(npos/2);
    tl= (nlatitude ->GetNumberEntry()->GetNumber());
  //cout<<pos<<endl;
}
void testframe::rsdo(Int_t rpos){
  double_t npos = rpos;
  nlength->Clear();
  nlength->SetNumber(npos);
  tr =(nlength   ->GetNumberEntry()->GetNumber());
  //cout<<pos<<endl;
}
void testframe::gsdo(Int_t gpos){
  double_t npos = gpos;
  ngravity->Clear();
  ngravity->SetNumber(npos);
  tg= (ngravity  ->GetNumberEntry()->GetNumber());
  //cout<<pos<<endl;
}
void testframe::ssdo(Int_t spos){
  double_t npos = spos;
  nspintime->Clear();
  nspintime->SetNumber(npos);
  ts =(nspintime ->GetNumberEntry()->GetNumber());
  //cout<<spos<<endl;
}
void testframe::tsdo(Int_t tpos){
  double_t npos = tpos;
  ntime->Clear();
  ntime->SetNumber(npos);
    tt=(ntime     ->GetNumberEntry()->GetNumber());
  //cout<<pos<<endl;
}
void testframe::vsdo(Int_t vpos){
  double_t npos = vpos;
  nviewpoint->Clear();
  nviewpoint->SetNumber(npos);

  
  // gPad->GetView()->SetRange(-v,-v,0,v,v,v);
  
  //cout<<pos<<endl;
}
void testframe::asdo(Int_t apos){
  double_t npos = apos;
  nviewangle->Clear();
  nviewangle->SetNumber(npos/2);
}
void testframe::rasdo(Int_t rapos){
  double_t npos = rapos;
  nradius->Clear();
  nradius->SetNumber(npos);
  trr=(nradius   ->GetNumberEntry()->GetNumber());
}
void testframe::rado(Long_t raval){
  radius ->SetPosition(nradius->GetNumberEntry()->GetNumber());
  trr=  (nradius   ->GetNumberEntry()->GetNumber());
}
void testframe::ado(Long_t aval){
  viewangle ->SetPosition(nviewangle->GetNumberEntry()->GetNumber()*2);
}
void testframe::vdo(Long_t vval){
  viewpoint ->SetPosition(nviewpoint->GetNumberEntry()->GetNumber());
}
void testframe::ldo(Long_t lval){
  latitude  ->SetPosition(nlatitude->GetNumberEntry()->GetNumber()*2);
  tl= (nlatitude ->GetNumberEntry()->GetNumber());
}
void testframe::rdo(Long_t rval){
  length    ->SetPosition(nlength->GetNumberEntry()->GetNumber());
  tr =(nlength   ->GetNumberEntry()->GetNumber());
}
void testframe::gdo(Long_t gval){
  gravity   ->SetPosition(ngravity->GetNumberEntry()->GetNumber());
  tg= (ngravity  ->GetNumberEntry()->GetNumber());
}
void testframe::sdo(Long_t sval){
  spintime  ->SetPosition(nspintime->GetNumberEntry()->GetNumber());
  ts =(nspintime ->GetNumberEntry()->GetNumber());
}
void testframe::tdo(Long_t tval){
  time      ->SetPosition(ntime->GetNumberEntry()->GetNumber());
  tt=(ntime     ->GetNumberEntry()->GetNumber());
}



//definition of functions for buttons
void testframe::fstart(){
  if (!isstart)
    {
      InitialCondition();
      start->SetText("&Pause");
      start->SetDown(kFALSE);
      sstart=kTRUE;
      isstart=kTRUE;
      senable=kFALSE;
    }
  else{
    if (!sstart)
      {
	start->SetText("&Pause");
	start->SetDown(kFALSE);
	sstart=kTRUE;
	isstart=kTRUE;
	senable=kFALSE;
	//EnableCondition();
	//pendulum();
	if(!gPad){teta = gPad->GetView()->GetLongitude();}
	timer->SetCommand("pendulum()");
	cout<<"pendulum is working."<<endl;
	timer->TurnOn();
      
      }
    else
      {
	start->SetText("&Resume");
	start->SetDown(kTRUE);
	sstart=kFALSE;
	start->SetToolTipText("You have to pause for changing properties");
	senable=kTRUE;
	cout<<"pendulum is paused"<<endl;
	//EnableCondition();
	timer->Stop();
      }
  }
}
void testframe::fclear(){
  
  foucault2d->Reset();



}
void testframe::freset(){
  if (isstart)
    {start->SetText("&Start");
      start->SetDown(kFALSE);
      sstart=kFALSE;
      start->SetToolTipText("");
      isstart = kFALSE;
      
    }
  //senable=kTRUE;
  //EnableCondition();
 
  //sspin=kFALSE;
  timer->Stop();
  foucault2d->Reset();
  gStyle->SetOptStat("n");
  //setspin();
  n = 2;
  txpos =   0;
  typos = -20;
  txvel =   2;
  tyvel =   0;
  fCanvas->GetCanvas()->Clear();
   fCanvas->GetCanvas()->Update();
  
  cout<<"pendulum is reseted!"<<endl;
  if (broke)
    {
      fpreset();
      broke = 0;
    }
}
void testframe::fpreset(){
  int pid =  fPreset->GetSelected()-1;
  //senable=kTRUE;
  //EnableCondition();
  nlatitude  ->SetNumber(set[pid].lat);
             ldo(2*set[pid].lat);
	     nlength    ->SetNumber(set[pid].len);
             rdo(set[pid].len);
  ngravity   ->SetNumber(set[pid].gra);
  gravity    ->SetPosition(ngravity->GetNumberEntry()->GetNumber());
  nspintime  ->SetNumber(set[pid].spi);
             sdo(set[pid].spi);
	     ntime      ->SetNumber(set[pid].tim);
             tdo(set[pid].tim);
	     //sspin=kTRUE;
	     //setspin();
  //fCanvas->GetCanvas()->Clear();
  //fCanvas->GetCanvas()->Update();
  
  cout<<"properties are resetted!"<<endl;
  
}
//void testframe::fsave(){
//  int pid =  fPreset->GetNumberOfEntries()-1;
  
  
//  cout<<"Properties are Saved!"<<endl;
  
//}
void testframe::fPresetdo(Int_t id){
  int pid = id -1;
  
  //senable=kTRUE;
  //EnableCondition();
  nlatitude  ->SetNumber(set[pid].lat);
             ldo(2*set[pid].lat);
	     //nlength    ->SetNumber(set[pid].len);
             //rdo(set[pid].len);
  ngravity   ->SetNumber(set[pid].gra);
  gravity    ->SetPosition(ngravity->GetNumberEntry()->GetNumber());
  nspintime  ->SetNumber(set[pid].spi);
             sdo(set[pid].spi);
	     nradius    ->SetNumber(set[pid].rad);
	     rado(set[pid].rad);
	     //ntime      ->SetNumber(set[pid].tim);
             //tdo(set[pid].tim);
  sspin=kTRUE;
  setspin();
  //fCanvas->GetCanvas()->Clear();
  //fCanvas->GetCanvas()->Update();
  
  cout<<"properties are resetted!"<<endl;
  
}


void testframe::EnableCondition(){
  nlatitude ->SetState(senable);
  latitude  ->SetState(senable);
  nlength   ->SetState(senable);
  length    ->SetState(senable);
  ngravity  ->SetState(senable);
  gravity   ->SetState(senable);
  nspintime ->SetState(senable);
  spintime  ->SetState(senable);
  ntime     ->SetState(senable);
  time      ->SetState(senable);
  if (senable){enable->SetState(kButtonEngaged);}
  else {enable->SetState(kButtonDisabled);}
  
}
void testframe::setspin(){
  if (sspin)
    {enable->SetOn(kFALSE);
      if(!gPad){
      if (fCombo->GetSelected()==2)
	{teta = gPad->GetView()->GetLongitude();}
      }
      sspin=kFALSE;
    }
  else
    {enable->SetOn();
      sspin=kTRUE;
    }
}
void testframe::nospin(){
  if (!saprox2)
    {
      nspintime ->SetState(saprox2);
      spintime  ->SetState(saprox2);
      saprox2=kTRUE;
    }
  else
    {
      nspintime ->SetState(saprox2);
      spintime  ->SetState(saprox2);
      saprox2=kFALSE;
    }
}

void testframe::rotate(){
 
    {
      Int_t irep;
      TView *view = gPad->GetView();
      gPad->GetView()->SetView(teta,(-(nviewangle->GetNumberEntry()->GetNumber())+90),view->GetPsi(),irep);
    }
  
}


// functions for 3d modeling

void testframe::graphdrawing(){
  vid = fCombo->GetSelected();
  foucault2d->Fill(xpos,ypos,1);
  if(vid==1){
    timer->SetTime(40);
   geomsetting();
   
  }
  if(vid==2){
  //foucault2d->Reset();
    timer->SetTime(0);
  gStyle->SetOptStat("n");
  foucault2d->Draw("colz2");
   
  fCanvas->GetCanvas()->Modified();
    fCanvas->GetCanvas()->Update();
  }
  if(vid==3){
    timer->SetTime(40);
    graphicdrawing();
  }
}

void testframe::graphicdrawing(){
  fCanvas->Clear();
  //TLine *line;
  //TEllipse *pend, *pend1, *pend2;
  //TArrow *garrow, *carrow, *varrow, *garrow1, *carrow1, *varrow1;
  //TText *legend1;
  //TText *legend2;
  // TText *legend3;
  v = nviewpoint->GetNumberEntry()->GetNumber();
  //fCanvas->Range(-v,v);
  
}



  

void testframe::geomsetting()
{
  //top->RemoveNode(pendulum);
  top->ClearNodes();
  v = nviewpoint->GetNumberEntry()->GetNumber();
  for(int i = 1; i<9; i++){
  double ttt= (45*i)*M_PI/180;
  if (i == 4)
    {box->SetLineColor(1);
      top->AddNode(box,i,new TGeoTranslation((v/2+20)*sin(ttt),(v/2+20)*cos(ttt),0));
      
    }
  //box->SetLineColor(kBlack);
  box->SetLineColor(1);
  top->AddNode(box,i,new TGeoTranslation((v/2+20)*sin(ttt),(v/2+20)*cos(ttt),0));
  }
  top->AddNode(box,10,new TGeoTranslation(0,0,0));
  
    r =  nlength->GetNumberEntry()->GetNumber();
    
  tube->ClearShape();
  tube = geom->MakeBox("tube", Tube, 0.1,0.1, r/2);
  pendulumbox-> ClearNodes();
  pendulumbox-> AddNode(sphere,0 , new TGeoTranslation(0,0,0));
  pendulumbox-> AddNode(tube,1 , new TGeoTranslation(0,0,r/2));
  pendulumbox-> AddNode(cone,2 , new TGeoTranslation(0,0,-3));
  TGeoTranslation t1(xpos,ypos,zpos);
  TGeoRotation r1("r1",tphi,ttheta,0);
  TGeoCombiTrans *c1 = new TGeoCombiTrans(t1,r1);
  
			     
  top->AddNodeOverlap(pendulumbox, 0, c1);
  //geom->SetVisLevel(2);
  top->Draw();
  
  
  gPad->GetView()->SetRange(-v,-v,0,v,v,v);
  rotate();
  gPad->Modified();
  gPad->Update();
  
  //fCanvas->GetCanvas()->Modified();
  //fCanvas->GetCanvas()->Update();
  
  
}
//Initial Condition Setting

void testframe::InitialCondition(){
  c1 = new TCanvas("c1","Foucault`s Pendulum",800,800);
  c1->Range(-1,-1,1,1);
  
  gl1.Draw();
  gl2.Draw();
  gl3.Draw();
  gl4.Draw();
  axisx = new TArrow(-1,0,1,0,.05,"|>");
      axisy = new TArrow(0,-1,0,1,.05,"|>");
      axisx->Draw();
      axisy->Draw();
      xaxis.Draw();
      yaxis.Draw();
  
  gPad->Modified();
  c1->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject *)", 0 ,0,
      "Check(Int_t,Int_t,Int_t,TObject *");
  
    }


testframe *fp;





//below functions is calculating functions.

void pendulum(){
  
  //  r =     nlength->GetNumberEntry()->GetNumber();
  //  l =   nlatitude->GetNumberEntry()->GetNumber();
  //  g =    ngravity->GetNumberEntry()->GetNumber();
  //  s =   nspintime->GetNumberEntry()->GetNumber();
  //  t =       ntime->GetNumberEntry()->GetNumber();
  r = fp->Getr();
  l = fp->Getl();
  g = fp->Getg();
  s = fp->Gets();
  t = fp->Gett();
  Double_t ra = fp->Getrad();
  xpos = fp->Getxpos();
  ypos = fp->Getypos();
  xvel = fp->Getxvel();
  yvel = fp->Getyvel();
  double eta = fp->Geteta();
  float phi = M_PI*l/180;
  Bool_t sspin = fp->Getsspin();
  //Bool_t sgrav = fp->Getsgrav();
  float omega = 2*M_PI/s;
  if (fp->Getaprox2()){omega = 0;}
  float timem = 0.005;
  double grav;
  if (!fp->Getapprox1()->GetState()){grav= g - ra*1000*omega*omega*cos(phi)*cos(phi);}
 else {grav =g;}
  if(!fp->Getaprox2()){
    if(fp->Getspin()->GetState()){
      eta = eta + 10*t*timem/s*360;
  
    }

  }
    for(int i = 10*t; i>0; i--){
      double del  = sqrt(xpos*xpos+ypos*ypos);
      double sinu = del/r;
      double acc  = grav*fabs(sinu)       ;//g(e)
      double xacc = -acc*xpos/del+ 2*omega*yvel*sin(phi);
      double yacc;
      //cout << fp->Getapprox1()->GetState()<<grav<<"    "<<g<<endl;
      if (fp->Getapprox1()->GetState()){
	yacc = -acc*ypos/del- 2*omega*xvel*sin(phi);
      }
      else{
	yacc = -acc*ypos/del- 2*omega*xvel*sin(phi)- ra*1000*omega*omega*cos(phi);
      }
      xvel = xvel + xacc*timem;
      yvel = yvel + yacc*timem;
      xpos = xpos + xvel*timem;
      ypos = ypos + yvel*timem;
      zpos = r - sqrt(r*r-xpos*xpos - ypos*ypos);
       double ttheta = asin(xpos/r)*180/M_PI;
       double tphi = -atan(xpos/ypos)*180/M_PI;
       //  printf("%10f,%10f \n\n", ttheta,tphi);
      
        if (del*del > r*r)
	  {
	     cout<<"Pendulum has broken!"<<endl;
	     broke =1;
	   fp->freset();
	   break;
	 }
  }
    double del  = sqrt(xpos*xpos+ypos*ypos);
     double ttheta = asin(del/r)*180/M_PI;
     double tphi;
     if(xpos>=0&&ypos>=0)
       { tphi =  -atan(xpos/ypos)*180/M_PI;
       ttheta = asin(del/r)*180/M_PI;}
     if(xpos>=0&&ypos<0)
       { tphi =  -atan(xpos/ypos)*180/M_PI;
       ttheta = -asin(del/r)*180/M_PI;}
     if(xpos<0&&ypos>=0)
       { tphi =  -atan(xpos/ypos)*180/M_PI;
       ttheta = asin(del/r)*180/M_PI;}
     if(xpos<0&&ypos<0)
       { tphi =  -atan(xpos/ypos)*180/M_PI;
       ttheta = -asin(del/r)*180/M_PI;}
     fp->Seteta(eta);
     fp->Setxpos(xpos);
     fp->Setypos(ypos);
     fp->Setxvel(xvel);
     fp->Setyvel(yvel);
     fp->Settheta(ttheta);
     fp->Setphi(tphi);
       //fp->graphsetting(xpos,ypos,zpos);
       n = n+ 1;
       fp->graphdrawing();
       //       fp->drawcanvas();
}

void Check(Int_t event,Int_t x,Int_t y,TObject *selected)
{
  if (x<0){x=0;}
  if (x>800){x=800;}
  if (y<0){y=0;}
  if (y>800){y=800;}
    
    //cout << event<<endl;
  if (event == 51)
    {
      fp->c1->Clear();
      xpo=x;
      ypo=y;
     axisx = new TArrow(-1,0,1,0,.05,"|>");
      axisy = new TArrow(0,-1,0,1,.05,"|>");
      axisx->Draw();
      axisy->Draw();
      xaxis.Draw();
      yaxis.Draw();
      xpo = (xpo-400)/400;
      ypo = (ypo-400)/400;
      el = new TEllipse(xpo,-ypo,.1,.1);
      el->SetFillColor(kBlack);
      el->Draw();
      //printf("%10f %10f %10d\n %10d %10d \n",xpo,ypo,event,x,y);
      fp->c1->Update();
      
    }
  if (event == 1)
    {
       fp->c1->Clear();
      xpo=x;
      ypo=y;
      xpos =x;
      ypos =y;
      xpo = (xpo-400)/400;
      ypo = (ypo-400)/400;
      axisx = new TArrow(-1,0,1,0,.05,"|>");
      axisy = new TArrow(0,-1,0,1,.05,"|>");
      axisx->Draw();
      axisy->Draw();
      xaxis.Draw();
      yaxis.Draw();
      el = new TEllipse(xpo,-ypo,.1,.1);
      el->SetFillColor(kBlack);
      el->Draw();
      //printf("%10f %10f %10d\n %10d %10d \n",xpo,ypo,event,x,y);
      fp->c1->Update();
    }
  if (event == 21)
    {
      fp->c1->Clear();
      xpoa=x;
      ypoa=y;
      xpoa = (xpoa-400)/400;
      ypoa = (ypoa-400)/400;
      axisx = new TArrow(-1,0,1,0,.05,"|>");
      axisy = new TArrow(0,-1,0,1,.05,"|>");
      axisx->Draw();
      axisy->Draw();
      a1= new TArrow(xpo,-ypo,xpoa,-ypoa,.07,"|>");
      a1->SetLineWidth(10);
      a1->SetFillColor(kRed);
      a1->SetLineColor(kRed);
      el = new TEllipse(xpo,-ypo,.1,.1);
      el->SetFillColor(kBlack);
      a1->Draw();
      el->Draw();
      xaxis.Draw();
      yaxis.Draw();
      //printf("%10f %10f %10d\n %10d %10d \n",xpo,ypo,event,x,y);
      fp->c1->Update();
    }
  if (event == 11)
    {
      double xpoa1,ypoa1,xv,yv;
      fp->c1->Clear();
      xpoa=x;
      ypoa=y;
      xpoa = (xpoa-400)/400;
      ypoa = (ypoa-400)/400;
      xpoa1 = x;
      ypoa1 = y;
      axisx = new TArrow(-1,0,1,0,.05,"|>");
      axisy = new TArrow(0,-1,0,1,.05,"|>");
      axisx->Draw();
      axisy->Draw();
      xaxis.Draw();
      yaxis.Draw();
      a1= new TArrow(xpo,-ypo,xpoa,-ypoa,.07,"|>");
      a1->SetLineWidth(10);
      a1->SetFillColor(kRed);
      a1->SetLineColor(kRed);
      el = new TEllipse(xpo,-ypo,.1,.1);
      el->SetFillColor(kBlack);
      xv = (xpoa1 - xpos)/40;
      yv = -(ypoa1 - ypos)/40;
      xpos = (xpos -400)/20;
      ypos = -(ypos -400)/20;
      
      a1->Draw();
      el->Draw();
      //printf("%10f %10f %10d\n %10d %10d \n",xpo,ypo,event,x,y);
      // cout << xpos << "  ,  " << ypos << "  ,  " <<xv<<  "  ,  " <<yv  << "  ,  " <<endl;
      fp->c1->Update();
      fp->c1->Disconnect();
      fp->Setxpos(xpos);
      fp->Setypos(ypos);
      fp->Setxvel(xv);
      fp->Setyvel(yv);
      fp->c1->Close();
      fp->GetTimer()->TurnOn();
      //char xpoc[10],ypoc[10],xvec[10],yvec[10];
      //sprintf(xpoc,"xpos = %g",xpos);
      //sprintf(ypoc,"ypos = %g",ypos);
      //sprintf(xvec,"xvel = %g",xv);
      //sprintf(yvec,"yvel = %g",yv);
	
      //TText xx1(0,-0,xpoc);
      //TText xx2(0,-0.7,ypoc);
      //TText xx3(0,-0.8,xvec);
      //TText xx4(0,-0.9,yvec);
      //xx1.Draw();
      //xx2.Draw();
      //xx3.Draw();
      //xx4.Draw();
      
      
    }
           
}

 
void foucault_pendulum()
{
  fp = new testframe();
}

	   
  
  
