/*
coded by Clayton Otey
supercollider wrapped by Victor Bombi
glitched by Nathan Ho
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "piano.h"
#include <string.h>
#include "dwglib/DWG.cpp"
///////////////////////////////////////////////////////
void Piano::xormask(int xormask) {
  for (int k = 0; k < nstrings; k++) {
    string[k]->fracdelay3.xormask = xormask;
  }
}
void Piano::trigger(float v){
    this->v0 = v;
    hammer->trigger(v);
}
long Piano :: go(float *out, int samples)
{
  long n = 0;
  float Zx2  = 2*Z;
  float facZ = Zx2/(Z*nstrings+Zb);

  for(int i=0;i<samples;i++) {
   // if(sample > this->samples)
   //   break;
    n++;
    float vstring = 0.0;
    for(int k=0;k<nstrings;k++) {
      vstring += string[k]->input_velocity();
    }

    float hload = hammer->load(vstring/nstrings);
    float load = 0;
    /*
    for(int k=0;k<nstrings;k++) {
      load += (2*Z*string[k]->go_hammer(hload/(2*Z)))/(Z*nstrings+Zb);
    }
    */
    for(int k=0;k<nstrings;k++) {
      load += string[k]->go_hammer(hload/(Zx2));
    }
    load *=facZ;
    float output = 0.0;
    for(int k=0;k<nstrings;k++) {
      output += string[k]->go_soundboard(load);
    }

    //output = soundboard->reverb(output);
    /*
    output += filter(output,&shaping1);
    output = filter(output,&shaping2);
    output += filter(output,&shaping3);
*/
    //t += dt;
    //sample++;
    out[i] = output * 100;
  }
  return n;
}
float F_NOTE_31 = 6.875 * pow(2.0,(31.0 + 3)/12.0);
float F_NOTE_41 = 6.875 * pow(2.0,(41.0 + 3)/12.0);



Piano :: Piano(Unit * unit):unit(unit)
{

}

float linearmap(float s,float e,float ds,float de,float v){
    return ((de-ds)*(v-s)/(e-s)) + ds;
}

float sigmoidal(float midi,float minV,float maxV,float ampL,float ampR){
    //float ampL = amp * center;
    //float ampR = amp * (1 - center);
    //Print("sigmoidal %g %g %g %g\n",minV,maxV, ampL, ampR);
    float i=linearmap(21,108,ampL,ampR,midi);

    float minV2 = 1/(1+exp(ampR));
    float maxV2 = 1/(1+exp(ampL));

    float escale = (maxV - minV)/(maxV2 - minV2);

    float offset = maxV - maxV2 * escale;

    float val = offset + escale/(1 + exp(i));
    return val;
}
void Piano :: init(float f, float Fs, float velocity, float minr,float maxr,float amprl,float amprr, float mult_radius_core_string, float minL,float maxL,float ampLl,float ampLr, float mult_density_string, float mult_modulus_string, float mult_impedance_bridge, float mult_impedance_hammer, float mult_mass_hammer, float mult_force_hammer, float mult_hysteresis_hammer, float mult_stiffness_exponent_hammer, float position_hammer, float mult_loss_filter,float detune,int hammer_type)
{
  //this->amp = amp * 100.0;
  float f0 = 27.5;
  float midinote = 12 * log(f/f0)/log(2) + 21;
  //float L = .04 + 1.4/(1+exp(-3.4+1.4*log(f/f0)));
  float L = sigmoidal(midinote,minL,maxL,ampLl,ampLr);
  //float r = .002*pow(1+.6*log(f/f0),-1.4);
  float r = sigmoidal(midinote,minr,maxr,amprl,amprr)*0.001;
  float rho = 7850.0;
  rho *= mult_density_string;
  float logff0 = log(f/f0)/log(4192/f0);
  logff0 = sc_max(0,logff0);
  float m = .06 - .058*pow(logff0,.1);
  m *= mult_mass_hammer;
  float alpha = 0.1e-4*logff0;
  alpha *= mult_hysteresis_hammer;
  float p = 2.0+1.0*logff0;
  p *= mult_stiffness_exponent_hammer;
  float K = 40.0 / pow(.7e-3,p);
  K *= mult_force_hammer;
  Zb = 4000.0 * mult_impedance_bridge;
  Zh = 1.0 * mult_impedance_hammer;
  float rcore = (r<.0006)?r:.0006;
  rcore *= mult_radius_core_string;
  float E = 200e9;
  E *= mult_modulus_string;
  float pos = position_hammer;

  this->Fs = Fs;
  this->v0 = velocity;
  //sample = 0;
  //t = 0.0;
  //dt = 1.0/Fs;
  this->f = f;

    if(f < F_NOTE_31)
        nstrings = 1;
    else if(f < F_NOTE_41)
        nstrings = 2;
    else
        nstrings = 3;

  float c1 = 0.25;
  float c3 = 5.85;
  c3 *= mult_loss_filter;
  //float c1b = 20.0;
  //float c3b = 20.0;
  float rhoL = PI*r*r*rho;
  float T = (2*L*f)*(2*L*f)*rhoL;
  Z = sqrt(T*rhoL);
  float B = (PI*PI*PI)*E*rcore*rcore*rcore*rcore/(4.0*L*L*T);

  for(int k=0;k<nstrings;k++) {
    string[k] = new dwgs(f*(1.0+TUNE[k]*detune),Fs,pos,c1,c3,B,Z,Zb+(nstrings-1)*Z,Zh,unit);
  }


  switch (hammer_type){
      case 1:
            hammer = new StulovHammer(f,Fs,m,K,p,Z,alpha,v0);
            break;
      case 2:
            hammer = new BanksHammer(f,Fs,m,K,p,Z,alpha,v0);
            break;
      default:
        hammer = new StulovHammer(f,Fs,m,K,p,Z,alpha,v0);
  }

//Print("f = %g, r = %g mm, L = %g, T = %g, hammer = %g, Z = %g, K = %g, B = %g midi%g alpha=%g p= %g m= %g\n",f,1000*r,L,T,pos,Z,K,B,midinote,alpha,p,m);
//Print("f = %g, r = %g, rold = %g mm, L = %g, Lold = %g, T = %g, Z = %g\n",f,1000*r,1000*rold,L,Lold,T,Z);
}

Piano :: ~Piano() {
  for(int k=0;k<nstrings;k++) {
    delete string[k];
  }

  delete hammer;

}

struct GlitchPiano : public Unit
{
    Piano piano;
    GlitchPiano(Unit *unit):piano(unit){};
    int xormaskpos;
};
extern "C"
{
    void GlitchPiano_next(GlitchPiano *unit, int inNumSamples);
    void GlitchPiano_Ctor(GlitchPiano* unit);
    void GlitchPiano_Dtor(GlitchPiano* unit);
}

void GlitchPiano_Ctor(GlitchPiano* unit) {
    int inpos=0;
    float freq = ZIN0(inpos++);
    float velocity = ZIN0(inpos++);
    float gate = ZIN0(inpos++);
    float minr = ZIN0(inpos++);
    float maxr = ZIN0(inpos++);
    float amprl = ZIN0(inpos++);
    float amprr = ZIN0(inpos++);
    float rcore = ZIN0(inpos++);
    float minl = ZIN0(inpos++);
    float maxl = ZIN0(inpos++);
    float ampll = ZIN0(inpos++);
    float amplr = ZIN0(inpos++);
    float rho = ZIN0(inpos++);
    float young = ZIN0(inpos++);
    float zb = ZIN0(inpos++);
    float zh = ZIN0(inpos++);
    float mh = ZIN0(inpos++);
    float k = ZIN0(inpos++);
    float alpha = ZIN0(inpos++);
    float p = ZIN0(inpos++);
    float ph = ZIN0(inpos++);
    float loss = ZIN0(inpos++);
    float detune = ZIN0(inpos++);
    int hammer_type = ZIN0(inpos++);
    unit->xormaskpos = inpos;

    gWorld = unit->mWorld;
    //unit->piano= (Piano*) new Piano(unit);
    new(unit) GlitchPiano(unit);
    unit->piano.init(freq,SAMPLERATE,velocity*10,minr,maxr,amprl,amprr,rcore,minl,maxl,ampll,amplr,rho,young,zb,zh,mh,k,alpha,p,ph,loss,detune,hammer_type);
    SETCALC(GlitchPiano_next);
}

void GlitchPiano_Dtor(GlitchPiano* unit) {
    unit->~GlitchPiano();

}

void GlitchPiano_next(GlitchPiano *unit, int inNumSamples) {

    float * out = OUT(0);
    float velocity = ZIN0(1);
    float gate = ZIN0(2);
    int xormask = ZIN0(unit->xormaskpos);

    unit->piano.xormask(xormask);

    if (gate > 0.0)
        unit->piano.trigger(velocity*10.0);
    unit->piano.go(out,inNumSamples);

}

PluginLoad(GlitchPiano){
    ft = inTable;
    DefineDtorUnit(GlitchPiano);
}
