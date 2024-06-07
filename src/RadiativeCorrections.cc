#include <cmath>
#include <TMath.h>
#include <RadiativeCorrections.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

#include <iostream>

/////////////////////////////////////////////////////////////
///////// Define functions for radiative corrections ////////
/////////////////////////////////////////////////////////////
double RadiativeCorrections::R(double rm)
{
    double R = (rm - (1. / 4.) * rm * rm - (3. / 4.)) * (std::log((Mv / me) * std::sqrt(1. - rm)) - (1. / 2.));
    return R; // term1 * term2;
}

double RadiativeCorrections::E_to_rm(double E)
{
    return E / (Mv / 2.);
}

// Define the main function
double RadiativeCorrections::VM_delta(double rm)
{
    double delta = -(4. / 137.) / PI * ((1. / 2. - log(Mv / me)) * (log(rm) + 1. / 3.) + (1. / 2.) * (TMath::DiLog(rm) - (pow(PI, 2)) / 6.) + 17. / 72. + R(rm));
    return delta;
}

double RadiativeCorrections::VM_delta_fun(double *x, double *par)
{
    return VM_delta(x[0]);
}

double RadiativeCorrections::constant_VM_delta(double rm)
{
    return VM_delta(0.999);
}

double RadiativeCorrections::d_VM_delta_fun(double *x, double *par)
{
    double h = 1e-5;
    return (VM_delta(x[0] + h) - VM_delta(x[0] - h)) / (2.0 * h);
}

double RadiativeCorrections::varying_VM_delta(double rm)
{
    return VM_delta(rm) - VM_delta(0.999);
}

double RadiativeCorrections::varying_VM_delta_fun(double *x, double *par)
{
    return varying_VM_delta(x[0]);
}
/////////////////////////////////////////////////////////////




double cs_emmision_photon(double *xx, double *par)
{

  double s_ll = par[1] * par[1];
  double Es = xx[0];

  double me = 0.00051;
  double alpha = 1. / 137.;
  double PI = 3.14159265358979312;
  double beta = sqrt(1 - (4. * me * me / s_ll));
  double F = 1. - (alpha * alpha / 3.) * pow((1 + ((1. + beta * beta) / (2. * beta)) * log((1. - beta) / (1. + beta))), 2); /////equal to 1

  double delta = (-alpha / PI) * (log((4. * Es * Es) / (s_ll)) * (1. + log((me * me) / (s_ll))) - (PI * PI / 3.));

  double delta_prime = (-alpha / PI) * (1. + log((me * me) / (s_ll))) * (2. / Es);

  double cs_photon = delta_prime * F * exp(delta);

  return cs_photon;
}

double radiative_correction_factor(double *xx, double *par)
{

  double s_ll = xx[0] * xx[0];
  double Delta_E = par[1];

  // Forumla (66) in https://journals.aps.org/prd/pdf/10.1103/PhysRevD.97.076012
  double me = 0.00051;
  double alpha = 1. / 137.;
  double PI = 3.14159265358979312;
  double beta = sqrt(1 - (4. * me * me / s_ll));
  double F = 1. - (alpha * alpha / 3.) * pow((1 + ((1. + beta * beta) / (2. * beta)) * log((1. - beta) / (1. + beta))), 2); /////equal to 1

  double argument_exponential = (-alpha / PI) * (log((4. * Delta_E * Delta_E) / (me * me)) + log((1. - beta) / (1. + beta))) * (1. + ((1. + beta * beta) / (2. * beta)) * log((1. - beta) / (1. + beta))); //-(alpha/PI)*( ( (1-beta)/beta )*log( (1.-beta)/(1.+beta) ) + ( (1+beta*beta)/(2.*beta) )*(4.*TMath::DiLog( (2.*beta)/(1+beta) ) - PI*PI) )

  double third_factor = 1. - (alpha / PI) * (((1 - beta) / beta) * log((1. - beta) / (1. + beta)) + ((1 + beta * beta) / (2. * beta)) * (4. * TMath::DiLog((2. * beta) / (1 + beta)) - PI * PI));

  // double rad_corr_factor = F*exp(argument_exponential)*third_factor;
  double delta = (-alpha / PI) * (log((4. * Delta_E * Delta_E) / (s_ll)) * (1. + log((me * me) / (s_ll))) - (PI * PI / 3.));

  double rad_corr_factor = F * exp(delta);

  return rad_corr_factor;
}

RadiativeCorrections::RadiativeCorrections(double cut_off_min_in)
{
  cut_off_min = cut_off_min_in;
  
  cs_function = new TF1("cs_emmision_photon_func", cs_emmision_photon, cut_off_min, cut_off_max, 2);

  cs_correction_factor_func = new TF1("cs_correction_factor_func", radiative_correction_factor, 1.0, 4.0, 2);
  cs_correction_factor_func->SetParameter(1, cut_off_min);

  fun_varying_VM_delta = new TF1("varying_VM_delta_fun", varying_VM_delta_fun, 0.0, 0.999, 0);
}

void RadiativeCorrections::Set_Inv_Mass(double Inv_Mass)
{
cut_off_max = Inv_Mass/2.0;
}

double  RadiativeCorrections::Compute_cs_correction_factor(double Inv_Mass)
{
  cs_correction_factor = cs_correction_factor_func->Eval(Inv_Mass);

  
//  TF1* cs_correction_factor_func1 = new TF1("cs_correction_factor_func", radiative_correction_factor, 1.0, 4.0, 2);
//  cs_correction_factor_func1->SetParameter(1, cut_off_max);

/*std::cout<<std::endl;
  std::cout<<cs_correction_factor<<std::endl;
  std::cout<<cs_correction_factor_func1->Eval(Inv_Mass)<<std::endl;*/
  //cut_off_max = min(Inv_Mass/2.0,cut_off_max); //Make sure the maximum energy emmited is equal to the momentum of on of the lepton in the CoM

  return cs_correction_factor;
};

double  RadiativeCorrections::Compute_cs_correction_factor_VM()
{
  
  double r_m_min = E_to_rm(cut_off_min);
  double int_cs_from_zero_to_E_s_min = fun_varying_VM_delta->Eval(r_m_min);

  return int_cs_from_zero_to_E_s_min;
};

//Old version with emmision along one lepton
/*void RadiativeCorrections::Soft_Photon_Emission(TLorentzVector &Electron, TLorentzVector &Positron, TLorentzVector &Rad_Photon, TLorentzVector &Rad_Photon1)
{

  /// use pointer instead to modify electron or positron
  double me = 0.00051;
  TRandom3 rand;
  rand.SetSeed(0); // to be checked

  double M_lepton_pair = (Electron + Positron).M();
  cut_off_max = M_lepton_pair/2.0; //Make sure the maximum energy emmited is equal to the momentum of on of the lepton in the CoM


  TF1 *cs_emmision_photon_func = new TF1("cs_emmision_photon_func", cs_emmision_photon, cut_off_min, cut_off_max, 2);
  cs_emmision_photon_func->SetParameter(1, M_lepton_pair);
  double E_s = cs_emmision_photon_func->GetRandom();

  ////Randomly choose between electron or positron to emit a photon////
  int Choice_of_lepton = rand.Integer(2);
  // std::cout<<Choice_of_lepton<<std::endl;
  TLorentzVector Lepton_Radiated = (Choice_of_lepton == 0) ? Electron : Positron;
  ////// get theta and phi of the lepton that emitted the photon
  double theta_radiated = Lepton_Radiated.Theta();
  double phi_radiated = Lepton_Radiated.Phi();

  Rad_Photon.SetXYZM(E_s * cos(phi_radiated) * sin(theta_radiated), E_s * sin(phi_radiated) * sin(theta_radiated), E_s * cos(theta_radiated), 0.0);

  double momentum_lepton_after_rad =  sqrt(pow(Lepton_Radiated.E()-E_s,2)-me*me);
  Lepton_Radiated.SetXYZM(momentum_lepton_after_rad * cos(phi_radiated) * sin(theta_radiated), momentum_lepton_after_rad * sin(phi_radiated) * sin(theta_radiated), momentum_lepton_after_rad * cos(theta_radiated), me);

  //Lepton_Radiated = Lepton_Radiated - Rad_Photon;

  if (Choice_of_lepton == 0)
    Electron = Lepton_Radiated;
  else
    Positron = Lepton_Radiated;

  return;
}*/



void RadiativeCorrections::Soft_Photon_Emission(TLorentzVector &Electron, TLorentzVector &Positron, TLorentzVector &Rad_Photon_1, TLorentzVector &Rad_Photon_2)
{

  /// use pointer instead to modify electron or positron
  double me = 0.00051;
  TRandom3 rand;
  rand.SetSeed(0); // to be checked

  double M_lepton_pair = (Electron + Positron).M();

  TF1 *cs_emmision_photon_func = new TF1("cs_emmision_photon_func", cs_emmision_photon, cut_off_min, cut_off_max, 2);
  cs_emmision_photon_func->SetParameter(1, M_lepton_pair);
  double E_s = cs_emmision_photon_func->GetRandom();

  ////// Randomly choose theta and phi of the radiated photon
  double theta_radiated_1 = Electron.Theta();
  double phi_radiated_1 = Electron.Phi();

  double theta_radiated_2 = Positron.Theta();
  double phi_radiated_2 = Positron.Phi();

  //Rad_Photon.SetXYZM(E_s * cos(phi_radiated) * sin(theta_radiated), E_s * sin(phi_radiated) * sin(theta_radiated), E_s * cos(theta_radiated), 0.0);

  double frac_energy = rand.Uniform(0.,1.); //Uniformly draw a fraction of energy radited by the electron

  Rad_Photon_1.SetXYZM(E_s * frac_energy * cos(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * sin(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * cos(theta_radiated_1), 0.0);
  Rad_Photon_2.SetXYZM(E_s * (1.-frac_energy) * cos(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * sin(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * cos(theta_radiated_2), 0.0);


  double momentum_electron_after_rad =  sqrt(pow(Electron.E()-(frac_energy)*E_s,2)-me*me);
  double momentum_positron_after_rad =  sqrt(pow(Positron.E()-(1.-frac_energy) *E_s,2)-me*me);

  Electron.SetXYZM(momentum_electron_after_rad * cos(phi_radiated_1) * sin(theta_radiated_1), momentum_electron_after_rad * sin(phi_radiated_1) * sin(theta_radiated_1), momentum_electron_after_rad * cos(theta_radiated_1), me);
  Positron.SetXYZM(momentum_positron_after_rad * cos(phi_radiated_2) * sin(theta_radiated_2), momentum_positron_after_rad * sin(phi_radiated_2) * sin(theta_radiated_2), momentum_positron_after_rad * cos(theta_radiated_2), me);

  //Electron = Electron - Rad_Photon_1;
  //Positron = Positron - Rad_Photon_2;

  //return the additional weight per event, in this case 0.0
  //return 0.0;
}


void RadiativeCorrections::Soft_Photon_Emission_VM(TLorentzVector &Electron, TLorentzVector &Positron, TLorentzVector &Rad_Photon_1, TLorentzVector &Rad_Photon_2)
{

  double me = 0.00051;
  TRandom3 rand;
  rand.SetSeed(0); // to be checked

  double weight = constant_VM_delta(0.0);
  TF1 *d_cs_emmision_photon_func_0 = new TF1("d_cs_emmision_photon_func_0", d_VM_delta_fun, E_to_rm(cut_off_min), 0.999, 0);
  double E_s = (Mv/2.)*d_cs_emmision_photon_func_0->GetRandom();

  ////// Randomly choose theta and phi of the radiated photon
  double theta_radiated_1 = Electron.Theta();
  double phi_radiated_1 = Electron.Phi();

  double theta_radiated_2 = Positron.Theta();
  double phi_radiated_2 = Positron.Phi();

  //Rad_Photon.SetXYZM(E_s * cos(phi_radiated) * sin(theta_radiated), E_s * sin(phi_radiated) * sin(theta_radiated), E_s * cos(theta_radiated), 0.0);

  double frac_energy = rand.Uniform(0.,1.); //Uniformly draw a fraction of energy radited by the electron

  Rad_Photon_1.SetXYZM(E_s * frac_energy * cos(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * sin(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * cos(theta_radiated_1), 0.0);
  Rad_Photon_2.SetXYZM(E_s * (1.-frac_energy) * cos(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * sin(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * cos(theta_radiated_2), 0.0);


  double momentum_electron_after_rad =  sqrt(pow(Electron.E()-(frac_energy)*E_s,2)-me*me);
  double momentum_positron_after_rad =  sqrt(pow(Positron.E()-(1.-frac_energy) *E_s,2)-me*me);

  Electron.SetXYZM(momentum_electron_after_rad * cos(phi_radiated_1) * sin(theta_radiated_1), momentum_electron_after_rad * sin(phi_radiated_1) * sin(theta_radiated_1), momentum_electron_after_rad * cos(theta_radiated_1), me);
  Positron.SetXYZM(momentum_positron_after_rad * cos(phi_radiated_2) * sin(theta_radiated_2), momentum_positron_after_rad * sin(phi_radiated_2) * sin(theta_radiated_2), momentum_positron_after_rad * cos(theta_radiated_2), me);

  //return weight;
}


RadiativeCorrections::~RadiativeCorrections()
{
}
