#include <cmath>
#include <TMath.h>
#include <RadiativeCorrections.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>

#include <iostream>

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

RadiativeCorrections::RadiativeCorrections(double cut_off_min_in, double cut_off_max_in)
{
  cut_off_min = cut_off_min_in;
  cut_off_max = cut_off_max_in;
  cs_function = new TF1("cs_emmision_photon_func", cs_emmision_photon, cut_off_min, cut_off_max, 2);

  cs_correction_factor_func = new TF1("cs_correction_factor_func", radiative_correction_factor, 1.0, 4.0, 2);
  cs_correction_factor_func->SetParameter(1, cut_off_min);
}

double  RadiativeCorrections::Compute_cs_correction_factor(double Inv_Mass)
{
  cs_correction_factor = cs_correction_factor_func->Eval(Inv_Mass);
  return cs_correction_factor;
};

//Old version with emmision along one lepton
/*void RadiativeCorrections::Soft_Photon_Emission(TLorentzVector &Electron, TLorentzVector &Positron, TLorentzVector &Rad_Photon)
{

  /// use pointer instead to modify electron or positron
  double me = 0.00051;
  TRandom3 rand;
  rand.SetSeed(0); // to be checked

  double M_lepton_pair = (Electron + Positron).M();

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
  Lepton_Radiated = Lepton_Radiated - Rad_Photon;

  if (Choice_of_lepton == 0)
    Electron = Lepton_Radiated;
  else
    Positron = Lepton_Radiated;

  return;
}*/


//////// V2 //////////////////
/*void RadiativeCorrections::Soft_Photon_Emission(TLorentzVector &Electron, TLorentzVector &Positron, TLorentzVector &Rad_Photon_1, TLorentzVector &Rad_Photon_2)
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
  double theta_radiated = rand.Uniform(0,PI);
  double phi_radiated = rand.Uniform(0,2.*PI);

  //Rad_Photon.SetXYZM(E_s * cos(phi_radiated) * sin(theta_radiated), E_s * sin(phi_radiated) * sin(theta_radiated), E_s * cos(theta_radiated), 0.0);

  double frac_energy = rand.Uniform(0.,1.); //Uniformly draw a fraction of energy radited by the electron

  Rad_Photon_1.SetXYZM(E_s * frac_energy * cos(phi_radiated) * sin(theta_radiated), E_s * frac_energy * sin(phi_radiated) * sin(theta_radiated), E_s * frac_energy * cos(theta_radiated), 0.0);
  Rad_Photon_2.SetXYZM(E_s * (1.-frac_energy) * cos(phi_radiated) * sin(theta_radiated), E_s * (1.-frac_energy) * sin(phi_radiated) * sin(theta_radiated), E_s * (1.-frac_energy) * cos(theta_radiated), 0.0);

  Electron = Electron - Rad_Photon_1;
  Positron = Positron - Rad_Photon_2;

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
  double theta_radiated_1 = rand.Uniform(0,PI);
  double phi_radiated_1 = rand.Uniform(0,2.*PI);

  double theta_radiated_2 = rand.Uniform(0,PI);
  double phi_radiated_2 = rand.Uniform(0,2.*PI);

  //Rad_Photon.SetXYZM(E_s * cos(phi_radiated) * sin(theta_radiated), E_s * sin(phi_radiated) * sin(theta_radiated), E_s * cos(theta_radiated), 0.0);

  double frac_energy = rand.Uniform(0.,1.); //Uniformly draw a fraction of energy radited by the electron

  Rad_Photon_1.SetXYZM(E_s * frac_energy * cos(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * sin(phi_radiated_1) * sin(theta_radiated_1), E_s * frac_energy * cos(theta_radiated_1), 0.0);
  Rad_Photon_2.SetXYZM(E_s * (1.-frac_energy) * cos(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * sin(phi_radiated_2) * sin(theta_radiated_2), E_s * (1.-frac_energy) * cos(theta_radiated_2), 0.0);

  Electron = Electron - Rad_Photon_1;
  Positron = Positron - Rad_Photon_2;

  return;
}


RadiativeCorrections::~RadiativeCorrections()
{
}
