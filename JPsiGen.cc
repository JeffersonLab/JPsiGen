#include <TF1.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <TRandom2.h>
#include "TTCSKine.h"
#include "KinFunctions.h"
#include "RadiativeCorrections.h"
#include <TLorentzVector.h>

using namespace std;
using namespace KinFuncs;

int main() {


    // ==================================
    // ==== Reading the input config file
    // ==================================

    ifstream inpconfig("GenOptions.dat");

    map<std::string, std::string> m_Settings;
    if (inpconfig.is_open()) {
        while (!inpconfig.eof()) {
            std::string Key;
            std::string Val;
            inpconfig>>Key;
            inpconfig>>Val;
            m_Settings[Key] = Val;
            //cout<<setw(10)<<Key<<setw(20)<<m_Settings[Key]<<endl;
        }
    } else {
        cout << "Can not open the file GenOptions.dat" << endl;
        cout << "So can not initialize settings " << endl;
        cout << "Exiting" << endl;
        exit(1);
    }

    int Nsim;
    double Eb;
    double t_lim;
    double Eg_min;
    double Eg_max;
    bool isLund;
    bool Rad_corr;
    double q2_cut;
    double l_targ;
    double tSlope;
    int n_perfile;
    int seed;
    double vz_max;
    double vz_min;
    double rad_cut_off_min;
    double rad_cut_off_max;


    for (map<std::string, std::string>::iterator it = m_Settings.begin(); it != m_Settings.end(); it++) {

        std::string key = (*it).first;
        std::string val = (*it).second;

        if (key.compare("Nsim") == 0) {
            Nsim = atoi(val.c_str());
        } else if (key.compare("NPerFile") == 0) {
            n_perfile = atoi(val.c_str());
        } else if (key.compare("Eb") == 0) {
            Eb = atof(val.c_str());
        } else if (key.compare("tLim") == 0) {
            t_lim = atof(val.c_str());
        } else if (key.compare("EgMin") == 0) {
            Eg_min = atof(val.c_str());
        } else if (key.compare("EgMax") == 0) {
            Eg_max = atof(val.c_str());
        } else if (key.compare("Q2Cut") == 0) {
            q2_cut = atof(val.c_str());
        } else if (key.compare("tSlope") == 0) {
            tSlope = atof(val.c_str());
        } else if (key.compare("lTarg") == 0) {
            l_targ = atof(val.c_str());
        } else if (key.compare("LUND") == 0) {
            isLund = atof(val.c_str());
        } else if (key.compare("RAD_CORR") == 0) {
            Rad_corr = atof(val.c_str());
        } else if (key.compare("Seed") == 0) {
            seed = atoi(val.c_str());
        } else if( key.compare("vzMax") == 0 ){
            vz_max = atof(val.c_str());
        }else if( key.compare("vzMin") == 0 ){
            vz_min = atof(val.c_str());
        }else if( key.compare("rad_cut_off_min") == 0 ){
            rad_cut_off_min = atof(val.c_str());
        }else if( key.compare("rad_cut_off_max") == 0 ){
            rad_cut_off_max = atof(val.c_str());
        }



    }

    cout << "Nsim = " << Nsim << endl;
    cout << "Eb = " << Eb << endl;
    cout << "t_lim = " << t_lim << endl;
    cout << "Eg_min = " << Eg_min << endl;
    cout << "Eg_max = " << Eg_max << endl;
    cout << "q2_cut = " << q2_cut << endl;
    cout << "vz_max = " << vz_max << endl;
    cout << "vz_min = " << vz_min << endl;
    cout << "tSlope = " << tSlope << endl;
    cout << "IsLund = " << isLund << endl;
    cout << "Rad_corr = " << Rad_corr << endl;
    cout << "rad_cut_off_min = " << rad_cut_off_min << endl;
    cout << "rad_cut_off_max = " << rad_cut_off_max << endl;
    cout<<"**************************************************"<<endl;
    cout<<"*******"<<" RandomSeedActuallyUsed: "<<seed<<" *******"<<endl;
    cout<<"**************************************************"<<endl;


    const double PI = 3.14159265358979312;
    const double radian = 57.2957795130823229;
    const double Mp = 0.9383;
    const double Me = 0.00051;
    const double MJPsi = 3.097;

    Eg_min = MJPsi * (MJPsi + 2 * Mp) / (2 * Mp); //Gev
    Eg_max = Eb; //GeV
    //  const double Minv_min = sqrt(Mp*Mp + 2*Mp*Eg_min ) - Mp;
    const double Q2min = 2 * Mp * Eg_min + t_lim - (Eg_min / Mp)*(2 * Mp * Mp - t_lim - sqrt(t_lim * t_lim - 4 * Mp * Mp * t_lim));
    const double Minv_min = sqrt(Q2min);
    const double SLAC_Fit_scale = 7.79117e-23;

    bool write_root = !isLund;

    TRandom2 rand;
    rand.SetSeed(seed);

    TF1 *f_JPsi_dsigm_dt = new TF1("f_JPsi_dsigm_dt", JPsi_dsigm_dt, 8.25, 25., 3);
    f_JPsi_dsigm_dt->SetParameter(1, SLAC_Fit_scale);
    f_JPsi_dsigm_dt->SetParameter(2, tSlope);

    TTCSKine tcs_kin1(Mp, Eb);

    TLorentzVector target(0., 0., 0., Mp);
    TLorentzVector Lcm;

    TFile *file_out;
    ofstream Lund_out;
    int file_number = 0;
    if (!isLund) {
        file_out = new TFile("JPsi_gen.root", "Recreate");
    } else {
        //Lund_out.open(Form("JPsi_gen_%d.txt", file_number), ofstream::out);
        Lund_out.open("JPsiGen.dat", ofstream::out);
    }

    TH2D *h_ph_h_ph_cm1 = new TH2D("h_ph_h_ph_cm1", "", 200, 0., 360., 200, 0., 360.);
    TH2D *h_th_g_th_cm1 = new TH2D("h_th_g_th_cm1", "", 200, 0., 180., 200, 0., 180.);

    //================= Definition of Tree Variables =================
    double Eg, Minv, t, Q2;
    double psf, crs_BH, crs_INT, crs_int, crs_JPsi;
    double psf_flux, flux_factor;
    TLorentzVector L_em, L_ep, L_prot, L_rad;
    TLorentzVector L_gprime;

    double px_prot, py_prot, pz_prot, E_prot;
    double px_ep, py_ep, pz_ep, E_ep;
    double px_em, py_em, pz_em, E_em;
    double px_rad, py_rad, pz_rad, E_rad;
    double E_rad_cm;
    double Inv_Mass;


    TTree *tr1 = new TTree("tr1", "TCS MC events");
    tr1->Branch("L_em", "TLorentzVector", &L_em, 3200, 99);
    tr1->Branch("L_ep", "TLorentzVector", &L_ep, 3200, 99);
    tr1->Branch("L_prot", "TLorentzVector", &L_prot, 3200, 99);
    tr1->Branch("Eg", &Eg, "Eg/D");
    tr1->Branch("Q2", &Q2, "Q2/D");
    tr1->Branch("t", &t, "t/D");
    tr1->Branch("psf", &psf, "psf/D");
    tr1->Branch("flux_factor", &flux_factor, "flux_factor/D");
    tr1->Branch("crs_JPsi", &crs_JPsi, "crs_JPsi/D");
    tr1->Branch("px_prot", &px_prot, "px_prot/D");
    tr1->Branch("py_prot", &py_prot, "py_prot/D");
    tr1->Branch("pz_prot", &pz_prot, "pz_prot/D");
    tr1->Branch("E_prot", &E_prot, "E_prot/D");
    tr1->Branch("px_ep", &px_ep, "px_ep/D");
    tr1->Branch("py_ep", &py_ep, "py_ep/D");
    tr1->Branch("pz_ep", &pz_ep, "pz_ep/D");
    tr1->Branch("E_ep", &E_ep, "E_ep/D");
    tr1->Branch("px_em", &px_em, "px_em/D");
    tr1->Branch("py_em", &py_em, "py_em/D");
    tr1->Branch("pz_em", &pz_em, "pz_em/D");
    tr1->Branch("E_em", &E_em, "E_em/D");
    tr1->Branch("px_rad", &px_rad, "px_rad/D");
    tr1->Branch("py_rad", &py_rad, "py_rad/D");
    tr1->Branch("pz_rad", &pz_rad, "pz_rad/D");
    tr1->Branch("E_rad", &E_rad, "E_rad/D");
    tr1->Branch("E_rad_cm", &E_rad_cm, "E_rad_cm/D");
    tr1->Branch("Inv_Mass", &Inv_Mass, "Inv_Mass/D");




    /////////////Set Rad Corr Parameters//////////////
    RadiativeCorrections Rad_corr_1(rad_cut_off_min, rad_cut_off_max);


    /////////////////////////////////////////////////


    for (int i = 0; i < Nsim; i++) {

        if (i % 50000 == 0) {
            cout.flush() << "Processed " << i << " events, approximetely " << double(100. * i / double(Nsim)) << "%\r";
        }

        //        if ((i + 1) % n_perfile == 0) {
        //            if (isLund) {
        //                Lund_out.close();
        //                file_number++;
        //                //Lund_out.open(Form("JPsi_gen_%d.txt", file_number), ofstream::out);
        //                Lund_out.open(Form("JPsi_gen_%d.txt", file_number), ofstream::out);
        //            }
        //        }

        Q2 = MJPsi*MJPsi;

        double psf_Eg = Eg_max - Eg_min;
        Eg = rand.Uniform(Eg_min, Eg_min + psf_Eg);
        flux_factor = N_EPA(Eb, Eg, q2_cut) + N_Brem(Eg, Eb);
        
        double s = Mp * Mp + 2 * Mp*Eg;
        double t_min = T_min(0., Mp*Mp, Q2, Mp*Mp, s);
        double t_max = T_max(0., Mp*Mp, Q2, Mp*Mp, s);
        double psf_t = t_min - TMath::Max(t_max, t_lim);

        if (t_min > t_lim) {
            t = rand.Uniform(t_min - psf_t, t_min);

            float fl_s = (float) s;
            float fl_t = (float) t;

            //cout<<"==== JPSi cross section   "<<jpsi_dsdt_(&fl_s, &fl_t)<<endl;

            //crs_JPsi = JPsi_dsdt(s, t);
            f_JPsi_dsigm_dt->SetParameter(0, Eg);
            crs_JPsi = f_JPsi_dsigm_dt->Eval(t);
            //cout<<"Q2 = "<<Q2<<"     JPsi dSdt = "<<crs_JPsi<<endl;


            double u = 2 * Mp * Mp + Q2 - s - t;
            double th_qprime = acos((s * (t - u) - Mp * Mp * (Q2 - Mp * Mp)) / sqrt(Lambda(s, 0, Mp * Mp) * Lambda(s, Q2, Mp * Mp))); //Byukling Kayanti (4.9)
            double th_pprime = PI + th_qprime;

            double Pprime = 0.5 * sqrt(Lambda(s, Q2, Mp * Mp) / s); // Momentum in c.m. it is the same for q_pr and p_pr

            Lcm.SetPxPyPzE(0., 0., Eg, Mp + Eg);
            L_prot.SetPxPyPzE(Pprime * sin(th_pprime), 0., Pprime * cos(th_pprime), sqrt(Pprime * Pprime + Mp * Mp));
            L_gprime.SetPxPyPzE(Pprime * sin(th_qprime), 0., Pprime * cos(th_qprime), sqrt(Pprime * Pprime + Q2));

            double psf_cos_th = 2.; // cos(th):(-1 : 1)
            double psf_phi_cm = 2 * PI;

            double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
            double sin_th = sqrt(1 - cos_th * cos_th);
            double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);

            double El = sqrt(Q2) / 2.; // Energy of lepton in the rest frame of qprime
            double Pl = sqrt(El * El - Me * Me);

            
            L_em.SetPxPyPzE(Pl * sin_th * cos(phi_cm), Pl * sin_th * sin(phi_cm), Pl*cos_th, El);
            L_ep.SetPxPyPzE(-Pl * sin_th * cos(phi_cm), -Pl * sin_th * sin(phi_cm), -Pl*cos_th, El);

            /*cout<<" "<<endl;
            cout<<"New Event "<<endl;
            cout<<"Before radiation"<<endl;
            cout<<"P electron"<<L_em.P()<<" E electron "<<L_em.E()<<endl;
            cout<<"P positron"<<L_ep.P()<<" E positron "<<L_ep.E()<<endl;
*/

            bool in_rad_tail = (rand.Uniform(0,1) > Rad_corr_1.Compute_cs_correction_factor(sqrt(Q2))); // randomly choose if the photon is above cut_off_min
            //if(Rad_corr && in_rad_tail)

            if(Rad_corr && in_rad_tail)
                Rad_corr_1.Soft_Photon_Emission(L_em, L_ep, L_rad);
            else 
                L_rad.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

            E_rad_cm = L_rad.E();

            //cout<<Rad_corr_1.Compute_cs_correction_factor(sqrt(Q2))<<endl;

            /*cout<<"After radiation"<<endl;
            cout<<"P electron"<<L_em.P()<<" E electron "<<L_em.E()<<endl;
            cout<<"P positron"<<L_ep.P()<<" E positron "<<L_ep.E()<<endl;
            cout<<"E photon"<<L_rad.P()<<endl;*/

            L_em.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_ep.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_rad.RotateY(th_qprime);

            L_em.Boost(L_gprime.BoostVector()); // Move to the CM Frame
            L_ep.Boost(L_gprime.BoostVector()); // Move to the CM Frame
            L_rad.Boost(L_gprime.BoostVector()); // Move to the CM Frame

            L_em.Boost(Lcm.BoostVector()); // Move to the Lab Frame
            L_ep.Boost(Lcm.BoostVector()); // Move to the Lab Frame
            L_rad.Boost(Lcm.BoostVector()); // Move to the Lab Frame



            L_gprime.Boost(Lcm.BoostVector());
            L_prot.Boost(Lcm.BoostVector());


            //cout<<"test = "<<tcs_kin1.GetMM2()<<endl;
            //cout<<"phi_cm = "<<phi_cm*TMath::RadToDeg()<<"  tcs_kine.GetPhi_cm()"<<tcs_kin1.GetPhi_cm()<<endl;
            //cout<<"phi_cm = "<<phi_cm*TMath::RadToDeg()<<"  tcs_kine.GetPhi_cm()"<<tcs_kin1.GetPhi_cm()<<endl;

            double psf_phi_lab = 2 * PI;
            double phi_rot = rand.Uniform(0., psf_phi_lab);

            L_prot.RotateZ(phi_rot);
            L_gprime.RotateZ(phi_rot);
            L_em.RotateZ(phi_rot);
            L_ep.RotateZ(phi_rot);
            L_rad.RotateZ(phi_rot);

            tcs_kin1.SetLemLepLp(L_em, L_ep, L_prot);

            h_ph_h_ph_cm1->Fill(phi_cm * TMath::RadToDeg(), tcs_kin1.GetPhi_cm());
            h_th_g_th_cm1->Fill(acos(cos_th) * TMath::RadToDeg(), tcs_kin1.GetTheta_cm());

            psf = psf_t * psf_Eg;// * psf_phi_lab * psf_cos_th * psf_phi_cm ;

            double eta = Q2 / (2 * (s - Mp * Mp) - Q2);

            double vz = rand.Uniform( vz_min, vz_max);

            px_prot = L_prot.Px();
            py_prot = L_prot.Py();
            pz_prot = L_prot.Pz();
            E_prot = L_prot.E();
            px_ep = L_ep.Px();
            py_ep = L_ep.Py();
            pz_ep = L_ep.Pz();
            E_ep = L_ep.E();
            px_em = L_em.Px();
            py_em = L_em.Py();
            pz_em = L_em.Pz();
            E_em = L_em.E();
            px_rad = L_rad.Px();
            py_rad = L_rad.Py();
            pz_rad = L_rad.Pz();
            E_rad = L_rad.E();

            Inv_Mass=(L_em+L_ep).M();

            double tot_weight = crs_JPsi * psf*flux_factor;

            if (write_root) {
                tr1->Fill();
            } else {
                // Writing Header
                Lund_out << 4 << setw(5) << 1 << setw(5) << 1 << setw(15) << psf << setw(15) << crs_JPsi << setw(15) << 0 << setw(15) << flux_factor << setw(15) << crs_JPsi << setw(15) << psf << setw(15) << tot_weight << endl;
                // Writing Proton
                Lund_out << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_prot << setw(15) << py_prot << setw(15) << pz_prot;
                Lund_out << setw(15) << L_prot.E() << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
                // Writing Electron
                Lund_out << 2 << setw(5) << -1 << setw(5) << 1 << setw(7) << 11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_ep << setw(15) << py_ep << setw(15) << pz_ep;
                Lund_out << setw(15) << L_em.E() << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
                // Writing Positron
                Lund_out << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << -11 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_em << setw(15) << py_em << setw(15) << pz_em;
                Lund_out << setw(15) << L_ep.E() << setw(15) << Me << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
                Lund_out << 4 << setw(5) << 0 << setw(5) << 1 << setw(7) << 22 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_rad << setw(15) << py_rad << setw(15) << pz_rad;
                Lund_out << setw(15) << L_rad.E() << setw(15) << 0.0 << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
            }
        } else {
            cout << " |t_min| > |t_lim|" << endl;
            cout << " t_min =  " << t_min << "   t_lim = " << t_lim << "  Eg = " << Eg << endl;
        }

    }



    if (write_root) {
        tr1->Write();
        h_ph_h_ph_cm1->Write();
        h_th_g_th_cm1->Write();
        file_out->Close();
    }
}
