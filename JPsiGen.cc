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
    double Eg_minUser; // Eg min defined by the user, however, 
                       // if it is below the J/Psi production threshold, then the threshold value will be used for Egmin
    double Eg_max;
    bool isLund;
    double q2_cut;
    double tSlope;
    int n_perfile;
    int seed;
    double vz_max;
    double vz_min;
    bool isFermi;

    int lepType; // Lepton Type: 0=e-/e+, 1=mu-/mu+ 

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
            Eg_minUser = atof(val.c_str());
        } else if (key.compare("EgMax") == 0) {
            Eg_max = atof(val.c_str());
        } else if (key.compare("Q2Cut") == 0) {
            q2_cut = atof(val.c_str());
        } else if (key.compare("tSlope") == 0) {
            tSlope = atof(val.c_str());
        } else if (key.compare("LUND") == 0) {
            isLund = atof(val.c_str());
        } else if (key.compare("Seed") == 0) {
            seed = atoi(val.c_str());
        } else if (key.compare("vzMax") == 0) {
            vz_max = atof(val.c_str());
        } else if (key.compare("vzMin") == 0) {
            vz_min = atof(val.c_str());
        } else if (key.compare("Fermi") == 0) {
            isFermi = atof(val.c_str());
        }else if ( key.compare("LepType") == 0 ){
            lepType = atoi(val.c_str());
        }

    }



    const double PI = 3.14159265358979312;
    const double radian = 57.2957795130823229;
    const double Mp = 0.9383;
    const double Me = 0.00051;
    const double Mmu = 0.105658;
    const double MJPsi = 3.097;
    std::map<int, double> m_leptonMass;
    m_leptonMass[0] = Me;
    m_leptonMass[1] = Mmu;
    
    std::map<int, int> m_posLeptonPID;
    std::map<int, int> m_negLeptonPID;
    m_posLeptonPID[0] = -11;
    m_posLeptonPID[1] = -13;
    m_negLeptonPID[0] = 11;
    m_negLeptonPID[1] = 13;
    
    std::map<int, std::string> m_leptonType;
    m_leptonType[0] = "e-/e+";
    m_leptonType[1] = "mu-/mu+";
    
    if (Eg_max > Eb) {
        cout<<" ******* W A R N I N G ******** "<<endl;
        cout<<" ** Eg_max is "<<Eg_max<<" GeV, while the beam energy is only "<<Eb<<" GeV. "<<endl;
        cout<<" ** As a Eg_max the beam energy "<<Eb<<" GeV will be used"<<endl;
        Eg_max = Eb; //GeV
    }

    if( !(lepType == 0 || lepType == 1) ){
        cout<<"Wrong Lepton Type is provided. The value is "<<lepType<<endl;
        cout<<"It should be 0 or 1"<<endl;
        cout<<"Exiting."<<endl;
        exit(1);
    }
        
    const double m_l = m_leptonMass[lepType];
    const int pid_pos_lep = m_posLeptonPID[lepType];
    const int pid_neg_lep = m_negLeptonPID[lepType];
    const std::string lep_name = m_leptonType[lepType];
    
    cout << "Nsim           = " << Nsim << endl;
    cout << "Eb             = " << Eb << endl;
    cout << "t_lim          = " << t_lim << endl;
    cout << "Eg_min         = " << Eg_minUser << endl;
    cout << "Eg_max         = " << Eg_max << endl;
    cout << "q2_cut         = " << q2_cut << endl;
    cout << "vz_max         = " << vz_max << endl;
    cout << "vz_min         = " << vz_min << endl;
    cout << "tSlope         = " << tSlope << endl;
    cout << "IsLund         = " << isLund << endl;
    cout << "LepType        = " << lep_name << endl;
    cout << "**************************************************" << endl;
    cout << "*******" << " RandomSeedActuallyUsed: " << seed << " *******" << endl;
    cout << "**************************************************" << endl;


    const double SLAC_Fit_scale = 7.79117e-23;

    bool write_root = !isLund;

    TRandom2 rand;
    rand.SetSeed(seed);

    TF1 *f_JPsi_dsigm_dt = new TF1("f_JPsi_dsigm_dt", JPsi_dsigm_dt, 8.25, 25., 3);
    f_JPsi_dsigm_dt->SetParameter(1, SLAC_Fit_scale);
    f_JPsi_dsigm_dt->SetParameter(2, tSlope);

    TF1 *f_FermiDistr = new TF1("f_FermiDistr", Fermi_Distribution, 0., 1, 0);

    TTCSKine tcs_kin1(Mp, Eb);

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

    TH1D *h_P_Fermi1 = new TH1D("h_P_Fermi1", "", 200, 0., 1.05);

    //================= Definition of Tree Variables =================
    double Eg, Minv, t, Q2;
    double psf, crs_BH, crs_INT, crs_int, crs_JPsi;
    double psf_flux, flux_factor;
    TLorentzVector L_lm, L_lp, L_prot;
    TLorentzVector L_ProtFermi, L_gamma;
    TLorentzVector L_gprime;

    double px_prot, py_prot, pz_prot, E_prot;
    double px_lp, py_lp, pz_lp, E_lp;
    double px_lm, py_lm, pz_lm, E_lm;


    TTree *tr1 = new TTree("tr1", "TCS MC events");
    tr1->Branch("L_lm", "TLorentzVector", &L_lm, 3200, 99);
    tr1->Branch("L_lp", "TLorentzVector", &L_lp, 3200, 99);
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
    tr1->Branch("px_lp", &px_lp, "px_lp/D");
    tr1->Branch("py_lp", &py_lp, "py_lp/D");
    tr1->Branch("pz_lp", &pz_lp, "pz_lp/D");
    tr1->Branch("px_lm", &px_lm, "px_lm/D");
    tr1->Branch("py_lm", &py_lm, "py_lm/D");
    tr1->Branch("pz_lm", &pz_lm, "pz_lm/D");

    for (int i = 0; i < Nsim; i++) {

        if (i % 50000 == 0) {
            cout.flush() << "Processed " << i << " events, approximetely " << double(100. * i / double(Nsim)) << "%\r";
        }

        Q2 = MJPsi*MJPsi;
        
        // Check if Fermi option is active, if so generrate Fermi momentum for proton,
        // Otherwise the proton is at rest
         // Let's take it in the range of 0 to 1 GeV
        double p_prot_Fermi = isFermi ? f_FermiDistr->GetRandom(0., 1.) : 0;

        h_P_Fermi1->Fill(p_prot_Fermi);

        double cosThFermi = rand.Uniform(-1., 1);
        double sinThFermi = sqrt(1. - cosThFermi * cosThFermi);
        double phiFermi = rand.Uniform(0, 2 * PI);

        double pxFermi = p_prot_Fermi * sinThFermi * cos(phiFermi);
        double pyFermi = p_prot_Fermi * sinThFermi * sin(phiFermi);
        double pzFermi = p_prot_Fermi*cosThFermi;
        double EFermi = sqrt(p_prot_Fermi * p_prot_Fermi + Mp * Mp);

        //***** The Eg Threshld in the Lab frae is calculated as  *****
        double EgMin_CM = (MJPsi * MJPsi + 2 * Mp * MJPsi) / (2 * (Mp + MJPsi));
        double Eg_thr = (MJPsi * MJPsi + 2 * Mp * MJPsi) / (2 * (EFermi - p_prot_Fermi * cosThFermi));
                
        double Eg_min = TMath::Max(Eg_minUser, Eg_thr);
        
        // When it happens that because of the Fermi momentum the Eg threshold becomes higher
        // than the Eg_max, then this is "not possible (or undesired) kinematics", so le't skip this event
        if( Eg_thr > Eg_max ){
            i = i - 1;
            cout<<"Eg threshold is higher than Eg_max. Will skip this event"<<endl;
            cout<<"P_Fermi = "<<p_prot_Fermi<<"   costThFermi = "<<cosThFermi<<"   Eg_Threshold = "<<Eg_thr<<endl;
            continue;
        }

        double psf_Eg = Eg_max - Eg_min;

        Eg = rand.Uniform(Eg_min, Eg_min + psf_Eg);
        flux_factor = N_EPA(Eb, Eg, q2_cut) + N_Brem(Eg, Eb);

        double s = Mp * Mp + 2 * Eg*(EFermi - p_prot_Fermi*cosThFermi );
        double t_min = T_min(0., Mp*Mp, Q2, Mp*Mp, s);
        double t_max = T_max(0., Mp*Mp, Q2, Mp*Mp, s);
        double psf_t = t_min - TMath::Max(t_max, t_lim);

        if (t_min > t_lim) {
            t = rand.Uniform(t_min - psf_t, t_min);
            
            // The x-sec is obtained in the frame where the proton is at rest, so we should move to that
            // frame, get the value of Eg, and estimate the cross-section.
            L_ProtFermi.SetPxPyPzE(pxFermi, pyFermi, pzFermi, EFermi);
            L_gamma.SetPxPyPzE(0, 0, Eg, Eg);
            L_gamma.Boost( -L_ProtFermi.BoostVector() );
            double Eg_ProtRestFrame = L_gamma.E();
            f_JPsi_dsigm_dt->SetParameter(0, Eg_ProtRestFrame);
            crs_JPsi = f_JPsi_dsigm_dt->Eval(t);

            double u = 2 * Mp * Mp + Q2 - s - t;
            double th_qprime = acos((s * (t - u) - Mp * Mp * (Q2 - Mp * Mp)) / sqrt(Lambda(s, 0, Mp * Mp) * Lambda(s, Q2, Mp * Mp))); //Byukling Kayanti (4.9)
            double th_pprime = PI + th_qprime;

            double Pprime = 0.5 * sqrt(Lambda(s, Q2, Mp * Mp) / s); // Momentum in c.m. it is the same for q_pr and p_pr

            // ** The LorentzVector of CM frame is equal L_gamma + L_proton_Fermi
            Lcm.SetPxPyPzE(pxFermi, pyFermi, pzFermi + Eg, EFermi + Eg);
            L_prot.SetPxPyPzE(Pprime * sin(th_pprime), 0., Pprime * cos(th_pprime), sqrt(Pprime * Pprime + Mp * Mp));
            L_gprime.SetPxPyPzE(Pprime * sin(th_qprime), 0., Pprime * cos(th_qprime), sqrt(Pprime * Pprime + Q2));

            double psf_cos_th = 2.; // cos(th):(-1 : 1)
            double psf_phi_cm = 2 * PI;

            double cos_th = rand.Uniform(-1., -1 + psf_cos_th);
            double sin_th = sqrt(1 - cos_th * cos_th);
            double phi_cm = rand.Uniform(0., 0. + psf_phi_cm);

            double El = sqrt(Q2) / 2.; // Energy of lepton in the rest frame of qprime
            double Pl = sqrt(El * El - m_l * m_l);

            L_lm.SetPxPyPzE(Pl * sin_th * cos(phi_cm), Pl * sin_th * sin(phi_cm), Pl*cos_th, El);
            L_lp.SetPxPyPzE(-Pl * sin_th * cos(phi_cm), -Pl * sin_th * sin(phi_cm), -Pl*cos_th, El);

            L_lm.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame
            L_lp.RotateY(th_qprime); // Rotate in order to get Z axis be antiparallel to the p_prime direction in the CM frame

            L_lm.Boost(L_gprime.BoostVector()); // Move to the CM Frame
            L_lp.Boost(L_gprime.BoostVector()); // Move to the CM Frame

            L_lm.Boost(Lcm.BoostVector()); // Move to the Lab Frame
            L_lp.Boost(Lcm.BoostVector()); // Move to the Lab Frame

            L_gprime.Boost(Lcm.BoostVector());
            L_prot.Boost(Lcm.BoostVector());

            double psf_phi_lab = 2 * PI;
            double phi_rot = rand.Uniform(0., psf_phi_lab);

            L_prot.RotateZ(phi_rot);
            L_gprime.RotateZ(phi_rot);
            L_lm.RotateZ(phi_rot);
            L_lp.RotateZ(phi_rot);
                        
            psf = psf_t*psf_Eg;

            double eta = Q2 / (2 * (s - Mp * Mp) - Q2);

            double vz = rand.Uniform(vz_min, vz_max);

            px_prot = L_prot.Px();
            py_prot = L_prot.Py();
            pz_prot = L_prot.Pz();
            E_prot = L_prot.E();
            px_lp = L_lp.Px();
            py_lp = L_lp.Py();
            pz_lp = L_lp.Pz();
            E_lp = L_lp.E();
            px_lm = L_lm.Px();
            py_lm = L_lm.Py();
            pz_lm = L_lm.Pz();
            E_lm = L_lm.E();


            double tot_weight = crs_JPsi * psf*flux_factor;

            if (write_root) {
                tr1->Fill();
            } else {
                // Writing Header
                Lund_out << 3 << setw(15) << 1 << setw(5) << 1 << setw(15) << psf << setw(15) << crs_JPsi << setw(15) << 0 << setw(15) << flux_factor << setw(15) << crs_JPsi << setw(15) << psf << setw(15) << tot_weight << endl;
                // Writing Proton
                Lund_out << 1 << setw(5) << 1 << setw(5) << 1 << setw(7) << 2212 << setw(5) << 0 << setw(5) << 0 << setw(15) << px_prot << setw(15) << py_prot << setw(15) << pz_prot;
                Lund_out << setw(15) << L_prot.E() << setw(15) << Mp << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
                // Writing Electron
                Lund_out << 2 << setw(5) << -1 << setw(5) << 1 << setw(7) << pid_neg_lep << setw(5) << 0 << setw(5) << 0 << setw(15) << px_lm << setw(15) << py_lm << setw(15) << pz_lm;
                Lund_out << setw(15) << L_lm.E() << setw(15) << m_l << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
                // Writing Positron
                Lund_out << 3 << setw(5) << 1 << setw(5) << 1 << setw(7) << pid_pos_lep << setw(5) << 0 << setw(5) << 0 << setw(15) << px_lp << setw(15) << py_lp << setw(15) << pz_lp;
                Lund_out << setw(15) << L_lp.E() << setw(15) << m_l << setw(15) << 0. << setw(15) << 0. << setw(15) << vz << endl;
            }
        } else {
            cout << " |t_min| > |t_lim|" << endl;
            cout << " t_min =  " << t_min << "   t_lim = " << t_lim << "  Eg = " << Eg << endl;
        }

    }



    if (write_root) {
        tr1->Write();
        h_P_Fermi1->Write();
        file_out->Close();
    }
}
