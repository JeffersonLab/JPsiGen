/* 
 * File:   RadiativeCorrections.h
 * Author: pierrec 
 *
 * Created on January 11, 2020, 4:38 PM
 */

#ifndef RADIATIVECORRECTIONS_H
#define RADIATIVECORRECTIONS_H

#include <TLorentzVector.h>
#include <TF1.h>
#include <TRandom3.h>


/*namespace RadiativeFuncs {

    void Soft_Photon_Emission(TLorentzVector&, TLorentzVector&, TLorentzVector&, double , double);

}*/


class RadiativeCorrections {
public:

    RadiativeCorrections(double cut_off_min);

    static double R(double rm);
    static double E_to_rm(double E);
    static double VM_delta(double rm);
    static double VM_delta_fun(double *x, double *par);
    static double constant_VM_delta(double rm);
    static double d_VM_delta_fun(double *x, double *par);
    static double varying_VM_delta(double rm);
    static double varying_VM_delta_fun(double *x, double *par);

    void Set_Inv_Mass(double InvMass);
    double Compute_cs_correction_factor(double Inv_Mass);
    double Compute_cs_correction_factor_VM();
    void Soft_Photon_Emission(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&);
    void Soft_Photon_Emission_VM(TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&);
    virtual ~RadiativeCorrections();

private:

    double cut_off_min;
    double cut_off_max;
    double cs_correction_factor;
    TF1* cs_function;
    TF1* cs_correction_factor_func;
    TF1 *fun_varying_VM_delta;
    static constexpr double me = 0.00051;
    static constexpr double PI = 3.14159265358979312;
    static constexpr double Mv = 3.0969;


};

#endif /* RADIATIVECORRECTIONS_H */