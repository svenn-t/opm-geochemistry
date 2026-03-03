/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <opm/simulators/geochemistry/Thermo/water.h>

water::water()
: v_(0.0)
, u_(0.0)
, s_(0.0)
, h_(0.0)
, cp_(0.0)
, cv_(0.0)
, w_(0.0)
, g_(0.0)
, G_(0.0)
, H_(0.0)
, denst_(0.0)
, alpha_(0.0)
, alpha_t_(0.0)
, beta_(0.0)
, Psat_(0.0)
//
, P_(0.0)
, T_(-273.15)
, R_(0.461526e3)  // specific gas constant of water [J/kg/K]
, Mw_(18.01528e-3)  // [kg/mol]
, Tcrit_(647.096)
, Pcrit_(22.064e6)
, rho_crit_(322.0)  // [kg/m^3]
{

}

/* Note: If T==T && P==P, the calculation is already done. Pressure is in Pascals, temperature in Kelvin. */
void water::gibbsIAPWS(double T, double P)
{
    if (T != T_ || P != P_){
        gibbsIAPWSlocal(T, P);
    }
}

void water::gibbsIAPWSlocal(double T, double P)
{
    T_ = T;
    P_ = P;
    if (T > 623.15 || T < 273.15)
    {
        std::cout << "Temperature: " << T << " outside range 273.15 to 623 K\n";
        exit(0);
    }
	Psat_ = PsatIAPWS(T);
	if (P < Psat_ || P> 100e6){
        std::cout << "Pressure:" << P << " is outside range from saturation pressure ";
        std::cout << Psat_ << " to 100MPa \n.";
		exit(0);
	}



    static constexpr std::array<int, 34> Ii = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2,
                                                   2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32 };

    static constexpr std::array<int, 34> Ji = { -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1,
                                                 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41 };

	static constexpr std::array<double, 34> ni =
    {
        1.46329712131670e-01, -8.45481871691140e-01, -3.75636036720400, 3.38551691683850, -9.57919633878720e-01,
        1.57720385132280e-01, -1.66164171995010e-02, 8.12146299835680e-04, 2.83190801238040e-04, -6.07063015658740e-04,
        -1.89900682184190e-02, -3.25297487705050e-02, -2.18417171754140e-02, -5.28383579699300e-05, -4.71843210732670e-04,
        -3.00017807930260e-04, 4.76613939069870e-05, -4.41418453308460e-06, -7.26949962975940e-16, -3.16796448450540e-05,
        -2.82707979853120e-06, -8.52051281201030e-10, -2.24252819080000e-06, -6.51712228956010e-07, -1.43417299379240e-13,
        -4.05169968601170e-07, -1.27343017416410e-09, -1.74248712306340e-10, -6.87621312955310e-19, 1.44783078285210e-20,
        2.63357816627950e-23, -1.19476226400710e-23, 1.82280945814040e-24, -9.35370872924580e-26
    };

	static constexpr double Ps = 16.53e6;  // 16.53 MPa
    static constexpr double Ts = 1386.0;  // K

    const double pi = P / Ps;
    const double pi_diff = 7.1 - pi;
    const double tau = Ts / T;
    const double tau_diff = tau - 1.222;

    double pi_pow = 0.0;
    double dpi_pow = 0.0;
    double ddpi_pow = 0.0;
    double tau_pow = 0.0;
    double dtau_pow = 0.0;
    double ddtau_pow = 0.0;

    double g = 0.0;  // dimensionless Gibbs free energy (in paper: lambda)
    // For the derivatives, we use shorter notation p~pi, and t~tau
    double dg_p = 0.0;
    double dg_pp = 0.0;
    double dg_t = 0.0;
    double dg_tt = 0.0;
    double dg_pt = 0.0;
    double dg_ptt = 0.0;

    // Equation (7) in IAPWS-97 paper ("the basic equation")
    for (std::size_t i=0; i < ni.size(); ++i)
    {
        //mypow(pi_diff, Ii[i], pi_pow, dpi_pow, ddpi_pow);
		//mypow(tau_diff, Ji[i], tau_pow, dtau_pow, ddtau_pow);
        nth_power(pi_diff, Ii[i], pi_pow, dpi_pow, ddpi_pow);
        nth_power(tau_diff, Ji[i], tau_pow, dtau_pow, ddtau_pow);

		g += ni[i] * pi_pow * tau_pow;
        dg_p += ni[i] * dpi_pow * tau_pow;
		dg_pp += ni[i] * ddpi_pow * tau_pow;
		dg_t += ni[i] * pi_pow * dtau_pow;
		dg_tt += ni[i] * pi_pow * ddtau_pow;
		dg_pt += ni[i] * dpi_pow * dtau_pow;
		dg_ptt += ni[i] * dpi_pow * ddtau_pow;
	}
    // Because the derivative of the kernel is -1 (chain rule, d/dpi)
    dg_p = -dg_p;
    dg_pt = -dg_pt;
    dg_ptt = -dg_ptt;

    // Table 3 of IAPWS-97 paper
	v_ = pi*dg_p*R_*T / P; // kg/m^3
	u_ = R_*T*(tau*dg_t - pi*dg_p);
	s_ = R_*(tau*dg_t - g);
	h_ = R_*T*tau*dg_t;
	cp_ = -tau*tau*dg_tt*R_;
	const double cvi = (dg_p - tau*dg_pt);
	cv_ = cp_ + R_*cvi *cvi/ dg_pp;
	w_ = R_*T*dg_p*dg_p / (cvi*cvi / tau / tau / dg_tt - dg_pp);
	w_ = sqrt(w_);

    g_ = g*R_*T;
    G_ = g_*Mw_;
    H_ = h_*Mw_;
    denst_ = 1.0 / v_;

    // Question: Where are these formulas from?
    // beta_  = water isothermal compressibility = -1/v_(dv_/dp)_T
    // alpha_ = water isobaric compressibility = 1/v_(dv_/dT)_p
    // calculated analytically by replacing specific volume with v_=RT/Ps gamma_pi in
    // Table 3 of IAPWS-97 paper
	alpha_ = (dg_p - tau*dg_pt) / T / dg_p;
	beta_ = -dg_pp / Ps / dg_p;
	alpha_t_ = (dg_p - tau*dg_pt);
	alpha_t_ *= alpha_t_;
	alpha_t_ = tau*tau*dg_ptt*dg_p - alpha_t_;
	alpha_t_ = alpha_t_ / (T*T*dg_p*dg_p);
}

/* Returns the saturation pressure for 273.15 K <= T <= 647.096 K. */
double water::PsatIAPWS(double T)
{
    if (T <273.15 || T > 647.096)
    {
        std::cout << "Temperature is outside valid range " << T << "\n";
        std::cout << "Only region 1 is implemented in GeoChemX, to be updated ... " << std::endl;
        exit(0);
}

    static constexpr std::array<double, 10> n = { 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
                                             -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
                                             -0.23855557567849, 0.65017534844798e3 };
    static constexpr double Ts = 1.0;

    const double t = T / Ts;
    const double Th = t + n[8] / (t - n[9]);
    const double Th2 = Th*Th;
    const double A = Th2 + n[0] * Th + n[1];
    const double B = n[2] * Th2 + n[3] * Th + n[4];
    const double C = n[5] * Th2 + n[6] * Th + n[7];

    double p = -B + sqrt(B*B - 4 * A*C);
    p = 2.0*C / p;

    return 1.0e6*p*p*p*p;
}

void water::printProperties() const
{
    printf("T[K]\tP[MPa]\tv[kg/m3]\th[kJ/kg]\tu[kJ/kg]\ts[kJ/kgK]\tcp[kJ/kgK]\tcv[kJ/kgK]\tw[m/s]\trho[kg/m3]\n");
    printf
    (
        //"{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\t{:4.8e}\n",
        "%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\t%4.8e\n",
        T_,
        1.0e-6*P_,
        v_,
        1.0e-3*h_,
        1.0e-3*u_,
        1.0e-3*s_,
        1.0e-3*cp_,
        1.0e-3*cv_,
        w_,
        denst_
    );
}

void water::mypow(double x, int n, double& xn, double& dxn, double& ddxn)
{
    if (n == 0)
    {
	xn = 1.0;
        dxn = 0.0;
        ddxn = 0.0;
    }
    else if (n == 1)
    {
        xn = x;
        dxn = 1.0;
        ddxn = 0.0;
    }
    else if (n == 2)
    {
        xn = x*x;
        dxn = 2.0*x;
        ddxn = 2.0;
    }
    else if (xn == 0.0)
    {
        xn = dxn = ddxn = 0.0;
    }
    else
    {
        const bool positive_power = (n>0);
        const int p = positive_power ? n : -n;
        double xp = 1.0;
        for (int i = 0; i < p; ++i)
        {
            xp *= x;
        }

        const auto dbl_n = static_cast<double>(n);
        xn = positive_power ? xp : 1.0 / xp;

        dxn = dbl_n*xn / x;
        ddxn = (dbl_n - 1.0)*dxn / x;
    }
}

