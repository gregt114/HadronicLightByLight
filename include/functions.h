
#include <cmath>

typedef long double double16; // Use 16 byte(128 bit) precision

const double16 pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062L;

const double16 m_mu = 0.1056583745;   // GeV / c^2 (muon mass)
const double16 mpion = 0.1349768;     // GeV / c^2 (neutral pion mass)

const double16 F_pi = 0.0924;         // GeV
const double16 Mv1 = 0.77549;         // GeV
const double16 Mv2 = 1.465;           // GeV
const double16 h1 = 0.0;

const double16 alpha = 0.0072973525693; // Fine Structure Constant

// Form Factor Constants (LMD+V)
const double16 deltasq = 0.2;         // GeV^2
const double16 h2 = -4.0*(pow(Mv1,2) + pow(Mv2,2)) + deltasq*(16.0/9);
const double16 h5 = 6.93;             // GeV^4
const double16 h7 = -3.0*pow(Mv1*Mv2,4) / (4*pow(pi*F_pi,2));   // GeV^6

// Form Factor Expansion Parameters
const double16 gamma_pi0 = 7.7291993365804157745e-09;            // GeV TODO
const double16 a0 = 1/pow(Mv1,2) + 1/pow(Mv2,2) + h5/h7;
const double16 b0 = 1/pow(Mv1,4) + 1/pow(Mv2,4) + 1/pow(Mv1*Mv2,2) + h5/h7*(1/pow(Mv1,2) + 1/pow(Mv2,2));
const double16 c0 = pow(pow(Mv1,-2) + pow(Mv2,-2), 2) + h2/h7 + 2*h5/h7 * (pow(Mv1,-2) + pow(Mv2,-2));
const double16 d0 = -4.59258;    // GeV^-6
const double16 e0 = -5.58268;    // GeV^-6


double16 Rmi(double16 q)
{
    return sqrt(1 + 4*pow(m_mu,2) / pow(q,2));
}


double16 X(double16 q1, double16 q2, double16 tau)
{
    double16 x = sqrt(1 - pow(tau, 2));

    double16 z = q1*q2/(4*pow(m_mu,2)) * (1 - Rmi(q1))*(1-Rmi(q2));

    return 1/(q1*q2*x) * atan(z*x / (1-z*tau));
}

// Used in calculating weighting functions (I1 in Nyffeler Paper)
double16 eye1(double16 q1, double16 q2, double16 tau)
{
    double16 q1sq = pow(q1,2);
    double16 q2sq = pow(q2,2);
    double16 dot = q1*q2*tau;
    double16 q3sq = q1sq + 2*dot + q2sq;

    double16 p1 = 1/q1sq;
    double16 p2 = 1/q2sq;
    double16 p3 = 1/q3sq;
    double16 m_musq = pow(m_mu,2);

    double16 result = X(q1, q2, tau);
    
    double16 factor = 8*p1*p2*dot;
    factor += -2*p1*p3*(pow(q2,4)/m_musq - 2*q2sq);
    factor += 4*p2*p3*q1sq - 4*p2;
    factor += -2*p1*(2 - q2sq/m_musq + 2*(dot/m_musq));
    factor += -2*p3*(4 + q1sq/m_musq -2*q2sq/m_musq);
    factor += 2/m_musq;
    
    result *= factor;
    result += -2*p1*p2*(1 + (1-Rmi(q1))*dot/m_musq);
    result += p1*p3*(2 - (1-Rmi(q1))*q2sq/m_musq);
    result += p2*p3*(2 + pow(1-Rmi(q1),2)*dot/m_musq);
    result += p1*(1-Rmi(q1))/m_musq;
    result += 3*p3*(1-Rmi(q1))/m_musq;

    return result;
}

// Used in calculating weighting functions (I2 in Nyffeler Paper)
double16 eye2(double16 q1, double16 q2, double16 tau)
{
    double16 q1sq = pow(q1,2);
    double16 q2sq = pow(q2,2);
    double16 dot = q1*q2*tau;
    double16 q3sq = q1sq + 2*dot + q2sq;

    double16 p1 = 1/q1sq;
    double16 p2 = 1/q2sq;
    double16 p3 = 1/q3sq;
    double16 m_musq = pow(m_mu,2);

    double16 result = X(q1, q2, tau);

    double16 factor = 4*p1*p2*dot;
    factor += 2*p1*p3*q2sq;
    factor += -2*p1 + 2*p2*p3*q1sq - 2*p2 - 4*p3 - 4/m_musq;
    result *= factor;
    
    result += -2*p1*p2 - 3*p1*(1-Rmi(q2))/(2*m_musq);
    result += -3*p2*(1-Rmi(q1))/(2*m_musq);
    result += -p3*(2 - Rmi(q1) - Rmi(q2))/(2*m_musq);
    result += p1*p3*(2 + 3*(1-Rmi(q2))*q2sq/(2*m_musq) + pow(1-Rmi(q2),2)*dot/(2*m_musq));
    result += p2*p3*(2 + 3*(1-Rmi(q1))*q1sq/(2*m_musq) + pow(1-Rmi(q1),2)*dot/(2*m_musq));

    return result;
}


// 1st weighting function
double16 w1(double16 q1, double16 q2, double16 tau)
{
    double16 result = (-2*pi/3) * sqrt(1-pow(tau,2));
    result *= pow(q1*q2,3) / (pow(q2,2) + pow(mpion,2)) * eye1(q1,q2,tau);
    return result;
}

// 2nd weighting function
double16 w2(double16 q1, double16 q2, double16 tau)
{
    double16 q3sq = pow(q1,2) + 2*q1*q2*tau + pow(q2,2);
    
    double16 result = (-2*pi/3) * sqrt(1-pow(tau,2));
    result *= pow(q1*q2,3) / (q3sq + pow(mpion,2)) * eye2(q1,q2,tau);   
    return result;
}


// LMD+V Form Factor
double16 form_factor(double16 q1sq, double16 q2sq)
{
    // Note that this function is F(q1^2, q2^2) rather than F(q1, q2)

    double16 numerator = q1sq*q2sq*(q1sq+q2sq);
    numerator += -h2*q1sq*q2sq + h5*(q1sq+q2sq) - h7;
    double16 denom = (q1sq+pow(Mv1,2))*(q1sq+pow(Mv2,2));
    denom *= (q2sq+pow(Mv1,2))*(q2sq+pow(Mv2,2));

    return (F_pi/3.0) * numerator / denom;
}

// Integrand1 using full form factor
double16 integrand1(double16 q1, double16 q2, double16 tau)
{
    return w1(q1,q2,tau)*form_factor(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau))*form_factor(pow(q2,2),0.0);
}

// Integrand2 using full form factor
double16 integrand2(double16 q1, double16 q2, double16 tau)
{
    return w2(q1,q2,tau)*form_factor(pow(q1,2), pow(q2,2))*form_factor((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0);
}

////////////////////////////////////////////////////////////////////////////////////

// Q^6 approximation of form factor
double16 form_factor_q6(double16 q1sq, double16 q2sq, double16 a=a0, double16 b=b0, double16 c=c0, double16 d=d0, double16 e=e0, double16 gamma=gamma_pi0)
{

    double16 factor = sqrt(4*gamma / (pi * pow(alpha,2) * pow(mpion,3)));
    double16 result = factor * (1 - a*(q1sq + q2sq) + b*(pow(q1sq,2) + pow(q2sq,2)) + c*q1sq*q2sq + d*(pow(q1sq,3)+pow(q2sq,3)) + e*(pow(q1sq,4)*q2sq + q1sq*pow(q2sq,2)));
    return result;
}

// Q^4 approximation of form factor
double16 form_factor_q4(double16 q1sq, double16 q2sq, double16 a=a0, double16 b=b0, double16 c=c0)
{

    double16 factor = sqrt(4*gamma_pi0 / (pi * pow(alpha,2) * pow(mpion,3)));
    double16 result = factor * (1 - a*(q1sq + q2sq) + b*(pow(q1sq,2) + pow(q2sq,2)) + c*q1sq*q2sq);
    return result;
}

// Integrand1 using Q^4 approximation
double16 integrand1_q4(double16 q1, double16 q2, double16 tau, double16 a=a0, double16 b=b0, double16 c=c0)
{
    return w1(q1,q2,tau)*form_factor_q4(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau), a, b, c)*form_factor_q4(pow(q2,2),0.0,a, b, c);
}

// Integrand2 using Q^4 approximation
double16 integrand2_q4(double16 q1, double16 q2, double16 tau, double16 a=a0, double16 b=b0, double16 c=c0)
{
    return w2(q1,q2,tau)*form_factor_q4(pow(q1,2), pow(q2,2), a, b, c)*form_factor_q4((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0, a, b, c);
}

// Integrand1 using Q^6 approximation
double16 integrand1_q6(double16 q1, double16 q2, double16 tau, double16 a=a0, double16 b=b0, double16 c=c0, double16 d=d0, double16 e=e0, double16 gamma=gamma_pi0)
{
    return w1(q1,q2,tau)*form_factor_q6(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau), a, b, c, d, e, gamma)*form_factor_q6(pow(q2,2),0.0, a, b, c, d, e, gamma);
}

// Integrand2 using Q^6 approximation
double16 integrand2_q6(double16 q1, double16 q2, double16 tau, double16 a=a0, double16 b=b0, double16 c=c0, double16 d=d0, double16 e=e0, double16 gamma=gamma_pi0)
{
    return w2(q1,q2,tau)*form_factor_q6(pow(q1,2), pow(q2,2), a, b, c, d, e, gamma)*form_factor_q6((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0, a, b, c, d, e, gamma);
}



