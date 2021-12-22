
#include <cmath>

const double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062L;

const double m_mu = 0.1056583745;   // GeV / c^2 (muon mass)
const double mpion = 0.1349768;     // GeV / c^2 (neutral pion mass)

const double F_pi = 0.0924;         // GeV
const double Mv1 = 0.77549;         // GeV
const double Mv2 = 1.465;           // GeV
const double h1 = 0.0;

const double alpha = 0.0072973525693; // Fine Structure Constant

// Form Factor Constants (LMD+V)
const double deltasq = 0.2;         // GeV^2
const double h2 = -4.0*(pow(Mv1,2) + pow(Mv2,2)) + deltasq*(16.0/9);
const double h5 = 6.93;             // GeV^4
const double h7 = -3.0*pow(Mv1*Mv2,4) / (4*pow(pi*F_pi,2));   // GeV^6

// Form Factor Expansion Parameters
const double gamma_pi0 = 7.7291993365804157745e-09;            // GeV TODO
const double a0 = 1/pow(Mv1,2) + 1/pow(Mv2,2) + h5/h7;
const double b0 = 1/pow(Mv1,4) + 1/pow(Mv2,4) + 1/pow(Mv1*Mv2,2) + h5/h7*(1/pow(Mv1,2) + 1/pow(Mv2,2));
const double c0 = pow(pow(Mv1,-2) + pow(Mv2,-2), 2) + h2/h7 + 2*h5/h7 * (pow(Mv1,-2) + pow(Mv2,-2));
const double d0 = -4.59258;    // GeV^-6
const double e0 = -5.58268;    // GeV^-6


double Rmi(double q)
{
    return sqrt(1 + 4*pow(m_mu,2) / pow(q,2));
}


double X(double q1, double q2, double tau)
{
    double x = sqrt(1 - pow(tau, 2));

    double z = q1*q2/(4*pow(m_mu,2)) * (1 - Rmi(q1))*(1-Rmi(q2));

    return 1/(q1*q2*x) * atan(z*x / (1-z*tau));
}

// Used in calculating weighting functions (I1 in Nyffeler Paper)
double eye1(double q1, double q2, double tau)
{
    double q1sq = pow(q1,2);
    double q2sq = pow(q2,2);
    double dot = q1*q2*tau;
    double q3sq = q1sq + 2*dot + q2sq;

    double p1 = 1/q1sq;
    double p2 = 1/q2sq;
    double p3 = 1/q3sq;
    double m_musq = pow(m_mu,2);

    double result = X(q1, q2, tau);
    
    double factor = 8*p1*p2*dot;
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
double eye2(double q1, double q2, double tau)
{
    double q1sq = pow(q1,2);
    double q2sq = pow(q2,2);
    double dot = q1*q2*tau;
    double q3sq = q1sq + 2*dot + q2sq;

    double p1 = 1/q1sq;
    double p2 = 1/q2sq;
    double p3 = 1/q3sq;
    double m_musq = pow(m_mu,2);

    double result = X(q1, q2, tau);

    double factor = 4*p1*p2*dot;
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
double w1(double q1, double q2, double tau)
{
    double result = (-2*pi/3) * sqrt(1-pow(tau,2));
    result *= pow(q1*q2,3) / (pow(q2,2) + pow(mpion,2)) * eye1(q1,q2,tau);
    return result;
}

// 2nd weighting function
double w2(double q1, double q2, double tau)
{
    double q3sq = pow(q1,2) + 2*q1*q2*tau + pow(q2,2);
    
    double result = (-2*pi/3) * sqrt(1-pow(tau,2));
    result *= pow(q1*q2,3) / (q3sq + pow(mpion,2)) * eye2(q1,q2,tau);   
    return result;
}


// LMD+V Form Factor
double form_factor(double q1sq, double q2sq)
{
    // Note that this function is F(q1^2, q2^2) rather than F(q1, q2)

    double numerator = q1sq*q2sq*(q1sq+q2sq);
    numerator += -h2*q1sq*q2sq + h5*(q1sq+q2sq) - h7;
    double denom = (q1sq+pow(Mv1,2))*(q1sq+pow(Mv2,2));
    denom *= (q2sq+pow(Mv1,2))*(q2sq+pow(Mv2,2));

    return (F_pi/3.0) * numerator / denom;
}

// Integrand1 using full form factor
double integrand1(double q1, double q2, double tau)
{
    return w1(q1,q2,tau)*form_factor(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau))*form_factor(pow(q2,2),0.0);
}

// Integrand2 using full form factor
double integrand2(double q1, double q2, double tau)
{
    return w2(q1,q2,tau)*form_factor(pow(q1,2), pow(q2,2))*form_factor((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0);
}

////////////////////////////////////////////////////////////////////////////////////

// Q^6 approximation of form factor
double form_factor_q6(double q1sq, double q2sq, double a=a0, double b=b0, double c=c0, double d=d0, double e=e0, double gamma=gamma_pi0)
{

    double factor = sqrt(4*gamma / (pi * pow(alpha,2) * pow(mpion,3)));
    double result = factor * (1 - a*(q1sq + q2sq) + b*(pow(q1sq,2) + pow(q2sq,2)) + c*q1sq*q2sq + d*(pow(q1sq,3)+pow(q2sq,3)) + e*(pow(q1sq,4)*q2sq + q1sq*pow(q2sq,2)));
    return result;
}

// Q^4 approximation of form factor
double form_factor_q4(double q1sq, double q2sq, double a=a0, double b=b0, double c=c0)
{

    double factor = sqrt(4*gamma_pi0 / (pi * pow(alpha,2) * pow(mpion,3)));
    double result = factor * (1 - a*(q1sq + q2sq) + b*(pow(q1sq,2) + pow(q2sq,2)) + c*q1sq*q2sq);
    return result;
}

// Integrand1 using Q^4 approximation
double integrand1_q4(double q1, double q2, double tau, double a=a0, double b=b0, double c=c0)
{
    return w1(q1,q2,tau)*form_factor_q4(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau), a, b, c)*form_factor_q4(pow(q2,2),0.0,a, b, c);
}

// Integrand2 using Q^4 approximation
double integrand2_q4(double q1, double q2, double tau, double a=a0, double b=b0, double c=c0)
{
    return w2(q1,q2,tau)*form_factor_q4(pow(q1,2), pow(q2,2), a, b, c)*form_factor_q4((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0, a, b, c);
}

// Integrand1 using Q^6 approximation
double integrand1_q6(double q1, double q2, double tau, double a=a0, double b=b0, double c=c0, double d=d0, double e=e0, double gamma=gamma_pi0)
{
    return w1(q1,q2,tau)*form_factor_q6(pow(q1,2), (pow(q1,2)+pow(q2,2)+2*q1*q2*tau), a, b, c, d, e, gamma)*form_factor_q6(pow(q2,2),0.0, a, b, c, d, e, gamma);
}

// Integrand2 using Q^6 approximation
double integrand2_q6(double q1, double q2, double tau, double a=a0, double b=b0, double c=c0, double d=d0, double e=e0, double gamma=gamma_pi0)
{
    return w2(q1,q2,tau)*form_factor_q6(pow(q1,2), pow(q2,2), a, b, c, d, e, gamma)*form_factor_q6((pow(q1,2)+pow(q2,2)+2*q1*q2*tau),0.0, a, b, c, d, e, gamma);
}



