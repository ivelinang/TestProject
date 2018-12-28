#include <math.h>
#include <iostream>
#include <assert.h>
#include <algorithm> 

const double almost_zero = 0.0000001;

bool is_almost_zero(double x)
{
	return std::abs(x) < almost_zero;
}

double x_of_z(double rho, double nu, double z)
{
	//x_of_z is unstable when nu is close to zero
	//x_of_z series expansion to 1st order in nu if nu is small
	const double nu_z = nu * z;
	if (std::abs(nu) < 0.00001)
		return (z + 0.5*z*nu_z*rho);
	else
		return(std::log((std::sqrt(1.0 - 2.0*rho*nu_z + nu_z * nu_z) - rho + nu_z) / (1.0 - rho)) / nu);

}

const double epsilonATM_SABR_SkewDamping = 0.01;
const double epsilonATM_SABR = 0.0005;
const double beta_threshold = 0.999;


double get_implied_ln_vol(double expiry, double moneyness, double alpha, double beta, double rho, double nu)
{
	if (std::abs(beta - 1.0) < almost_zero && is_almost_zero(nu))
		return alpha;

	const double u = moneyness - 1.0;
	const double one_sub_beta = 1.0 - beta;
	const double nu_squared = nu * nu;
	const double rho_nu = rho * nu;
	const double fudge_factor = 2.0;
	const double log_fudge_factor = std::log(fudge_factor);
	const double log_moneyness = log(moneyness);

	if (std::abs(u) < epsilonATM_SABR) //non ATM formula blows up if F and K close but not equal
	{
		if (is_almost_zero(one_sub_beta))
		{
			const double z_con = -log_fudge_factor / alpha;
			const double xofz_con = x_of_z(rho, nu, z_con);

			const double theta2 = std::log(xofz_con*sqrt(sqrt(1.0 - z_con * (2.0*rho_nu - nu_squared * z_con))) / z_con);
			const double corrfactnum = alpha * alpha * (-3.0 * alpha * (2.0 - moneyness) + 3.0*(4.0 - 2.0*moneyness)*rho_nu);
			const double corrfactden = 24.0 * (alpha * (2.0 - moneyness) + (1.0 - moneyness)*rho_nu);
			const double corrfact = corrfactnum / corrfactden + theta2 / (xofz_con*xofz_con);

			const double corrfact2 = 1.0 - 2.0*expiry*corrfact;
			const double multiplier = (1.0 - log_moneyness * rho_nu / (2.0*alpha));
			return (alpha / (multiplier * sqrt(corrfact2)));

		}
		else
		{
			const double fudge_factor_pow = std::exp(0.5*one_sub_beta * log_fudge_factor);
			const double z_con = (1.0 - fudge_factor_pow * fudge_factor_pow) / (alpha * one_sub_beta);
			const double zo_con = -fudge_factor_pow * log_fudge_factor / alpha;
			const double xofz_con = x_of_z(rho, nu, z_con);
			const double zofzo_on = x_of_z(rho, nu, zo_con);
			const double theta2 = std::log(xofz_con* (sqrt(sqrt(1.0 - z_con * (2.0 * rho_nu - nu_squared * z_con))) / z_con));

			const double corrfactnum = alpha * alpha * (1.0* alpha * (-2.0 + (-2.0 + beta)*beta)*(2.0 - moneyness) + 3.0*beta *((3.0 + beta) - (1.0 + beta)*(moneyness)*rho_nu));
			const double corrfactden = 24.0 *(alpha * (1.0 + beta - beta * moneyness) + (1.0 - moneyness)*rho_nu);
			const double corrfact = corrfactnum / corrfactden + theta2 / (xofz_con * xofz_con);
			const double corrfact2 = 1.0 - 2.0* expiry * corrfact;
			const double multiplier = (1.0 + log_moneyness * (alpha * one_sub_beta - rho_nu) / (2.0*alpha));
			return (alpha / (multiplier * sqrt(corrfact2)));
		}
	}
	else
	{
		if (is_almost_zero(one_sub_beta))
		{
			const double z = -log_moneyness / alpha;
			const double z_con = -log_fudge_factor / alpha;
			const double xofz = x_of_z(rho, nu, z);
			const double xofz_con = x_of_z(rho, nu, z_con);

			const double log_theta1_over_psi = 0.5*log_moneyness + log_fudge_factor - log(1.0 + moneyness);
			const double theta2 = log(xofz_con*sqrt(sqrt(1.0 - z_con * (2.0 - rho_nu - nu_squared * z_con))) / z_con);
			const double theta3 = 0.25 * rho_nu * alpha * z *z;

			const double temp = 1.0 - 2.0 * expiry * ((log_theta1_over_psi + theta3) / (xofz * xofz) + theta2 / (xofz_con*xofz));
			assert(temp > 0.0, "negative value"); //should return nan
			return (-log_moneyness / (xofz * sqrt(temp)));

		}
		else
		{
			const double fudge_factor_pow = exp(0.5*one_sub_beta * log_fudge_factor);
			const double moneyness_to_minus_half_beta = exp(-0.5*beta*log_moneyness);
			const double z = (1.0 - moneyness * moneyness_to_minus_half_beta*moneyness_to_minus_half_beta) / (alpha*one_sub_beta);
			const double z_con = (1.0 - fudge_factor_pow * fudge_factor_pow) / (alpha * one_sub_beta);
			const double zo = -log_moneyness * moneyness_to_minus_half_beta*sqrt(moneyness) / alpha;
			const double zo_con = -fudge_factor_pow * log_fudge_factor / alpha;
			const double xofz = x_of_z(rho, nu, z);
			const double xofzo = x_of_z(rho, nu, zo);
			const double xofz_con = x_of_z(rho, nu, z_con);
			const double xofzo_con = x_of_z(rho, nu, zo_con);

			const double psi = -0.5*log_moneyness*(1.0 + moneyness);
			const double theta1 = alpha * z / moneyness_to_minus_half_beta;
			const double theta2 = log(xofz_con*sqrt(sqrt(1.0 - z_con * (2.0*rho_nu - nu_squared * z_con))) / z_con);
			const double theta3 = 0.25 * rho_nu * alpha * beta * std::pow(0.5*(1.0 + moneyness), -one_sub_beta)*z*z;
			double temp = 1.0 - 2.0*expiry*((log(theta1 / psi) + theta3) / (xofzo*xofzo) + theta2 / (xofzo_con * xofzo_con));
			assert(temp > 0.0, "negative value");
			
			return (-log_moneyness / (xofz * sqrt(temp)));
		}
	}
	
}

double get_implied_ln_vol_skew_dampening(double expiry, double strike, double forward, double alpha, double beta, double rho, double nu)
{
	const double moneyness = strike / forward;
	const double one_sub_beta = 1.0 - beta;
	const double alpha_f_pow_b_m_one = alpha * pow(forward, beta - 1.0);
	assert(is_almost_zero(one_sub_beta), "beta = 1 not handled by SkewDamping yet");
	const double nu_squared = nu * nu;
	const double rho_nu = rho * nu;
	const double dampingLmabda = std::max(0.1*rho_nu*expiry, 0.0);
	const double u = moneyness - 1.0;
	const double log_fudge_factor = log(2.0);

	if (is_almost_zero(u))
	{
		const double theta1 = beta * (beta - 2.0) / 24.0;
		const double one_over_alpha = 1.0 / alpha_f_pow_b_m_one;
		const double theta3 = 0.25*beta*rho_nu*one_over_alpha;

		const double psi = 1.0 / 12.0;
		const double c_to_halfonesubbeta = exp(0.5*one_sub_beta*log_fudge_factor);
		const double zbar = (1.0 - c_to_halfonesubbeta * c_to_halfonesubbeta) / (alpha_f_pow_b_m_one*one_sub_beta);
		const double z0bar = -c_to_halfonesubbeta * log_fudge_factor * one_over_alpha;
		const double theta2 = log((x_of_z(rho, nu, zbar) * sqrt(sqrt(1.0 - zbar * (2.0*rho_nu - zbar * nu_squared)))) / zbar);
		const double x_of_z_z0bar = x_of_z(rho, nu, z0bar);
		const double term_under_root_sign = 1.0 - 2.0 *expiry*((theta1 + theta3 - psi)* alpha_f_pow_b_m_one * alpha_f_pow_b_m_one + theta2 / (x_of_z_z0bar*x_of_z_z0bar));
		const double theVol = alpha_f_pow_b_m_one / std::sqrt(term_under_root_sign);
		return theVol;

	}
	else if (std::abs(u) < epsilonATM_SABR_SkewDamping)
	{
		const double u_squared = u * u;
		const double u_cubed = u * u_squared;
		const double u_fourthed = u * u_cubed;

		double theta1 = 0.0;
		double theta3 = 0.0;
		{
			const double a1 = -0.5 *beta;
			const double a2 = -a1 * (1.0 + beta) / 3.0;
			const double a3 = -a2 * (2.0 + beta) / 4.0;
			const double a4 = -a3 * (3.0 + beta) / 5.0;
			const double a5 = -a4 * (4.0 + beta) / 6.0;
			const double a6 = -a5 * (5.0 + beta) / 7.0;

			const double b1 = -a1;
			const double b2 = b1 * (beta - 2.0) / 4.0;
			const double b3 = b2 * (beta - 4.0) / 6.0;
			const double b4 = b3 * (beta - 6.0) / 8.0;
			const double b5 = b4 * (beta - 8.0) / 10.0;
			const double b6 = b5 * (beta - 10.0) / 12.0;

			const double C = (a2 + a1 * b1 + b2) + (a3 + a2 * b1 + a1 * b2 + b3)*u + (a4 + a3 * b1 + a2 * b2 + a1 * b3 + b4)*u_squared + (b5 + a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1 + a5)*u_cubed + (b6 + a1 * b5 + a2 * b4 + a3 * b3 + a4 * b2 + a5 * b1 + a6)*u_fourthed;
			theta1 = C - 0.5*C*C*u_squared + C * C*C*u_fourthed / 3.0;

			const double q = -(1.0 / alpha_f_pow_b_m_one)*(1.0 + a1 * u + a2 * u_squared + a3 * u_cubed + a4 * u_fourthed);
			const double q_squared = q * q;
			const double A = q_squared * (1.0 - 0.25*dampingLmabda*q_squared*u_squared + 0.125*dampingLmabda*dampingLmabda*q_squared*q_squared*u_fourthed);
			theta3 = 0.25*rho_nu*alpha_f_pow_b_m_one*beta*A*pow(0.5*(1.0 + moneyness), -one_sub_beta);
		}
		const double rho_sq = rho * rho;
		double xz02 = 0.0;
		{

			const double A = -pow(forward*strike, 0.5*one_sub_beta) / alpha;
			const double z0 = A * (1.0 - 0.5*u + u_squared / 3.0 - 0.25*u_cubed + 0.2*u_fourthed);
			const double z0_sq = z0 * z0;
			const double z0_fr = z0_sq*z0_sq;
			xz02 = (z0 + rho_nu * 0.5*(z0_sq*u) + 0.5*(nu_squared)*(rho_sq - (1.0 / 3.0))*(z0*z0_sq*u_squared) + 0.125*(5.0*rho_sq - 3.0)*rho_nu*nu_squared*(z0_fr*u_cubed) + (1.0 / 40.0)*(nu_squared*nu_squared)*
					(3.0 - 30.0*rho_sq + 35.0*rho_sq*rho_sq)*(z0*z0_fr*u_fourthed));
			xz02 *= xz02;

		}
		double xzed = 0.0;
		{

			const double A = 1.0 / alpha_f_pow_b_m_one;
			const double z = A * (-1.0 + 0.5*beta*u - (1.0 / 6.0)*beta*(beta + 1.0)*u_squared + (1.0 / 24.0)*(beta*(beta + 1.0)*(beta + 2.0))*u_cubed - (1.0 / 120.0)*(beta*(beta + 1.0)*(beta + 2.0)*(beta + 3.0))*u_fourthed);
			const double z_sq = z * z;
			const double z_fr = z_sq * z_sq;
			//const double z_fr = z_sq * z_sq;
			xzed = (z + rho_nu * 0.5*(z_sq + u) + 0.5*(nu_squared)*(rho_sq - (1.0 / 3.0))*(z*z_sq*u_squared) + 0.125*(5.0*rho_sq - 3.0)*rho_nu*nu_squared*(z_fr*u_cubed) + (1.0 / 40.0)*(nu_squared*nu_squared)*(3.0 - 30.0*rho_sq + 35.0*rho_sq*rho_sq)*(z*z_fr*u_fourthed));

		}

		const double psi = 1.0 / 12.0 - (1.0 / 12.0)*u + (103.0 / 1440.0)*u_squared - (43.0 / 720.0)*u_cubed + (9071.0 / 181440.0)*u_fourthed;
		const double c_to_halfonesubbeta = exp(0.5*one_sub_beta*log_fudge_factor);
		const double zbar = (1.0 - c_to_halfonesubbeta * c_to_halfonesubbeta) / (alpha_f_pow_b_m_one*one_sub_beta);
		const double z0bar = -c_to_halfonesubbeta * log_fudge_factor / alpha_f_pow_b_m_one;
		const double theta2 = log((x_of_z(rho, nu, zbar)*sqrt(sqrt(1.0 - zbar * (2.0*rho_nu - zbar * nu_squared)))) / zbar);
		const double x_of_z_z0bar = x_of_z(rho, nu, z0bar);
		const double term_under_root_sign = 1.0 - 2.0*expiry*((theta1 + theta3 - psi) / xz02 + theta2 / (x_of_z_z0bar*x_of_z_z0bar));
		const double log_foverk = -(1.0 - 0.5*u + u_squared / 3.0 - 0.25*u_cubed + 0.2*u_fourthed);
		const double theVol = log_foverk / (xzed*sqrt(term_under_root_sign));
		return theVol;
			

	}
	else
	{
		const double a_f_pow_one_sub_beta = alpha_f_pow_b_m_one * one_sub_beta;
		const double log_moneyness = log(moneyness);
		const double moneyness_to_minus_half_beta = exp(-0.5*beta*log_moneyness);
		const double z = (1.0 - moneyness * moneyness_to_minus_half_beta*moneyness_to_minus_half_beta) / a_f_pow_one_sub_beta;
		const double dampedZsq = is_almost_zero(dampingLmabda) ? (z*z) : ((2.0 / dampingLmabda)*(sqrt(1.0 + dampingLmabda * z*z) - 1.0));
		const double fudge_factor_pow = exp(0.5*one_sub_beta*log_fudge_factor);
		const double z_con = (1.0 - fudge_factor_pow * fudge_factor_pow) / a_f_pow_one_sub_beta;
		const double zo = -sqrt(moneyness)*moneyness_to_minus_half_beta*log_moneyness / alpha_f_pow_b_m_one;
		const double zo_con = -fudge_factor_pow * log_fudge_factor / alpha_f_pow_b_m_one;
		const double xofz = x_of_z(rho, nu, z);
		const double xofzo = x_of_z(rho, nu, zo);
		const double xofz_con = x_of_z(rho, nu, z_con);
		const double xofzo_con = x_of_z(rho, nu, zo_con);

		const double f_m_k = forward - strike;
		const double psi = log(-0.5*log_moneyness*(forward + strike) / f_m_k);
		const double theta1 = log(alpha*z*pow(forward*strike, 0.5*beta) / f_m_k);
		const double theta2 = log(xofz_con*sqrt(sqrt(1.0 - z_con * (2.0*rho_nu - nu_squared * z_con))) / z_con);
		const double theta3 = 0.25*rho_nu*alpha_f_pow_b_m_one*beta*pow(0.5*(1.0 + moneyness), -one_sub_beta)*dampedZsq;
		const double temp = 1.0 - 2.0*expiry*((theta1 + theta3 - psi) / (xofzo*xofzo) + theta2 / (xofzo_con*xofzo_con));
		assert(temp > 0.0, "negative value");
		return (-log_moneyness / (xofz*sqrt(temp)));
	}
}

double option_price(double expiry,
	double discount,
	double forward,
	double strike,
	bool isCall,
	double alpha,
	double beta,
	double rho,
	double nu,
	bool is_skew_dampening)
{
	assert(std::abs(rho) < 1.0, "rho>1.0 not supported");

	if (forward < 0.0)
	{
		forward *= -1;
		strike *= -1;
		isCall = !isCall;
	}

	//for beta close to one we go for beta=1.0 which is stable. Beta not being fitted 1e-3 is largely sufficient
	if (beta > beta_threshold)
		beta = 1.0;
	double ln_vol;
	if (strike <= 0.0 || alpha < 1e-10 || is_almost_zero(forward))
	{
		ln_vol = 0.0;
	}
	else if (is_skew_dampening)
	{
		ln_vol = get_implied_ln_vol_skew_dampening(expiry, strike, forward, alpha, beta, rho, nu);
	}
	else
	{
		ln_vol = get_implied_ln_vol(expiry, strike / forward, alpha*pow(forward, beta - 1.0), beta, rho, nu);
	}
	return 1.0; //log normal model price TODO
}

const double epsilonATM_SABR_Normal = 0.00001;

double get_implied_normal_vol(double expiry, double moneyness, double alpha, double rho, double nu)
{
	assert(std::abs(rho)<1.0, "rho>=1 not supported");
	if (std::abs(alpha) < 1e-10)
		return 0;
	if (std::abs(nu) < 1e-5)
		return alpha;

	const double z = nu * moneyness / alpha;
	//for around the money case
	if (std::abs(moneyness) < epsilonATM_SABR_Normal)
	{
		//we go up to 3rd order to avoid second derivatives by FD to be 0

		const double normal_vol = alpha * (1.0 + (2.0 - 3.0*rho*rho) / 24.0*nu*nu*expiry)*(1.0 - rho * z / 2.0 + (2.0 - 3.0*rho*rho)*z*z / 12.0 + (5.0*rho - 6.0*rho*rho*rho)*z*z*z / 24.0);
		return normal_vol;
	}
	else
	{
		const double expr = 1.0 - z * (2.0*rho - z);

	}
}

