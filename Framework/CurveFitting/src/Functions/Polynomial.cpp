#include "MantidCurveFitting/Functions/Polynomial.h"
#include "MantidAPI/FunctionFactory.h"
#include <boost/lexical_cast.hpp>

using namespace Mantid::Kernel;

using namespace Mantid::API;

namespace Mantid {
namespace CurveFitting {
namespace Functions {

using namespace CurveFitting;

DECLARE_FUNCTION(Polynomial)

/** Constructor
 */
Polynomial::Polynomial() : m_n(0) {}

void Polynomial::init() {
  declareParameter("A0", 0, "Coefficient for the constant term");
}

/** Evaluate the function at the supplied points
 *  @param out The y-values (output)
 *  @param xValues The points to evaluate at
 *  @param nData The number of points
 */
void Polynomial::function1D(double *out, const double *xValues,
                            const size_t nData) const {
  // 1. Use a vector for all coefficient
  std::vector<double> coeff(m_n + 1, 0.0);
  for (int i = 0; i < m_n + 1; ++i)
    coeff[i] = getParameter(i);

  // 2. Calculate
  for (size_t i = 0; i < nData; ++i) {
    double x = xValues[i];
    double temp = coeff[0];
    double nx = x;
    for (int j = 1; j <= m_n; ++j) {
      temp += coeff[j] * nx;
      nx *= x;
    }
    out[i] = temp;
  }
}

/**
* Evaluate the Jacobian at the supplied points.
* @param out The Jacobian (output)
* @param xValues The points to evaluate at
* @param nData The number of points.
*/
void Polynomial::functionDeriv1D(API::Jacobian *out, const double *xValues,
                                 const size_t nData) {
  for (size_t i = 0; i < nData; i++) {
    double x = xValues[i];
    double nx = 1;
    for (int j = 0; j <= m_n; ++j) {
      out->set(i, j, nx);
      nx *= x;
    }
  }
}

/** Get Attribute names
 * @return A list of attribute names (identical to Polynomial)
*/
std::vector<std::string> Polynomial::getAttributeNames() const { return {"n"}; }

/** Get Attribute
 * @param attName :: Attribute name. If it is not "n" exception is thrown.
 * @return a value of attribute attName
 * (identical to Polynomial)
 */
API::IFunction::Attribute
Polynomial::getAttribute(const std::string &attName) const {
  if (attName == "n") {
    return Attribute(m_n);
  }

  throw std::invalid_argument("Polynomial: Unknown attribute " + attName);
}

/** Set Attribute
 * @param attName :: The attribute name. If it is not "n" exception is thrown.
 * @param att :: An int attribute containing the new value. The value cannot be
 * negative.
 * (identical to Polynomial)
 */
void Polynomial::setAttribute(const std::string &attName,
                              const API::IFunction::Attribute &att) {
  if (attName == "n") {
    // set the polynomial order
    if (m_n >= 0) {
      clearAllParameters();
    }
    m_n = att.asInt();
    if (m_n < 0) {
      throw std::invalid_argument(
          "Polynomial: polynomial order cannot be negative.");
    }
    for (int i = 0; i <= m_n; ++i) {
      std::string parName = "A" + std::to_string(i);
      declareParameter(parName);
    }
  }
}

/** Check if function has a given attribute
 *  @param attName The attribute name to check for
 *  @return True if the function has that attribute
 */
bool Polynomial::hasAttribute(const std::string &attName) const {
  return attName == "n";
}

} // namespace Functions
} // namespace CurveFitting
} // namespace Mantid
