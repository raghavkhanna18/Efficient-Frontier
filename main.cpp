#include <iostream>
#include <ql/quantlib.hpp>
#include <boost/format.hpp>
#include <functional>
#include <numeric>
#include <fstream>

using namespace QuantLib;

double calculatePortfolioReturn (double proportionA, double expectedReturnA,
                                 double expectedReturnB) {
  return proportionA * expectedReturnA +
         (1 - proportionA) * expectedReturnB;
}

Volatility
calculatePortfolioRisk (double proportionA, Volatility volatilityA,
                        Volatility volatilityB,
                        double covarianceAB) {
  return std::sqrt (std::pow (proportionA, 2) * std::pow (volatilityA, 2) +
                    std::pow (1 - proportionA, 2) *
                    std::pow (volatilityB, 2) +
                    (2 * proportionA * (1 - proportionA) * covarianceAB));
}

int main () {

  Matrix covarianceMatrix (4, 4);

//row 1
  covarianceMatrix[0][0] = .40;  //Equity1-Equity1
  covarianceMatrix[0][1] = .03;  //Equity1-Equity2
  covarianceMatrix[0][2] = .02;  //Equity1-Equity3
  covarianceMatrix[0][3] = .06;  //Equity1-Equity4
//row 2
  covarianceMatrix[1][0] = .03;  //Equity2-Equity1
  covarianceMatrix[1][1] = .20;  //Equity2-Equity2
  covarianceMatrix[1][2] = .01;  //Equity2-Equity3
  covarianceMatrix[1][3] = -.06; //Equity2-Equity4
//row 3
  covarianceMatrix[2][0] = .02;  //Equity3-Equity1
  covarianceMatrix[2][1] = .01;  //Equity3-Equity2
  covarianceMatrix[2][2] = .30;  //Equity3-Equity3
  covarianceMatrix[2][3] = .03;  //Equity3-Equity4
//row 4
  covarianceMatrix[3][0] = .06;  //Equity4-Equity1
  covarianceMatrix[3][1] = -.06; //Equity4-Equity2
  covarianceMatrix[3][2] = .03;  //Equity4-Equity3
  covarianceMatrix[3][3] = .15;  //Equity4-Equity4

  std::cout << "Covariance matrix of returns: " << std::endl;
  std::cout << covarianceMatrix << std::endl;

//portfolio return vector         
  Matrix portfolioReturnVector (4, 1);
  portfolioReturnVector[0][0] = .19; //Equity1
  portfolioReturnVector[1][0] = .11; //Equity2
  portfolioReturnVector[2][0] = .07; //Equity3
  portfolioReturnVector[3][0] = .08; //Equity4

  std::cout << "Portfolio return vector" << std::endl;
  std::cout << portfolioReturnVector << std::endl;

// Constant
  Rate c = .05;

// Portfolio return vector minus constant rate
  Matrix portfolioReturnVectorMinusC (4, 1);
  for (int i = 0; i < 4; ++i) {
    portfolioReturnVectorMinusC[i][0] = portfolioReturnVector[i][0] - c;
  }

  std::cout
    << boost::format ("Portfolio return vector minus constantrate (c = %f)") % c
    << std::endl;
  std::cout << portfolioReturnVectorMinusC << std::endl;

// Inverse of covariance matrix
  const Matrix &inverseOfCovarienceMatrix = inverse (covarianceMatrix);

// Z vectors
  const Matrix &portfolioAz = inverseOfCovarienceMatrix * portfolioReturnVector;
  std::cout << "Portfolio A z vector" << std::endl;
  std::cout << portfolioAz << std::endl;
  double sumOfPortfolioAz = 0.0;
  std::for_each (portfolioAz.begin (), portfolioAz.end (), [&] (Real n) {
      sumOfPortfolioAz += n;
  });

  const Matrix &portfolioBz =
    inverseOfCovarienceMatrix * portfolioReturnVectorMinusC;
  std::cout << "Portfolio B z vector" << std::endl;
  std::cout << portfolioBz << std::endl;
  double sumOfPortfolioBz = 0.0;
  std::for_each (portfolioBz.begin (), portfolioBz.end (), [&] (Real n) {
      sumOfPortfolioBz += n;
  });

// Portfolio weights
  Matrix weightsPortfolioA (4, 1);
  for (int i = 0; i < 4; ++i) {
    weightsPortfolioA[i][0] = portfolioAz[i][0] / sumOfPortfolioAz;
  }

  std::cout << "Portfolio A weights" << std::endl;
  std::cout << weightsPortfolioA << std::endl;

  Matrix weightsPortfolioB (4, 1);
  for (int i = 0; i < 4; ++i) {
    weightsPortfolioB[i][0] = portfolioBz[i][0] / sumOfPortfolioBz;
  }

  std::cout << "Portfolio B weights" << std::endl;
  std::cout << weightsPortfolioB << std::endl;

// Portfolio risk and return
  const Matrix &expectedReturnPortfolioAMatrix =
    transpose (weightsPortfolioA) * portfolioReturnVector;
  double expectedReturnPortfolioA = expectedReturnPortfolioAMatrix[0][0];
  const Matrix &variancePortfolioAMatrix =
    transpose (weightsPortfolioA) * covarianceMatrix * weightsPortfolioA;
  double variancePortfolioA = variancePortfolioAMatrix[0][0];
  double stdDeviationPortfolioA = std::sqrt (variancePortfolioA);
  std::cout << boost::format ("Portfolio A expected return: %f") %
               expectedReturnPortfolioA << std::endl;
  std::cout << boost::format ("Portfolio A variance: %f") % variancePortfolioA
            << std::endl;
  std::cout << boost::format ("Portfolio A standard deviation: %f") %
               stdDeviationPortfolioA << std::endl;

  const Matrix &expectedReturnPortfolioBMatrix =
    transpose (weightsPortfolioB) * portfolioReturnVector;
  double expectedReturnPortfolioB = expectedReturnPortfolioBMatrix[0][0];
  const Matrix &variancePortfolioBMatrix =
    transpose (weightsPortfolioB) * covarianceMatrix * weightsPortfolioB;
  double variancePortfolioB = variancePortfolioBMatrix[0][0];
  double stdDeviationPortfolioB = std::sqrt (variancePortfolioB);
  std::cout << boost::format ("Portfolio B expected return: %f") %
               expectedReturnPortfolioB << std::endl;
  std::cout << boost::format ("Portfolio B variance: %f") % variancePortfolioB
            << std::endl;
  std::cout << boost::format ("Portfolio B standard deviation: %f") %
               stdDeviationPortfolioB << std::endl;

// Covariance and correlation of returns
  const Matrix &covarianceABMatrix =
    transpose (weightsPortfolioA) * covarianceMatrix * weightsPortfolioB;
  double covarianceAB = covarianceABMatrix[0][0];
  double correlationAB =
    covarianceAB / (stdDeviationPortfolioA * stdDeviationPortfolioB);
  std::cout
    << boost::format ("Covariance of portfolio A and B: %f") % covarianceAB
    << std::endl;
  std::cout
    << boost::format ("Correlation of portfolio A and B: %f") % correlationAB
    << std::endl;

// Generate envelope set of portfolios
  double startingProportion = -.40;
  double increment = .10;
  std::map<double, std::pair<Volatility, double> > mapOfProportionToRiskAndReturn;
  std::map<Volatility, double> mapOfVolatilityToReturn;
  for (int i = 0; i < 21; ++i) {
    double proportionA = startingProportion + i * increment;
    Volatility risk_frontier = calculatePortfolioRisk (proportionA,
                                                       stdDeviationPortfolioA,
                                                       stdDeviationPortfolioB,
                                                       covarianceAB);
    double returnEF = calculatePortfolioReturn (proportionA,
                                                expectedReturnPortfolioA,
                                                expectedReturnPortfolioB);
    mapOfProportionToRiskAndReturn[proportionA] = std::make_pair (risk_frontier,
                                                                  returnEF);
    mapOfVolatilityToReturn[risk_frontier] = returnEF;
  }

// Write data to a file for plotting latter
  std::ofstream envelopeSetFile;
  envelopeSetFile.open ("./envelope.dat", std::ios::out);
  for (std::map<double, std::pair<Volatility, double> >::const_iterator i = mapOfProportionToRiskAndReturn.begin ();
       i != mapOfProportionToRiskAndReturn.end (); ++i) {
    envelopeSetFile << boost::format ("%f %f %f") % i->first % i->second.first %
                       i->second.second << std::endl;
  }
  envelopeSetFile.close ();

// Find minimum risk portfolio on efficient frontier
  std::pair<Volatility, double> minimumVariancePortfolioRiskAndReturn = *mapOfVolatilityToReturn.begin ();
  Volatility minimumRisk = minimumVariancePortfolioRiskAndReturn.first;
  double maximumReturn = minimumVariancePortfolioRiskAndReturn.second;
  std::cout << boost::format ("Maximum portfolio return for risk of %f is %f") %
               minimumRisk % maximumReturn << std::endl;

// Generate efficient frontier
  std::map<Volatility, double> efficientFrontier;
  for (std::map<double, std::pair<Volatility, double> >::const_iterator i = mapOfProportionToRiskAndReturn.begin ();
       i != mapOfProportionToRiskAndReturn.end (); ++i) {
    efficientFrontier[i->second.first] = i->second.second;
    if (i->second.first == minimumRisk) break;
  }

// Write efficient frontier to file
  std::ofstream efficient_frontier_file;
  efficient_frontier_file.open ("./efficient_frontier.dat", std::ios::out);
  for (std::map<Volatility, double>::const_iterator i = efficientFrontier.begin ();
       i != efficientFrontier.end (); ++i) {
    efficient_frontier_file << boost::format ("%f %f") % i->first % i->second
                            << std::endl;
  }
  efficient_frontier_file.close ();
  return 0;
}

