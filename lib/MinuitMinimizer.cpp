/*
 * MinuitMinimizer.cpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2016 Antonin Portelli
 *
 * LatAnalyze 3 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatAnalyze 3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatAnalyze 3.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <LatAnalyze/MinuitMinimizer.hpp>
#include <LatAnalyze/includes.hpp>
#include <Minuit2/Minuit2Minimizer.h>
#include <Math/Functor.h>

using namespace std;
using namespace Latan;

static constexpr double       initErr = 0.1;
static constexpr unsigned int maxTry  = 10u;

/******************************************************************************
 *                 MinuitMinimizer implementation                          *
 ******************************************************************************/
// constructors ////////////////////////////////////////////////////////////////
MinuitMinimizer::MinuitMinimizer(const Algorithm algorithm)
{
    setAlgorithm(algorithm);
}

MinuitMinimizer::MinuitMinimizer(const Index dim, const Algorithm algorithm)
: Minimizer(dim)
{
    setAlgorithm(algorithm);
}

// access //////////////////////////////////////////////////////////////////////
MinuitMinimizer::Algorithm MinuitMinimizer::getAlgorithm(void) const
{
    return algorithm_;
}

void MinuitMinimizer::setAlgorithm(const Algorithm algorithm)
{
    algorithm_ = algorithm;
}

// minimization ////////////////////////////////////////////////////////////////
const DVec & MinuitMinimizer::operator()(const DoubleFunction &f)
{
    using namespace ROOT;
    using namespace Minuit2;
    
    DVec                        &x = getState();
    int                         printLevel;
    EMinimizerType              minuitAlg;
    unique_ptr<Math::Minimizer> min;
    
    // convert Latan parameters to Minuit parameters
    switch (getVerbosity())
    {
        case Verbosity::Silent:
            printLevel = 0;
            break;
        case Verbosity::Normal:
            printLevel = 2;
            break;
        case Verbosity::Debug:
            printLevel = 3;
            break;
    }
    switch (getAlgorithm())
    {
        case Algorithm::Migrad:
            minuitAlg = kMigrad;
            break;
        case Algorithm::Simplex:
            minuitAlg = kSimplex;
            break;
        case Algorithm::Combined:
            minuitAlg = kCombined;
            break;
    }
    
    // resize minimizer state to match function number of arguments
    if (f.getNArg() != x.size())
    {
        resize(f.getNArg());
    }
    
    // create and set minimizer
    min.reset(new Minuit2Minimizer(minuitAlg));
    min->SetMaxFunctionCalls(getMaxIteration());
    min->SetTolerance(getPrecision());
    min->SetPrintLevel(printLevel);
    
    // set function and variables
    Math::Functor minuitF(f, x.size());
    string        name;
    double        val, step;
    
    min->SetFunction(minuitF);
    for (Index i = 0; i < x.size(); ++i)
    {
        name = "x_" + strFrom(i);
        val  = x(i);
        step = (fabs(x(i)) != 0.) ? initErr*fabs(x(i)) : 1.;
        if (hasHighLimit(i) and !hasLowLimit(i))
        {
            min->SetUpperLimitedVariable(i, name, val, step, getHighLimit(i));
        }
        else if (!hasHighLimit(i) and hasLowLimit(i))
        {
            min->SetLowerLimitedVariable(i, name, val, step, getLowLimit(i));
        }
        else if (hasHighLimit(i) and hasLowLimit(i))
        {
            min->SetLimitedVariable(i, name, val, step, getLowLimit(i),
                                    getHighLimit(i));
        }
        else
        {
            min->SetVariable(i, name, val, step);
        }
    }
    
    // minimize
    int          status;
    unsigned int n = 0;
    
    if (getVerbosity() >= Verbosity::Normal)
    {
        cout << "========== Minuit pre-minimization " << endl;
    }
    min->SetStrategy(0);
    min->Minimize();
    do
    {
        n++;
        if (getVerbosity() >= Verbosity::Normal)
        {
            cout << "========== Minuit minimization, try #" << n << endl;
        }
        min->SetStrategy(2);
        min->Minimize();
        status = min->Status();
    } while (status and (n < maxTry));
    if (getVerbosity() >= Verbosity::Normal)
    {
        cout << "==============================" << endl;
    }
    switch (status)
    {
        case 1:
            LATAN_WARNING("invalid minimum: covariance matrix was made positive");
            break;
        case 2:
            LATAN_WARNING("invalid minimum: Hesse analysis is not valid");
            break;
        case 3:
            LATAN_WARNING("invalid minimum: requested precision not reached");
            break;
        case 4:
            LATAN_WARNING("invalid minimum: iteration limit reached");
            break;
    }
    
    // save and return result
    for (Index i = 0; i < x.size(); ++i)
    {
        x(i) = min->X()[i];
    }
    
    return x;
}
