#include <iostream>
#include <cmath>
#include <LatAnalyze/CompiledModel.hpp>
#include <LatAnalyze/MinuitMinimizer.hpp>
#include <LatAnalyze/RandGen.hpp>
#include <LatAnalyze/XYSampleData.hpp>

using namespace std;
using namespace Latan;

const Index  nPoint1 = 10, nPoint2 = 10;
const Index  nSample = 1000;
const double xErr = .1, yErr   = .1;
const double exactPar[2] = {0.5,5.};
const double dx1 = 10.0/static_cast<double>(nPoint1);
const double dx2 = 5.0/static_cast<double>(nPoint2);

int main(void)
{
    // generate fake data
    XYSampleData data(nSample);
    RandGen      rg;
    double       xBuf[2];
    DoubleModel  f([](const double *x, const double *p)
                   {return p[1]*exp(-x[0]*p[0])+x[1];}, 2, 2);
    
    data.addXDim("x", nPoint1);
    data.addXDim("off", nPoint2);
    data.addYDim("y");
    for (Index s = central; s < nSample; ++s)
    {
        for (Index i1 = 0; i1 < nPoint1; ++i1)
        {
            xBuf[0]          = i1*dx1;
            data.x(i1, 0)[s] = rg.gaussian(xBuf[0], xErr);
            for (Index i2 = 0; i2 < nPoint2; ++i2)
            {
                xBuf[1]                           = i2*dx2;
                data.x(i2, 1)[s]                  = xBuf[1];
                data.y(data.dataIndex(i1, i2))[s] =
                    rg.gaussian(f(xBuf, exactPar), yErr);
            }
        }
    }
    data.assumeXExact(true, 1);
    cout << data << endl;
    
    // fit
    DVec init = DVec::Constant(2, 0.1);
    DMat err;
    SampleFitResult p;
    MinuitMinimizer minimizer;

    p   = data.fit(minimizer, init, f);
    err = p.variance().cwiseSqrt();
    cout << "a= " << p[central](0) << " +/- " << err(0) << endl;
    cout << "b= " << p[central](1) << " +/- " << err(1) << endl;
    cout << "chi^2/ndof= " << p.getChi2PerDof() << endl;
    cout << "p-value= " << p.getPValue() << endl;
    
    return EXIT_SUCCESS;
}
