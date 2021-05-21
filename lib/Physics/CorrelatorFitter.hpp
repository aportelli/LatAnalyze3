/*
 * CorrelatorFitter.hpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2020 Antonin Portelli
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

#ifndef Latan_CorrelatorFitter_hpp_
#define Latan_CorrelatorFitter_hpp_

#include <LatAnalyze/Global.hpp>
#include <LatAnalyze/Functional/Model.hpp>
#include <LatAnalyze/Numerical/FFT.hpp>
#include <LatAnalyze/Statistics/XYSampleData.hpp>

BEGIN_LATAN_NAMESPACE

/******************************************************************************
 *                         Correlator types & models                          *
 ******************************************************************************/
enum class CorrelatorType {undefined, exp, cosh, sinh, linear, cst};

namespace CorrelatorModels
{
    struct ModelPar
    {
        CorrelatorType type;
        Index          nState;
    };

    DoubleModel makeExpModel(const Index nState);
    DoubleModel makeCoshModel(const Index nState, const Index nt);
    DoubleModel makeSinhModel(const Index nState, const Index nt);
    DoubleModel makeConstModel(void);
    DoubleModel makeLinearModel(void);
    ModelPar    parseModel(const std::string s);
    DoubleModel makeModel(const ModelPar par, const Index nt);
    DVec        parameterGuess(const DMatSample &corr, const ModelPar par);
};

/******************************************************************************
 *                         Correlator utilities                               *
 ******************************************************************************/
namespace CorrelatorUtils
{
    DMatSample shift(const DMatSample &c, const Index ts);
    DMatSample fold(const DMatSample &c, const CorrelatorModels::ModelPar &par);
    DMatSample fourierTransform(const DMatSample &c, FFT &fft, 
                                const unsigned int dir = FFT::Forward);
};

/******************************************************************************
 *                       Correlator fit utility class                         *
 ******************************************************************************/
class CorrelatorFitter
{
public:
    // constructors
    CorrelatorFitter(const DMatSample &corr);
    CorrelatorFitter(const std::vector<DMatSample> &corr);
    // destructor
    virtual ~CorrelatorFitter(void) = default;
    // access
    XYSampleData & data(void);
    void setCorrelator(const DMatSample &corr);
    void setCorrelators(const std::vector<DMatSample> &corr);
    const DMatSample & getCorrelator(const Index i = 0) const;
    const std::vector<DMatSample> & getCorrelators(void) const;
    void setModel(const DoubleModel &model, const Index i = 0);
    const DoubleModel & getModel(const Index i = 0) const;
    void setFitRange(const Index tMin, const Index tMax, const Index i = 0);
    void setCorrelation(const bool isCorrelated, const Index i = 0, 
                        const Index j = 0);
    DMat getVarianceMatrix(void) const;
    void setThinning(const Index thinning, const Index i = 0);
    // fit functions
    SampleFitResult fit(Minimizer &minimizer, const DVec &init);
    SampleFitResult fit(std::vector<Minimizer *> &minimizer, const DVec &init);
private:
    // internal function to refresh fit ranges
    void refreshRanges(void);
private:
    Index                                nt_{0};
    std::unique_ptr<XYSampleData>        data_;
    std::vector<DoubleModel>             model_;
    std::vector<std::pair<Index, Index>> range_;
    std::vector<Index>                   thinning_;
};

END_LATAN_NAMESPACE

#endif // Latan_CorrelatorFitter_hpp_