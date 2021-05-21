/*
 * XYSampleData.cpp, part of LatAnalyze 3
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

#include <LatAnalyze/Statistics/XYSampleData.hpp>
#include <LatAnalyze/includes.hpp>
#include <LatAnalyze/Core/Math.hpp>

using namespace std;
using namespace Latan;

/******************************************************************************
 *                      SampleFitResult implementation                        *
 ******************************************************************************/
double SampleFitResult::getChi2(const Index s) const
{
    return chi2_[s];
}

const DSample & SampleFitResult::getChi2(const PlaceHolder ph __dumb) const
{
    return chi2_;
}

double SampleFitResult::getChi2PerDof(const Index s) const
{
    return chi2_[s]/getNDof();
}

DSample SampleFitResult::getChi2PerDof(const PlaceHolder ph __dumb) const
{
    return chi2_/getNDof();
}

double SampleFitResult::getNDof(void) const
{
    return static_cast<double>(nDof_);
}

Index SampleFitResult::getNPar(void) const
{
    return nPar_;
}

double SampleFitResult::getPValue(const Index s) const
{
    return Math::chi2PValue(getChi2(s), getNDof());
}

double SampleFitResult::getCcdf(const Index s) const
{
    return Math::chi2Ccdf(getChi2(s), getNDof());
}

const DoubleFunction & SampleFitResult::getModel(const Index s,
                                                 const Index j) const
{
    return model_[j][s];
}

const DoubleFunctionSample & SampleFitResult::getModel(
                                                    const PlaceHolder ph __dumb,
                                                    const Index j) const
{
    return model_[j];
}

FitResult SampleFitResult::getFitResult(const Index s) const
{
    FitResult fit;
    
    fit = (*this)[s];
    fit.chi2_ = getChi2();
    fit.nDof_ = static_cast<Index>(getNDof());
    fit.model_.resize(model_.size());
    for (unsigned int k = 0; k < model_.size(); ++k)
    {
        fit.model_[k] = model_[k][s];
    }
    
    return fit;
}

// IO //////////////////////////////////////////////////////////////////////////
void SampleFitResult::print(const bool printXsi, ostream &out) const
{
    char  buf[256];
    Index pMax = printXsi ? size() : nPar_;
    DMat  err = this->variance().cwiseSqrt();
    
    sprintf(buf, "chi^2/dof= %.1e/%d= %.2e -- chi^2 CCDF= %.2e -- p-value= %.2e",
        getChi2(), static_cast<int>(getNDof()), getChi2PerDof(), getCcdf(), 
        getPValue());
    out << buf << endl;
    for (Index p = 0; p < pMax; ++p)
    {
        sprintf(buf, "%8s= % e +/- %e", parName_[p].c_str(),
                (*this)[central](p), err(p));
        out << buf << endl;
    }
}

/******************************************************************************
 *                       XYSampleData implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
XYSampleData::XYSampleData(const Index nSample)
: nSample_(nSample)
{}

// data access /////////////////////////////////////////////////////////////////
DSample & XYSampleData::x(const Index r, const Index i)
{
    checkXIndex(r, i);
    scheduleDataInit();
    scheduleComputeVarMat();
    if (xData_[i][r].size() == 0)
    {
        xData_[i][r].resize(nSample_);
    }
    
    return xData_[i][r];
}

const DSample & XYSampleData::x(const Index r, const Index i) const
{
    checkXIndex(r, i);
    
    return xData_[i][r];
}

const DMatSample & XYSampleData::x(const Index k)
{
    checkDataIndex(k);
    
    updateXMap();
    
    return xMap_.at(k);
}

DSample & XYSampleData::y(const Index k, const Index j)
{
    checkYDim(j);
    if (!pointExists(k, j))
    {
        registerDataPoint(k, j);
    }
    scheduleDataInit();
    scheduleComputeVarMat();
    if (yData_[j][k].size() == 0)
    {
        yData_[j][k].resize(nSample_);
    }
    
    return yData_[j][k];
}

const DSample & XYSampleData::y(const Index k, const Index j) const
{
    checkPoint(k, j);
    
    return yData_[j].at(k);
}

void XYSampleData::setUnidimData(const DMatSample &xData,
                                 const vector<const DMatSample *> &v)
{
    FOR_STAT_ARRAY(xData, s)
    FOR_VEC(xData[central], r)
    {
        x(r, 0)[s] = xData[s](r);
        for (unsigned int j = 0; j < v.size(); ++j)
        {
            y(r, j)[s] = (*(v[j]))[s](r);
        }
    }
}

const DMat & XYSampleData::getXXVar(const Index i1, const Index i2)
{
    checkXDim(i1);
    checkXDim(i2);
    computeVarMat();
    
    return data_.getXXVar(i1, i2);
}

const DMat & XYSampleData::getYYVar(const Index j1, const Index j2)
{
    checkYDim(j1);
    checkYDim(j2);
    computeVarMat();
    
    return data_.getYYVar(j1, j2);
}

const DMat & XYSampleData::getXYVar(const Index i, const Index j)
{
    checkXDim(i);
    checkYDim(j);
    computeVarMat();
    
    return data_.getXYVar(i, j);
}

DVec XYSampleData::getXError(const Index i)
{
    checkXDim(i);
    computeVarMat();
    
    return data_.getXError(i);
}

DVec XYSampleData::getYError(const Index j)
{
    checkYDim(j);
    computeVarMat();
    
    return data_.getYError(j);
}

bool XYSampleData::checkFit()
{
    return goodFit_;
}

void XYSampleData::checkChi2PerDof(double Chi2PerDof)
{
    if(Chi2PerDof >= 2 or Chi2PerDof < 0 or isnan(Chi2PerDof)) 
    {
        goodFit_ = false;
    }
}

// get total fit variance matrix and its pseudo-inverse ////////////////////////
const DMat & XYSampleData::getFitVarMat(void)
{
    computeVarMat();
    
    return data_.getFitVarMat();
}

const DMat & XYSampleData::getFitVarMatPInv(void)
{
    computeVarMat();
    
    return data_.getFitVarMatPInv();
}

// set data to a particular sample /////////////////////////////////////////////
void XYSampleData::setDataToSample(const Index s)
{
    if (initData_ or (s != dataSample_))
    {
        for (Index i = 0; i < getNXDim(); ++i)
        for (Index r = 0; r < getXSize(i); ++r)
        {
            data_.x(r, i) = xData_[i][r][s];
        }
        for (Index j = 0; j < getNYDim(); ++j)
        for (auto &p: yData_[j])
        {
            data_.y(p.first, j) = p.second[s];
        }
        dataSample_ = s;
        initData_   = false;
    }
}

// get internal XYStatData /////////////////////////////////////////////////////
const XYStatData & XYSampleData::getData(void)
{
    setDataToSample(central);
    computeVarMat();
    
    return data_;
}

// fit /////////////////////////////////////////////////////////////////////////
SampleFitResult XYSampleData::fit(std::vector<Minimizer *> &minimizer,
                                  const DVec &init,
                                  const std::vector<const DoubleModel *> &v)
{
    computeVarMat();
    
    SampleFitResult result;
    FitResult       sampleResult;
    DVec            initCopy = init;
    
    result.resize(nSample_);
    result.chi2_.resize(nSample_);
    result.model_.resize(v.size());
    double chi2PerDof;
    goodFit_ = true;
    FOR_STAT_ARRAY(result, s)
    {
        if(goodFit_)
        {
            setDataToSample(s);
            if (s == central)
            {
                sampleResult = data_.fit(minimizer, initCopy, v);
                initCopy     = sampleResult.segment(0, initCopy.size());
                chi2PerDof   = sampleResult.getChi2PerDof();
                checkChi2PerDof(chi2PerDof);
            }
            else
            {
                sampleResult = data_.fit(*(minimizer.back()), initCopy, v);
                chi2PerDof   = sampleResult.getChi2PerDof();
            }
            result[s]       = sampleResult;
            result.chi2_[s] = sampleResult.getChi2();
            for (unsigned int j = 0; j < v.size(); ++j)
            {
                result.model_[j].resize(nSample_);
                result.model_[j][s] = sampleResult.getModel(j);
            }
        }
    }
    result.nPar_    = sampleResult.getNPar();
    result.nDof_    = sampleResult.nDof_;
    result.parName_ = sampleResult.parName_;
    
    return result;
}

SampleFitResult XYSampleData::fit(Minimizer &minimizer,
                                  const DVec &init,
                                  const std::vector<const DoubleModel *> &v)
{
    vector<Minimizer *> mv{&minimizer};
    
    return fit(mv, init, v);
}

// residuals ///////////////////////////////////////////////////////////////////
XYSampleData XYSampleData::getResiduals(const SampleFitResult &fit)
{
    XYSampleData res(*this);
    
    for (Index j = 0; j < getNYDim(); ++j)
    {
        const DoubleFunctionSample &f = fit.getModel(_, j);
        
        for (auto &p: yData_[j])
        {
            res.y(p.first, j) -= f(x(p.first));
        }
    }
    
    return res;
}

XYSampleData XYSampleData::getPartialResiduals(const SampleFitResult &fit,
                                               const DVec &ref, const Index i)
{
    XYSampleData res(*this);
    DMatSample   buf(nSample_);
    
    buf.fill(ref);
    for (Index j = 0; j < getNYDim(); ++j)
    {
        const DoubleFunctionSample &f = fit.getModel(_, j);
        
        for (auto &p: yData_[j])
        {
            FOR_STAT_ARRAY(buf, s)
            {
                buf[s](i) = x(p.first)[s](i);
            }
            res.y(p.first, j) -= f(x(p.first)) - f(buf);
        }
    }
    
    return res;
}

// buffer list of x vectors ////////////////////////////////////////////////////
void XYSampleData::scheduleXMapInit(void)
{
    initXMap_ = true;
}

void XYSampleData::updateXMap(void)
{
    if (initXMap_)
    {
        for (Index s = central; s < nSample_; ++s)
        {
            setDataToSample(s);
            for (auto k: getDataIndexSet())
            {
                if (s == central)
                {
                    xMap_[k].resize(nSample_);
                }
                xMap_[k][s] = data_.x(k);
            }
        }
        initXMap_ = false;
    }
}

// schedule data initilization from samples ////////////////////////////////////
void XYSampleData::scheduleDataInit(void)
{
    initData_ = true;
}

// variance matrix computation /////////////////////////////////////////////////
void XYSampleData::scheduleComputeVarMat(void)
{
    computeVarMat_ = true;
}

void XYSampleData::computeVarMat(void)
{
    if (computeVarMat_)
    {
        // initialize data if necessary
        setDataToSample(central);

        // compute relevant sizes
        Index size = 0, ySize = 0;
        
        for (Index j = 0; j < getNYDim(); ++j)
        {
            size += getYSize(j);
        }
        ySize = size;
        for (Index i = 0; i < getNXDim(); ++i)
        {
            size += getXSize(i);
        }
        
        // compute total matrix
        DMatSample z(nSample_, size, 1);
        DMat       var;
        Index      a;
        
        FOR_STAT_ARRAY(z, s)
        {
            a = 0;
            for (Index j = 0; j < getNYDim(); ++j)
            for (auto &p: yData_[j])
            {
                z[s](a, 0) = p.second[s];
                a++;
            }
            for (Index i = 0; i < getNXDim(); ++i)
            for (Index r = 0; r < getXSize(i); ++r)
            {
                z[s](a, 0) = xData_[i][r][s];
                a++;
            }
        }
        var = z.varianceMatrix();
        
        // assign blocks to data
        Index a1, a2;
        
        a1 = ySize;
        for (Index i1 = 0; i1 < getNXDim(); ++i1)
        {
            a2 = ySize;
            for (Index i2 = 0; i2 < getNXDim(); ++i2)
            {
                data_.setXXVar(i1, i2,
                               var.block(a1, a2, getXSize(i1), getXSize(i2)));
                a2 += getXSize(i2);
            }
            a1 += getXSize(i1);
        }
        a1 = 0;
        for (Index j1 = 0; j1 < getNYDim(); ++j1)
        {
            a2 = 0;
            for (Index j2 = 0; j2 < getNYDim(); ++j2)
            {
                data_.setYYVar(j1, j2,
                               var.block(a1, a2, getYSize(j1), getYSize(j2)));
                a2 += getYSize(j2);
            }
            a1 += getYSize(j1);
        }
        a1 = ySize;
        for (Index i = 0; i < getNXDim(); ++i)
        {
            a2 = 0;
            for (Index j = 0; j < getNYDim(); ++j)
            {
                data_.setXYVar(i, j,
                               var.block(a1, a2, getXSize(i), getYSize(j)));
                a2 += getYSize(j);
            }
            a1 += getXSize(i);
        }
        computeVarMat_ = false;
    }
    if (initVarMat())
    {
        data_.copyInterface(*this);
        scheduleFitVarMatInit(false);
    }
}

// create data /////////////////////////////////////////////////////////////////
void XYSampleData::createXData(const string name, const Index nData)
{
    data_.addXDim(nData, name);
    xData_.push_back(vector<DSample>(nData));
}

void XYSampleData::createYData(const string name)
{
    data_.addYDim(name);
    yData_.push_back(map<Index, DSample>());
}
