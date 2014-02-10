/*
 * Io.hpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2014 Antonin Portelli
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

#ifndef Latan_Io_hpp_
#define	Latan_Io_hpp_

#include <fstream>
#include <map>
#include <stack>
#include <sstream>
#include <string>
#include <latan/Global.hpp>
#include <latan/IoObject.hpp>
#include <latan/Mat.hpp>
#include <latan/ParserState.hpp>
#include <latan/RandGen.hpp>
#include <latan/Sample.hpp>

BEGIN_NAMESPACE

/******************************************************************************
 *                          Generic datafile class                            *
 ******************************************************************************/
typedef std::map<std::string, IoObject *> IoDataTable;

class File
{
public:
    class Mode
    {
    public:
        enum
        {
            null   = 0,
            write  = 1 << 0,
            read   = 1 << 1,
            append = 1 << 2
        };
    };
public:
    // constructors
    File(void);
    File(const std::string &name, const unsigned int mode);
    // destructor
    virtual ~File(void);
    // access
    std::string  getName(void) const;
    unsigned int getMode(void) const;
    template <typename IoT>
    const IoT &  read(const std::string &name);
    virtual void save(const DMat &m, const std::string &name)               = 0;
    virtual void save(const RandGen::State &state, const std::string &name) = 0;
    // tests
    virtual bool isOpen(void) const = 0;
    // IO
    virtual void close(void)                                            = 0;
    virtual void open(const std::string &name, const unsigned int mode) = 0;
protected:
    // data access
    void         deleteData(void);
    template <typename IoT>
    const IoT& getData(const std::string &name) const;
    // IO
    virtual void load(const std::string &name) = 0;
    // error checking
    void checkWritability(void);
protected:
    std::string  name_;
    unsigned int mode_;
    IoDataTable  data_;
};

// Template implementations
template <typename IoT>
const IoT& File::read(const std::string &name)
{
    load(name);
    
    return getData<IoT>(name);
}

template <typename IoT>
const IoT& File::getData(const std::string &name) const
{
    IoDataTable::const_iterator i = data_.find(name);
    
    if (i != data_.end())
    {
        return dynamic_cast<const IoT&>(*(i->second));
    }
    else
    {
        LATAN_ERROR(Definition, "no data with name '" + name + "'");
    }
}

/******************************************************************************
 *                          ASCII datafile class                              *
 ******************************************************************************/
class AsciiFile: public File
{
public:
    class AsciiParserState: public ParserState<IoDataTable>
    {
    public:
        // constructor
        explicit AsciiParserState(std::istream *stream, std::string *name,\
                                  IoDataTable *data);
        // destructor
        virtual ~AsciiParserState(void);
        // public members
        std::stack<DMat>   dMatBuf;
        std::stack<double> doubleBuf;
        std::stack<int>    intBuf;
    private:
        // allocation/deallocation functions defined in IoAsciiLexer.lpp
        virtual void initScanner(void);
        virtual void destroyScanner(void);
    };
public:
    // constructors
    AsciiFile(void);
    AsciiFile(const std::string &name, const unsigned int mode);
    // destructor
    virtual ~AsciiFile(void);
    // access
    virtual void save(const DMat &m, const std::string &name);
    virtual void save(const RandGen::State &state, const std::string &name);
    // tests
    virtual bool isOpen(void) const;
    // IO
    virtual void close(void);
    virtual void open(const std::string &name, const unsigned int mode);
private:
    // IO
    virtual void load(const std::string &name);
    // parser
    void parse(void);
private:
    std::fstream      fileStream_;
    bool              isParsed_;
    AsciiParserState* state_;
};

END_NAMESPACE

#endif
