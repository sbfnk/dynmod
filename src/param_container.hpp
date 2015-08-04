// -*- compile-command: "cd .. ; make -k; cd -"; -*-
/*! \file param_container.hpp
  \brief Header file for the ParamContainer class
*/

#ifndef PARAMCONTAINER_HPP
#define PARAMCONTAINER_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <fstream>

#include "parameter.hpp"

namespace po = boost::program_options;

/*! Base class Class for parameters */
class ParamContainer
{
public:

    //! Parameter information container
    struct ParamInfo {
        //! Constructor
        /*!
          \param[in] o Option name.
          \param[in] d Option descripton.
          \param[in] p Reference to option variable.
        */
        ParamInfo(std::string o, std::string d, Parameter* p, bool l = false) :
            option(o), description(d), param(p), logTransform(l) {;}

        //!< Command-line / config file option associated with the parameter
        std::string option;
        //!< Plain-text description of the parameter
        std::string description;
        //!< Pointer to the variable storing the parameter
        Parameter* param;
        //!< Do we log-transform the parameter?
        bool logTransform;
    };

    //! Data container
    /*!
      This holds time-varying data of type T, such as birth rates,
      population sizes, vaccination rates etc.
    */
    template <class T>
    class Data
    {
    public:
        //! Constructor
        Data(): lastTime(-1), currentValue(-1), counter(0) {}
        //! Accessor for data values at a given time.
        /*
          \param[in] time The time at which we want the data value.
          \param[in] prevTime If this is set to something >=0, getValue
          will interpolate the value of the data between prevTime and
          time, assuming it behaves linearly between the two.
        */
        T getValue(double time, double prevTime = -1.);
        //! Accessor for the last data value.
        T getLastValue() { return values.back(); }
        //! Set values at given times.
        void setValues(std::vector<double> const& newTimes,
                       std::vector<T> const& newValues);
        //! Add a value to the back of the data, with the given time.
        void addValue(double newTime, T newValue);
        //! Add a value to the back of the data, increasing the last time
        //! by the difference between the last two times.
        void addLastValue(T newValue);
    private:
        //! Time of the last query for data. This exists purely for
        //! performance reasons, so we don't have to perform any
        //! calculations if the time hasn't changed
        double lastTime;
        double currentValue; //!< Value of the last query
        size_t counter; //!< Counter for the time/value pair of the last query
        std::vector<double> times; //!< Data times
        std::vector<T> values; //!< Data values
    };

    //! Constructor.
    /*!
      This is where classes derived from Model define the parameters and command
      line options.

      \param[in] optionName Descriptive name of the options (e.g.,
      "Model parameters")
    */
    ParamContainer(std::string optionName = "") :
        options(po::options_description("\nOptions ("+optionName+")"))
    {;}
    //! Destructor
    virtual ~ParamContainer()
    {;}

    //! Read data table
    /*!
      Reads data from a CSV file, in simple (line,column) format,
      optionally reading in column headers, too.

      \param[in] fileName File that contains the table
      \param[out] data A container for the data in (header,column) format
    */
    template <class C, class R, class T>
    bool ReadTable(std::string fileName, std::vector<C> &columns,
                   std::vector<R> &rows, std::vector<std::vector<T> > &data,
                   bool transpose = false);

    //! Read command line parameters.
    /*!
      This should be called after the table has been read and command
      line parameters of the model haven been assigned. It initialises
      the model parameter variables with the values found in the command
      line parameters.

      \param[in] vm The map of command line parameters
    */
    bool ReadParams(po::variables_map &vm);

    //! Read parameters from input values.
    /*!
      Reads in the specified parameters value-by-value

      \param[in] readParams Vector of parameters that are to be read in.
      \param[in] values Vector of parameter values
    */
    int ReadInputValues(std::vector<std::string> readParams,
                        std::vector<double> values);


    //! Set all parameters to random values
    /*!
      Initialises all the parameters to random values

      \param[in] seed Random number generator seed
    */
    void setRandomParams(int& seed);

    //! Accessor for options
    const po::options_description getOptions() const
    { return options; }

    //! Accessor for name
    const std::string getName() const
    { return name; }

    //! Set a parameter value
    bool setParam(std::string name, double value);

    //! Print all parameter values
    void Print() const;

protected:

    /*! \brief Parameters

      A map of command line options to the species parameters
    */
    std::vector<ParamInfo> params;
    //! Command line options
    po::options_description options;

    std::string name;

};

template <class C, class R, class T>
bool ParamContainer::ReadTable
(std::string fileName, std::vector<C> &columns,
 std::vector<R> &rows, std::vector<std::vector<T> > &data,
 bool transpose)
{
    // tokenizer for reading csv file
    typedef boost::tokenizer<boost::escaped_list_separator<char> > Tokenizer;
    boost::escaped_list_separator<char> sep('\\', ',', '\"');

    std::ifstream in(fileName.c_str());
    if (!in.is_open())
    {
        std::cerr << "ERROR: Could not open " << fileName << std::endl;
        return false;
    }

    std::string line;
    bool firstLine = true;

    // read csv file line-by-line
    while (std::getline(in,line))
    {
        std::vector<std::string> lineVector;
        Tokenizer tok(line, sep);
        lineVector.assign(tok.begin(), tok.end());

        // store colum names if requested and we're in the first line
        if (firstLine)
        {
            firstLine = false;
            for (std::vector<std::string>::iterator it = (lineVector.begin() + 1);
                 it != lineVector.end(); it++)
            {
                C column;
                std::stringstream ss(*it);
                ss >> column;
                columns.push_back(column);
                if (transpose)
                    data.push_back(std::vector<T>());
            }
        }
        else
        {
            R row;
            std::stringstream ss(lineVector[0]);
            ss >> row;
            rows.push_back(row);

            if (!transpose)
                data.push_back(std::vector<T>());

            for (size_t i = 1; i < lineVector.size(); ++i)
            {
                T t;
                std::stringstream ss(lineVector[i]);
                ss >> t;
                if (transpose)
                    data[i-1].push_back(t);
                else
                    data.back().push_back(t);
            }
        }
    }

    in.close();
    return true;
}

template <class T>
void ParamContainer::Data<T>::setValues(std::vector<double> const& newTimes,
                                        std::vector<T> const& newValues)
{
    times = newTimes;
    values = newValues;

    if (times.size() != values.size())
        values.resize(times.size(), 0);
}

template <class T>
void ParamContainer::Data<T>::addValue(double newTime, T newValue)
{
    times.push_back(newTime);
    values.push_back(newValue);
}

template <class T>
void ParamContainer::Data<T>::addLastValue(T newValue)
{
    // we have at least two times, so the time for the new value will be
    // equally spaced to the last time as the two times before
    if (times.size() > 1)
    {
        times.push_back(times.back() + times.back() - times[times.size()-2]);
        values.push_back(newValue);
    }
    // we have exactly one time -- we increase it by one
    else if (times.size() > 0)
    {
        times.push_back(times.back() + 1);
        values.push_back(newValue);
    }
    // we have no time yet -- set it to zero
    else
    {
        times.push_back(.0);
        values.push_back(newValue);
    }
}

template <class T>
T ParamContainer::Data<T>::getValue(double time, double prevTime)
{
    if (time == lastTime) return currentValue; // shortcut

    if (time < lastTime) counter = 0; // reset counter if we're asking
    // for an earlier time

    // get counter that is closest to the requested time
    while ((counter + 1) < (times.size() - (prevTime >= 0 ? 1 : 0)) &&
           ((time - times[counter + 1]) > times[counter + 1] - lastTime))
        ++counter;

    // do we want an interpolated value?
    if ((prevTime >= 0) & (time > lastTime))
    {
        currentValue = (time - prevTime) /
            (times[counter + 1] - times[counter]) * values[counter];
    }
    else
        currentValue = values[counter];

    lastTime = time;
    return currentValue;
}

#endif
