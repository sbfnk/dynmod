/*! \file tsir.cpp
  \brief Main program for the TSIR model
*/

#include <boost/thread/thread.hpp>

#include "getseed.hpp"
#include "simulation.hpp"

/*! Main program. */
int main(int argc, char* argv[])
{

  // define tokenizer, to be used later for splitting a string of
  // comma-separated strings
  typedef boost::tokenizer<boost::escaped_list_separator<char> >
    Tokenizer;
  boost::escaped_list_separator<char> sep('\\', ',', '\"');

  //! This defines the main variables.

  std::string outFile; //!< output file

  unsigned int verbose = 0; //!< be verbose

  unsigned int nsim = 1;         //!< number of simulations to perform
  unsigned int simNb = 0;       //!< simulation number (if only one)
  int seed = 0; //!< simulation random seed

  //! Number of threads for parallelisation. If not given, it will use
  //! the maximum available.
  unsigned int nThreads = 0;

  //! File containing correlations between non-zero case numbers in
  //! cities, to compare the data to to estimate theta
  std::string correlationsFile;

  //! base directory
  std::string baseDir =
    std::string(getenv("HOME")) + "/code/measles/data/";

  std::vector<size_t> selectCities; //!< cities to be recorded

  //! Whether we're wrapped into an R function, i.e. expecting an
  //! input file and calculating summary statistics etc.
  bool wrapped = false;

  //! Then it reads in the command line parameters.

  //! main options
  po::options_description generic
    ("Usage: tsir [options]... \n\nOptions");

  generic.add_options()
    ("help,h",
     "produce help message")
    ("verbose,v",
     "produce verbose output")
    ("very-verbose,V",
     "produce very verbose output")
    ("config,f",
     po::value<std::string>()->default_value
     (baseDir + "config"),
     "config file")
    ;

  //! config options
  po::options_description config
    ("Configuration");

  config.add_options()
    ("output-file,o", po::value<std::string>(),
     "output file (\"-\" for output to stdout")
    ("noheader,d",
     "do not print header")
    ("print,p",
     "print results")
    ("nsim,n", po::value<unsigned int>()->default_value(1),
     "number of simulations")
    ("nthread,t", po::value<unsigned int>()->default_value(0),
     "number of threads")
    ("simnb", po::value<unsigned int>()->default_value(0),
     "run number (if only one simulation)")
    ("record,r", po::value<double>()->default_value(1.),
     "timestep of recording data")
    ("timestep", po::value<std::string>()->default_value("1"),
     "timestep")
    ("tstart,s", po::value<double>()->default_value(0.),
     "starting time")
    ("tend,e", po::value<double>()->default_value(100.),
     "end time")
    ("select,c", po::value<std::string>()->default_value("all"),
     "select a subset of cities")
    ("correlations", po::value<std::string>(),
     "get correlations from file (and calculate loss function)")
    ("introductions", po::value<double>()->default_value(.0),
     "frequency of introductions")
    ("future-vaccination", po::value<double>()->default_value(-1.),
     "future vaccination coverage (-1 for last known coverage)")
    ("vary-vaccination", po::value<int>()->default_value(0),
     "vary future vaccination coverage")
    ("input", po::value<std::string>(),
     "variables to read from input file(s)")
    ("number", po::value<unsigned int>()->default_value(0),
     "number of the input file (e.g., 1 for input1, 2 for input2 etc.)")
    ("summary", po::value<std::string>(),
     "summary statistics to calculate")
    ("seed", po::value<int>()->implicit_value(0),
     "random seed to use (0 for system random number generator)")
    ;

  po::positional_options_description input_number;
  input_number.add("number", -1);

  //! command line options
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  //! config file options
  po::options_description config_file_options;
  config_file_options.add(config);

  //! stores the options
  po::variables_map vm;

  try
    {
      po::store(po::command_line_parser(argc, argv).options(cmdline_options).
                allow_unregistered().run(), vm);
    }
  catch (std::exception& e)
    {
      std::cerr << "Error parsing command line parameters: " << e.what()
                << std::endl;
      return 1;
    }
  po::notify(vm);

  // help option
  if (vm.count("help"))
    {
      std::cout << generic << std::endl;
      return 0;
    }

  // assign number of threads for parallelisation
  nThreads = vm["nthread"].as<unsigned int>();
  if (nThreads == 0)
    nThreads = boost::thread::hardware_concurrency();

  // verbose options
  if (vm.count("verbose"))
    verbose = 1;
  if (vm.count("very-verbose"))
    {
      if (vm["nthread"].as<unsigned int>() == 0)
        nThreads = 1;
      verbose = 2;
    }

  if (verbose)
    {
      std::cout << "Parallelising to " << nThreads << " thread(s)."
                << std::endl;
    }

  // Read config file
  if (vm.count("config"))
    {
      std::ifstream in(vm["config"].as<std::string>().c_str());
      if (!in.is_open())
        {
          std::cerr << "WARNING: Could not open config file "
                    << vm["config"].as<std::string>() << std::endl;
        }
      else
        {
          try
            {
              po::store(po::parse_config_file(in, config_file_options), vm);
              po::notify(vm);
            }
          catch ( const boost::program_options::error& e )
            {
              std::cerr << e.what() << std::endl;
            }
          in.close();
        }
    }

  // Set output file
  std::streambuf *buf = 0;
  std::ofstream of;

  if (vm.count("output-file"))
    {
      std::string outFile = vm["output-file"].as<std::string>();
      if (outFile == "-")
        buf = std::cout.rdbuf();
      else
        {
          of.open(outFile.c_str());
          buf = of.rdbuf();
        }
    }

  std::ostream out(buf);

  // Number of simulations
  nsim = vm["nsim"].as<unsigned int>();
  // Simulation number (if only one)
  simNb = vm["simnb"].as<unsigned int>();
  if (nsim > 1 && simNb > 0)
    {
      std::cerr << "WARNING: simnb given with more than one simulation. "
                << "Ignoring." << std::endl;
    }

  if (verbose)
    std::cerr << "Running " << nsim << " simulations." << std::endl;

  // Read command line model options
  GravityModel model(baseDir, verbose);
  cmdline_options.add(model.getOptions());
  try
    {
      po::store
        (po::command_line_parser(argc, argv).options(cmdline_options).
         positional(input_number).run(), vm);
    }
  catch (std::exception& e)
    {
      std::cerr << "Error parsing command line parameters: " << e.what()
                << std::endl;
      return 1;
    }
  po::notify(vm);

  // Read model parameters, including data files
  if (!model.InitParams(vm))
    return 1;
  if (vm.count("input")) {

    wrapped = true;

    // Tokenizer for reading a comma-separated list of variables to
    // read from input file
    Tokenizer tok(vm["input"].as<std::string>(), sep);
    std::vector<std::string> lineVector;
    lineVector.assign(tok.begin(), tok.end());

    std::stringstream fileName;

    // add inputNb to filename if >0
    if (vm["number"].as<unsigned int>() > 0)
      fileName << "input" << vm["number"].as<unsigned int>();
    else
      fileName << "input";

    // read parameters from file, line-by-line
    std::ifstream in(fileName.str().c_str());

    std::string line;
    std::vector<double> values;

    while (std::getline(in, line))
      {
        std::istringstream iss(line);
        double lineValue;
        iss >> lineValue;
        values.push_back(lineValue);
      }

    if (values.size() < lineVector.size())
      {
        std::cerr << "ERROR: not enough lines in input file "
                  << fileName.str() << " to obtain all the parameters "
                  << "specified" << std::endl;
        return -1;
      }

    seed = model.ReadInputValues(lineVector, values);
  }

  if (vm.count("seed") && vm["seed"].as<int>() != 0)
    seed = vm["seed"].as<int>();

  if (seed == 0)
    seed = getSeed();

  //! Main simulation variable.
  Simulation sim(model, nThreads, verbose);

  // assign simulation parameters from comand line parameters

  if (vm["record"].as<double>() >=0 )
    sim.setRecordStep(vm["record"].as<double>());
  else
    {
      std::cerr << "ERROR: record has to be >=0" << std::endl;
      return -1;
    }

  sim.setStartTime(vm["tstart"].as<double>());
  sim.setEndTime(vm["tend"].as<double>());

  if (vm.count("timestep"))
    sim.setTimeStep(vm["timestep"].as<std::string>());

  if (vm.count("introductions"))
    sim.setIntroductionFrequency(vm["introductions"].as<double>());

  if (vm.count("future-vaccination"))
    sim.setFutureVaccinationRate(vm["future-vaccination"].as<double>());

  if (vm.count("vary-vaccination"))
    sim.setVaryVaccination();

  if (vm.count("select"))
    {
      if (vm["select"].as<std::string>() == "all")
        for (size_t i = 0; i < model.cities.size(); ++i)
          selectCities.push_back(i);
      else
        {
          // Tokenizer for reading a comma-separated list of cities to record
          Tokenizer tok(vm["select"].as<std::string>(), sep);
          std::vector<std::string> lineVector;
          lineVector.assign(tok.begin(), tok.end());

          for (std::vector<std::string>::iterator it = lineVector.begin();
               it != lineVector.end(); it++)
            {
              size_t i = 0;
              while (model.cities[i].name != (*it) && i < model.cities.size())
                ++i;
              if (i < model.cities.size())
                selectCities.push_back(i);
            }
        }
    }

  //! Vector of target correlations. If we're trying
  // to match the data we compare the simulated correlations (all
  // cities with with London) to what's been measured
  std::vector<double> targetCorrelations;

  // Write header
  if (out)
    {
      std::stringstream header;
      header << "run,time";
      for (std::vector<size_t>::iterator it = selectCities.begin();
           it != selectCities.end(); it++)
        {
          header << ",S." << model.cities[*it].name << ",I."
                 << model.cities[*it].name;
        }
      header << std::endl;

      out << header.str();
    }

  //! Then it runs the simulations.
  std::stringstream stream; //! output stream

  for (size_t current_sim = 0; current_sim < nsim; ++current_sim)
    {
      if (verbose)
        std::cout << "Sim #" << current_sim << std::endl;

      sim.run(seed, simNb > 0 ? simNb : current_sim );

      if (out)
        out << sim;
    }

  if (of.is_open())
    of.close();

  if (vm.count("summary"))
    {
      // split comma-separated list of summary statistics to
      // calculate into vector of strings
      Tokenizer tok(vm["summary"].as<std::string>(), sep);
      std::vector<std::string> stats;
      stats.assign(tok.begin(), tok.end());

      // calculate summary statistics
      std::vector<double> summary_stats;
      sim.getSummaryStats(stats, summary_stats);

      // if we're wrapped inside an R function, write summary
      // statistics to output file
      if (wrapped) {
        std::stringstream fileName;

        // add inputNb to filename if >0
        if (vm["number"].as<unsigned int>() > 0)
          fileName << "output" << vm["number"].as<unsigned int>();
        else
          fileName << "output";

        std::ofstream output(fileName.str().c_str());
        bool first = true;
        for (std::vector<double>::iterator it = summary_stats.begin();
             it != summary_stats.end(); it++)
          {
            {
              if (first)
                first = false;
              else
                output << " ";
              output << *it;
              std::cout << *it << std::endl;
            }
          }
        output << std::endl;
        output.close();
      }

      // if we're verbose or not wrapped inside an R function, print
      // summary statistics
      if (!wrapped || verbose)
        {
          if (summary_stats.size() > 0)
            {
              std::cout << "Summary stats: ";
              for (size_t i = 0; i < summary_stats.size(); ++i)
                {
                  std::cout << summary_stats[i] << " ";
                }
              std::cout << std::endl;
            }
        }
    }
}
