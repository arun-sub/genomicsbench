//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <iostream>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <cmath>
#include <execinfo.h>

#include "vertex_index.h"
#include "sequence_container.h"
#include "config.h"
#include "logger.h"
#include "utils.h"
#include "memory_info.h"
#include "time.h"
#include "sys/time.h"

#include <getopt.h>

// #define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

bool parseArgs(int argc, char** argv, std::string& readsFasta, 
			   std::string& logFile,
			   int& kmerSize, bool& debug, size_t& numThreads, int& minOverlap, 
			   std::string& configPath, int& minReadLength, bool& unevenCov,
			   bool& highmem)
{
	auto printUsage = []()
	{
		std::cerr << "Usage: kmer-cnt "
				  << " --reads path --out-asm path --config path [--genome-size size]\n"
				  << "\t\t[--min-read length] [--log path] [--treads num] [--extra-params]\n"
				  << "\t\t[--kmer size] [--meta] [--min-ovlp size] [--debug] [-h]\n\n"
				  << "Required arguments:\n"
				  << "  --reads path\tcomma-separated list of read files\n"
				  << "  --config path\tpath to the config file\n\n"
				  << "Optional arguments:\n"
				  << "  --kmer size\tk-mer size [default = 15] \n"
				  << "  --min-ovlp size\tminimum overlap between reads "
				  << "[default = 5000] \n"
				  << "  --debug \t\tenable debug output "
				  << "[default = false] \n"
				  << "  --log log_file\toutput log to file "
				  << "[default = not set] \n"
				  << "  --threads num_threads\tnumber of parallel threads "
				  << "[default = 1] \n"
				  << "  --highmem \t\tuse more memory for faster kmer counting "
				  << "[default = false] \n";
	};
	
	int optionIndex = 0;
	static option longOptions[] =
	{
		{"reads", required_argument, 0, 0},
		{"config", required_argument, 0, 0},
		{"min-read", required_argument, 0, 0},
		{"log", required_argument, 0, 0},
		{"threads", required_argument, 0, 0},
		{"kmer", required_argument, 0, 0},
		{"min-ovlp", required_argument, 0, 0},
		{"debug", no_argument, 0, 0},
		{"highmem", no_argument, 0, 0},
		{0, 0, 0, 0}
	};

	int opt = 0;
	while ((opt = getopt_long(argc, argv, "h", longOptions, &optionIndex)) != -1)
	{
		switch(opt)
		{
		case 0:
			if (!strcmp(longOptions[optionIndex].name, "kmer"))
				kmerSize = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-read"))
				minReadLength = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "threads"))
				numThreads = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "min-ovlp"))
				minOverlap = atoi(optarg);
			else if (!strcmp(longOptions[optionIndex].name, "log"))
				logFile = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "debug"))
				debug = true;
			else if (!strcmp(longOptions[optionIndex].name, "meta"))
				unevenCov = true;
			else if (!strcmp(longOptions[optionIndex].name, "reads"))
				readsFasta = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "config"))
				configPath = optarg;
			else if (!strcmp(longOptions[optionIndex].name, "highmem"))
				highmem = true;
			break;

		case 'h':
			printUsage();
			exit(0);
		}
	}
	if (readsFasta.empty() ||
		configPath.empty())
	{
		printUsage();
		return false;
	}

	return true;
}
const char *__parsec_roi_begin(const char *s, int *beg, int *end)
{
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = 0x7fffffff;
        return s + strlen(s);
    }
    return NULL;
}

const char *__parsec_roi_end(const char *s, int *beg, int *end)
{
    const char *colon = strrchr(s, ':');
    if (colon == NULL) {
        *beg = 0; *end = 0x7fffffff;
        return s + strlen(s);
    }
    return NULL;
}
int main(int argc, char** argv)
{
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	#ifdef NDEBUG
	signal(SIGSEGV, segfaultHandler);
	std::set_terminate(exceptionHandler);
	#endif

	int kmerSize = -1;
	int minReadLength = 0;
	size_t genomeSize = 0;
	int minOverlap = 5000;
	bool debugging = false;
	bool unevenCov = false;
	size_t numThreads = 1;
	std::string readsFasta;
	std::string logFile;
	std::string configPath;
	bool highmem = false;

	if (!parseArgs(argc, argv, readsFasta, logFile,
				   kmerSize, debugging, numThreads, minOverlap, configPath, 
				   minReadLength, unevenCov, highmem)) return 1;

	Logger::get().setDebugging(debugging);
	if (!logFile.empty()) Logger::get().setOutputFile(logFile);
	Logger::get().debug() << "Build date: " << __DATE__ << " " << __TIME__;
	std::ios::sync_with_stdio(false);

	Logger::get().debug() << "Total RAM: " 
		<< getMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Available RAM: " 
		<< getFreeMemorySize() / 1024 / 1024 / 1024 << " Gb";
	Logger::get().debug() << "Total CPUs: " << std::thread::hardware_concurrency();

	Config::load(configPath);

	if (kmerSize == -1)
	{
		kmerSize = Config::get("kmer_size");
	}
	Parameters::get().numThreads = numThreads;
	Parameters::get().kmerSize = kmerSize;
	Parameters::get().minimumOverlap = minOverlap;
	Parameters::get().unevenCoverage = unevenCov;
	Logger::get().debug() << "Running with k-mer size: " << 
		Parameters::get().kmerSize; 
	// Logger::get().debug() << "Running with minimum overlap " << minOverlap;
	// Logger::get().debug() << "Metagenome mode: " << "NY"[unevenCov];

	//TODO: unify minimumOverlap ad safeOverlap concepts
	Parameters::get().minimumOverlap = 1000;

	SequenceContainer readsContainer;
	std::vector<std::string> readsList = splitString(readsFasta, ',');
	Logger::get().info() << "Reading sequences";
	try
	{
		//only use reads that are longer than minOverlap,
		//or a specified threshold (used for downsampling)
		minReadLength = std::max(minReadLength, minOverlap);
		for (auto& readsFile : readsList)
		{
			readsContainer.loadFromFile(readsFile, minReadLength);
		}
	}
	catch (SequenceContainer::ParseException& e)
	{
		Logger::get().error() << e.what();
		return 1;
	}
	readsContainer.buildPositionIndex();
	VertexIndex vertexIndex(readsContainer,
							(int)Config::get("assemble_kmer_sample"),
							highmem);
	vertexIndex.outputProgress(true);

	/*int64_t sumLength = 0;
	for (auto& seq : readsContainer.iterSeqs())
	{
		sumLength += seq.sequence.length();
	}
	int coverage = sumLength / 2 / genomeSize;
	Logger::get().debug() << "Expected read coverage: " << coverage;*/

	// const int MIN_FREQ = 2;
	// static const float SELECT_RATE = Config::get("meta_read_top_kmer_rate");
	// static const int TANDEM_FREQ = Config::get("meta_read_filter_kmer_freq");

	//Building index
	struct timeval start_time, end_time;
	double runtime = 0;
	gettimeofday(&start_time, NULL);
	const char *roi_q;
	int roi_i, roi_j;
	char roi_s[20] = "chr22:0-5";
#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
	roi_q = __parsec_roi_begin(roi_s, &roi_i, &roi_j);
	bool useMinimizers = Config::get("use_minimizers");
	if (useMinimizers)
	{
		const int minWnd = Config::get("minimizer_window");
		vertexIndex.buildIndexMinimizers(/*min freq*/ 1, minWnd);
	}
	else	//indexing using solid k-mers
	{
		vertexIndex.countKmers();
		// vertexIndex.buildIndexUnevenCoverage(MIN_FREQ, SELECT_RATE, 
		//									 TANDEM_FREQ);
	}
	roi_q = __parsec_roi_end(roi_s, &roi_i, &roi_j);
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	gettimeofday(&end_time, NULL);
	runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;
	Logger::get().debug() << "Peak RAM usage: " 
		<< getPeakRSS() / 1024 / 1024 / 1024 << " Gb";
	fprintf(stderr, "Kernel time: %.3f sec\n", runtime * 1e-6);
	return 0;
}
