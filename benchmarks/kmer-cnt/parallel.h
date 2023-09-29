//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include <functional>
#include <atomic>
#include <thread>

#include "progress_bar.h"

//simple thread pool implementation
//updateFun should be thread-safe!
template <class T>
void processInParallel(const std::vector<T>& scheduledTasks,
					   std::function<void(const T&)> updateFun,
					   size_t maxThreads, bool progressBar)
{
	if (scheduledTasks.empty()) return;

	std::atomic<size_t> jobId(0);
	ProgressPercent progress(scheduledTasks.size());
	if (progressBar) progress.advance(0);

	#pragma omp parallel for
	for (size_t i = 0; i < std::min(maxThreads, scheduledTasks.size()); ++i)
	{
		bool finished = false;
		while (!finished)
		{
			size_t expected = 0;
			while(true)
			{
				expected = jobId;
				if (jobId == scheduledTasks.size()) 
				{
					finished = true;
					break;
				}
				if (jobId.compare_exchange_weak(expected, expected + 1))
				{
					break;
				}
			}
			if (!finished) {
				updateFun(scheduledTasks[expected]);
				if (progressBar) progress.advance();
			}
		}
	}
}

