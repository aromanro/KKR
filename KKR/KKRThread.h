#pragma once

#include <thread>
#include <vector>
#include <atomic>

#include "Options.h"

class KKRFrame;

class KKRThread
{
public:
	KKRThread(const Options& options, KKRFrame* frame);
	~KKRThread();


	void Start();
	void join();

	void Terminate();


	const Options& m_options;

	KKRFrame* m_frame;

	std::vector<std::vector<double>> results;

protected:
	void Calculate();

	std::thread mThread;
	std::atomic_bool terminate;
};

