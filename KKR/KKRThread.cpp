#include "KKRThread.h"

#include "KKRFrame.h"

KKRThread::KKRThread(const Options& options, KKRFrame* frame)
	: m_options(options), m_frame(frame), terminate(false)
{
}


KKRThread::~KKRThread()
{
	join();
}

void KKRThread::Start()
{
	mThread = std::thread([this]() {
		Calculate();
		});
}

void KKRThread::join()
{
	if (mThread.joinable()) mThread.join();
}


void KKRThread::Calculate()
{
	++m_frame->runningThreads;

	results = m_frame->bandStructure.Compute(terminate, m_options);

	--m_frame->runningThreads;
}

void KKRThread::Terminate()
{
	terminate = true;
}
