#pragma once

#include <string>
#include <vector>

#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/fileconf.h>

class Options
{
public:
	Options();

	void Load();
	void Save();

	void Open();
	void Close();

	int nrThreads;

	int nrPoints;

	int pathNo;

	std::vector<std::vector<std::string>> paths;

protected:
	wxFileConfig *m_fileconfig;
};

