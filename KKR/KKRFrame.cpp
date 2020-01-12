#include "KKRFrame.h"

#include "OptionsFrame.h"

#include "wx/aboutdlg.h"
#include "wx/statline.h"
#include "wx/generic/aboutdlgg.h"

#include <vtkAutoInit.h>


VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);

#define MY_VTK_WINDOW 102

#define ID_CALCULATE 105

wxBEGIN_EVENT_TABLE(KKRFrame, wxFrame)
EVT_MENU(ID_CALCULATE, KKRFrame::OnCalculate)
EVT_UPDATE_UI(ID_CALCULATE, KKRFrame::OnUpdateCalculate)
EVT_MENU(wxID_EXIT, KKRFrame::OnExit)
EVT_CLOSE(KKRFrame::OnClose)
EVT_MENU(wxID_PREFERENCES, KKRFrame::OnOptions)
EVT_MENU(wxID_ABOUT, KKRFrame::OnAbout)
EVT_TIMER(101, KKRFrame::OnTimer)
EVT_ERASE_BACKGROUND(KKRFrame::OnEraseBackground)
wxEND_EVENT_TABLE()


KKRFrame::KKRFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
	: wxFrame(NULL, wxID_ANY, title, pos, size),
	timer(this, 101),
	theThread(computeOptions, this),
	runningThreads(0)
{
	wxMenu* menuFile = new wxMenu;

	menuFile->Append(ID_CALCULATE, "C&alculate\tCtrl+a", "Starts computing");
	menuFile->Append(wxID_SEPARATOR);
	menuFile->Append(wxID_EXIT);

	wxMenu* menuView = new wxMenu;
	menuView->Append(wxID_PREFERENCES);

	wxMenu* menuHelp = new wxMenu;
	menuHelp->Append(wxID_ABOUT);

	wxMenuBar* menuBar = new wxMenuBar;
	menuBar->Append(menuFile, "&File");
	menuBar->Append(menuView, "&View");
	menuBar->Append(menuHelp, "&Help");

	SetMenuBar(menuBar);

	CreateStatusBar();
	SetStatusText("Welcome to KKR!");

	m_pVTKWindow = new wxVTKRenderWindowInteractor(this, MY_VTK_WINDOW);
	m_pVTKWindow->UseCaptureMouseOn();
	//m_pVTKWindow->DebugOn();
	m_pVTKWindow->DebugOff();

	ConstructVTK();

	std::vector<std::vector<double>> empty_results;
	std::vector<unsigned int> empty_pos;
	std::vector<std::string> empty_strings;
	ConfigureVTK("", empty_results, empty_pos, empty_strings);

	currentOptions.Load();
}


KKRFrame::~KKRFrame()
{
	DestroyVTK();
	if (m_pVTKWindow) m_pVTKWindow->Delete();
}

void KKRFrame::ConstructVTK()
{
	pRenderer = vtkRenderer::New();
	pContextView = vtkContextView::New();

	vtkRenderWindow* pRenderWindow = m_pVTKWindow->GetRenderWindow();
	pRenderWindow->AddRenderer(pRenderer);
	pContextView->SetInteractor(pRenderWindow->GetInteractor());
	//pContextView->GetInteractor()->Initialize();

	pChart = vtkChartXY::New();
	pChart->SetRenderEmpty(true);
	pContextView->GetScene()->AddItem(pChart);
}

void KKRFrame::DestroyVTK()
{
	if (pChart) pChart->Delete();
	if (pRenderer) pRenderer->Delete();
	if (pContextView) pContextView->Delete();
}


void KKRFrame::OnCalculate(wxCommandEvent& /*event*/)
{
	Compute();
}

void KKRFrame::OnUpdateCalculate(wxUpdateUIEvent& event)
{
	event.Enable(isFinished());
}


void KKRFrame::ConfigureVTK(const std::string& name, const std::vector<std::vector<double>>& results, std::vector<unsigned int>& symmetryPointsPositions, const std::vector<std::string>& symmetryPointsLabels)
{
	pChart->ClearPlots();

	if (name.empty()) pChart->SetTitle("Band");
	else
	{
		std::string Name = name;
		Name += " Band";
		pChart->SetTitle(Name.c_str());
	}

	pChart->SetAutoAxes(false);

	pChart->GetAxis(vtkAxis::BOTTOM)->SetTitle("k - Symmetry Points");
	pChart->GetAxis(vtkAxis::LEFT)->SetTitle("Energy");

	if (results.empty()) return;


	size_t maxNumPoints = 0;
	for (int i = 0; i < results.size(); ++i)
		maxNumPoints = std::max(results[i].size(), maxNumPoints);

	// set up the data table

	vtkNew<vtkTable> table;

	vtkNew<vtkFloatArray> arrX;
	arrX->SetName("X");
	table->AddColumn(arrX.GetPointer());

	for (int i = 0; i < maxNumPoints; ++i)
	{
		vtkNew<vtkFloatArray> arrY;
		arrY->SetName((std::string("Y") + std::to_string(i)).c_str());
		table->AddColumn(arrY.GetPointer());
	}

	table->SetNumberOfRows(results.size());

	// set values for X axis column
	for (int i = 0; i < results.size(); ++i)
		table->SetValue(i, 0, i);

	for (int i = 0; i < results.size(); ++i)
	{
		for (int j = 0; j < results[i].size(); ++j)
			table->SetValue(i, j + 1ULL, results[i][j]);

		for (int j = results[i].size(); j < maxNumPoints; ++j)
			table->SetValue(i, j + 1ULL, results[i].empty() ? 0 : results[i].back());
	}


	for (int i = 0; i < maxNumPoints; ++i)
	{
		vtkPlot* points = pChart->AddPlot(vtkChart::POINTS);
		points->SetInputData(table.GetPointer(), 0, i + 1LL);
		points->SetColor(0, 0, 255, 255);
		points->SetWidth(1.0);
		dynamic_cast<vtkPlotPoints*>(points)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
	}

	vtkAxis* xAxis = pChart->GetAxis(vtkAxis::BOTTOM);

	vtkNew<vtkDoubleArray> posArray;
	vtkNew<vtkStringArray> labelsArray;

	posArray->SetName("X");
	labelsArray->SetName("X");

	posArray->SetNumberOfValues(symmetryPointsPositions.size() + 1);
	labelsArray->SetNumberOfValues(symmetryPointsPositions.size() + 1);

	int index = 0;
	for (unsigned int pos : symmetryPointsPositions)
	{
		posArray->SetValue(index, pos);
		labelsArray->SetValue(index, symmetryPointsLabels[index]);
		++index;
	}
	// add the last

	posArray->SetValue(index, results.size() - 1ULL);
	labelsArray->SetValue(index, symmetryPointsLabels[index]);

	xAxis->SetCustomTickPositions(posArray.GetPointer(), labelsArray.GetPointer());	
}

bool KKRFrame::isFinished() const
{
	return 0 == runningThreads;
}


void KKRFrame::OnOptions(wxCommandEvent& /*event*/)
{
	OptionsFrame* optionsFrame = new OptionsFrame("Options", this);
	optionsFrame->options = currentOptions;
	if (wxID_OK == optionsFrame->ShowModal())
	{
		currentOptions = optionsFrame->options;
		currentOptions.Save();
	}

	delete optionsFrame;
}


void KKRFrame::Compute()
{
	if (!isFinished()) return;

	wxBeginBusyCursor();

	computeOptions = currentOptions;

	SetTitle("Computing - KKR");

	bandStructure.Initialize(computeOptions.paths[computeOptions.pathNo], computeOptions.nrPoints);


	unsigned int nrThreads = computeOptions.nrThreads;
	if (0 == nrThreads) computeOptions.nrThreads = nrThreads = 1;

	theThread.Start();

	timer.Start(100);
}


void KKRFrame::OnTimer(wxTimerEvent& WXUNUSED(event))
{
	if (isFinished())
	{
		timer.Stop();
		StopThreads();

		SetTitle("Finished - KKR");

		Refresh();
	}
}

void KKRFrame::OnEraseBackground(wxEraseEvent& event)
{
	event.Skip(false);
}


void KKRFrame::OnExit(wxCommandEvent& /*event*/)
{
	StopThreads(true);

	Close(true);
}

void KKRFrame::OnClose(wxCloseEvent& event)
{
	StopThreads(true);

	event.Skip();
}


void KKRFrame::OnAbout(wxCommandEvent& /*event*/)
{
	wxAboutDialogInfo info;

	info.SetName("KKR");

	static const int majorVer = 1;
	static const int minorVer = 0;
	wxString verStr = wxString::Format("%d.%d", majorVer, minorVer);
	info.SetVersion(verStr, wxString::Format("Version %s", verStr));

	info.SetDescription("   KKR Application   ");
	info.SetLicense("GNU GPL v3.0, see LICENSE file for details");

	info.AddDeveloper("Adrian Roman");

	info.SetWebSite("https://github.com/aromanro/KKR", "GitHub repository");


	wxAboutBox(info, this);
}


void KKRFrame::StopThreads(bool cancel)
{
	if (cancel)
		theThread.Terminate();

	theThread.join();

	if (!cancel)
		bandStructure.results = theThread.results;

	SetTitle("Finished - KKR");

	if (!cancel)
		ConfigureVTK("Cu", bandStructure.results, bandStructure.symmetryPointsPositions, bandStructure.GetPath());

	if (wxIsBusy()) wxEndBusyCursor();
}