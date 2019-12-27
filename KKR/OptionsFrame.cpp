
#define wxNEEDS_DECL_BEFORE_TEMPLATE

#include <wx/statline.h>
#include <wx/valgen.h>
#include <wx/valnum.h>


#include "KKRApp.h"
#include "OptionsFrame.h"
#include "KKRFrame.h"

#define ID_NRTHREADS 101
#define ID_NRPOINTS 103
#define ID_PATH 104

wxDECLARE_APP(KKRApp);

OptionsFrame::OptionsFrame(const wxString& title, wxWindow* parent)
	: wxDialog(parent, wxID_ANY, title, wxDefaultPosition, wxSize(400, 200))
{
	CreateControls();

	Centre();
}

void OptionsFrame::CreateControls()
{
	// box to contain them all
	wxBoxSizer* vbox = new wxBoxSizer(wxVERTICAL);
	SetSizer(vbox);

	// box with margin to contain option controls
	wxBoxSizer* boxSizer = new wxBoxSizer(wxVERTICAL);
	vbox->Add(boxSizer, 0, wxALIGN_CENTER_HORIZONTAL | wxGROW | wxALL, 5);

	// ******************************************************************
	// now option controls

	// nr threads 

	wxBoxSizer* box = new wxBoxSizer(wxHORIZONTAL);
	boxSizer->Add(box, 0, wxGROW | wxALL, 5);

	wxStaticText* label = new wxStaticText(this, wxID_STATIC, "&Threads:", wxDefaultPosition, wxSize(60, -1), wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL | wxALL, 5);

	wxString str = wxString::Format(wxT("%i"), options.nrThreads);
	wxTextCtrl* nrThreadsCtrl = new wxTextCtrl(this, ID_NRTHREADS, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(nrThreadsCtrl, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);


	box->Add(5, 5, 1, wxALIGN_CENTER_VERTICAL | wxALL, 5); // pushes to the right

	// next line

	box = new wxBoxSizer(wxHORIZONTAL);
	boxSizer->Add(box, 0, wxGROW | wxALL, 5);

	// nr points

	label = new wxStaticText(this, wxID_STATIC, "Nr. &Points:", wxDefaultPosition, wxSize(60, -1), wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL | wxALL, 5);

	str = wxString::Format(wxT("%i"), options.nrPoints);
	wxTextCtrl* nrPointsCtrl = new wxTextCtrl(this, ID_NRPOINTS, str, wxDefaultPosition, wxSize(60, -1), 0);
	box->Add(nrPointsCtrl, 0, wxGROW | wxALIGN_CENTER_VERTICAL | wxALL, 5);


	box->Add(5, 5, 1, wxALIGN_CENTER_VERTICAL | wxALL, 5); // pushes to the right


	// path

	label = new wxStaticText(this, wxID_STATIC, "&Path:", wxDefaultPosition, wxDefaultSize, wxALIGN_RIGHT);
	box->Add(label, 0, wxALIGN_RIGHT | wxALIGN_CENTER_VERTICAL | wxALL, 5);

	wxString* pathStrings = new wxString[options.paths.size()];

	for (int i = 0; i < options.paths.size(); ++i)
	{
		wxString p;
		for (int j = 0; j < options.paths[i].size(); ++j)
		{
			p += options.paths[i][j];
			if (j < options.paths[i].size() - 1) p += "->";
		}

		pathStrings[i] = p;
	}

	wxChoice* pathChoice = new wxChoice(this, ID_PATH, wxDefaultPosition, wxSize(160, -1), options.paths.size(), pathStrings, 0);
	pathChoice->SetSelection(options.pathNo);
	box->Add(pathChoice, 0, wxALIGN_CENTER_VERTICAL | wxALL, 5);

	delete[] pathStrings;

	// ******************************************************************
	// setting validators

	// there is a bug in validator ranges, so set the range from 1
	// and validate for the minimum value elsewhere
	wxIntegerValidator<int> val1(&options.nrThreads, wxNUM_VAL_DEFAULT);
	val1.SetRange(1, 64);
	nrThreadsCtrl->SetValidator(val1);

	wxIntegerValidator<int> val2(&options.nrPoints, wxNUM_VAL_DEFAULT);
	//val2.SetMin(100);
	val2.SetRange(1, 1000);
	nrPointsCtrl->SetValidator(val2);


	pathChoice->SetValidator(wxGenericValidator(&options.pathNo));


	// ******************************************************************

	// divider line

	wxStaticLine* line = new wxStaticLine(this, wxID_STATIC, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL);
	boxSizer->Add(line, 0, wxGROW | wxALL, 5);

	// bottom box with ok & cancel buttons

	wxBoxSizer* hbox = new wxBoxSizer(wxHORIZONTAL);


	wxButton* okButton = new wxButton(this, wxID_OK, "Ok", wxDefaultPosition, wxSize(70, 30));
	wxButton* closeButton = new wxButton(this, wxID_CANCEL, "Cancel", wxDefaultPosition, wxSize(70, 30));

	hbox->Add(okButton, 1);
	hbox->Add(closeButton, 1, wxLEFT, 5);

	vbox->Add(hbox, 0, wxALIGN_CENTER | wxTOP | wxBOTTOM, 10);
}


bool OptionsFrame::TransferDataFromWindow()
{
	if (!wxDialog::TransferDataFromWindow()) return false;

	wxTextCtrl* nrPointsCtrl = (wxTextCtrl*)FindWindow(ID_NRPOINTS);
	wxString str = nrPointsCtrl->GetValue();
	long int val = 0;
	if (!str.ToLong(&val)) return false;
	if (val < 100 || val > 1000)
	{
		wxMessageBox("Please enter between 100 and 1000 points", "Validation", wxOK | wxICON_INFORMATION, this);

		return false;
	}

	return true;
}

