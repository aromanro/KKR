#include "KKRApp.h"
#include "KKRFrame.h"


wxIMPLEMENT_APP(KKRApp);


bool KKRApp::OnInit()
{
	if (!wxApp::OnInit())
		return false;

	frame = new KKRFrame("KKR", wxPoint(50, 50), wxSize(1024, 800));
	frame->Show(true);

	return true;
}

