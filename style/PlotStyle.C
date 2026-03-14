{
	TStyle * style = new TStyle("minerva", "MINERvA publication plot style");

	// First, copy the stuff from the ROOT "Plain" style.
	style->SetCanvasBorderMode(0);
	style->SetPadBorderMode(0);
	style->SetPadColor(0);
	style->SetCanvasColor(0);
	style->SetTitleColor(1);
	style->SetStatColor(0);

	// Set the size of the default canvas: 600x500 looks almost square.
	style->SetCanvasDefH(500);
	style->SetCanvasDefW(600);
	style->SetCanvasDefX(10);
	style->SetCanvasDefY(10);

	// Color Scheme - Black & White
	style->SetPalette(1);
	style->SetFrameBorderMode(0);

	// Line Widths
	style->SetFrameLineWidth(1);
	style->SetLineWidth(1);
	style->SetHistLineWidth(2);
	style->SetLegendBorderSize(1);

	// Marker Styles
	style->SetMarkerStyle(20);

	// Stats
	style->SetOptStat(0000);
	style->SetOptFit(0000);

	// Margins
	style->SetPadTopMargin(0.09);
	style->SetPadBottomMargin(0.15);
	style->SetPadLeftMargin(0.15);
	style->SetPadRightMargin(0.15);

	// Titles
	style->SetOptTitle(0);  // titles are off by default...
	// ...but if they are turned on, they shouldn't look hideous.
	style->SetTitleBorderSize(0);
	style->SetTitleFillColor(kWhite);
	style->SetTitleAlign(23);    // centered
	style->SetTitleX(0.5);       // the center should be in the middle
	style->SetTitleW(0.8);       // don't use the full width so that the centering is visible

	// Errors
	style->SetEndErrorSize(3);
	style->SetErrorX(0.5);

	style->SetAxisColor(1, "XYZ");
    style->SetTickLength(0.03, "XYZ");
    style->SetNdivisions(510, "XYZ");
    style->SetPadTickX(1);
    style->SetPadTickY(1);

	// Finally...
	gROOT->SetStyle("minerva");
	gROOT->ForceStyle();
}
