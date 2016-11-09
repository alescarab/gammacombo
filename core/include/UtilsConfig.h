/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: June 2014
 *
 **/

#ifndef UtilsConfig_h
#define UtilsConfig_h
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <cassert>

#include "TString.h"

using namespace std;

namespace Utils
{
	enum config {
		year2014,
		babar,
		babar2007,
		babar2008,
		babar2010,
		babar2012,
		babar_dpi0,
		babar_dg,
		belle,
		belle2005cleo2009,
		belle2006,
		belle2007,
		belle2009,
		belle2012,
		belle2012preliminary,
		belle2013,
		belle2014,
		belle2013cleo2014,
		belle_dpi0,
		belle_dg,
    belle_md,
    belle_1ab,
    belle_5ab,
    belle_50ab,
    belle2_1ab,
    belle2_5ab,
    belle2_50ab,
    belle2_1ab_had,
    belle2_5ab_had,
    belle2_50ab_had,
    belle2_1ab_sl,
    belle2_5ab_sl,
    belle2_50ab_sl,
    belle2_1ab_WA,
    belle2_5ab_WA,
    belle2_50ab_WA,
    belle2_1ab_Dstar,
    belle2_5ab_Dstar,
    belle2_50ab_Dstar,
    belle2_1ab_D,
    belle2_5ab_D,
    belle2_50ab_D,
    belle2_1ab_tagged,
    belle2_5ab_tagged,
    belle2_50ab_tagged,
    belle2_1ab_untagged,
    belle2_5ab_untagged,
    belle2_50ab_untagged,
    cdf,
		cdf2007,
		cdf2012,
		cdf2013,
    check,
		ckm2014,
		ckm2015,
		cleo,
		cleo2001,
		cleo2009,
		cleo2012,
		cleo2014,
		cleo2015,
    cleo2016,
		cleoFullDP,
		combos2008,
    combpipi,
		default_config,
		excludeKdDdK3pi,
		exclusive2014,
		inclusive2014,
		focus2000,
		hfag,
		hfagFP2014,
		hfagLP2011,
		hfagCHARM2015,
		highrb,
		highstattoy,
		lhcb,
		lhcb_WA,
    lhcb_md,
    lhcb_mi,
		lhcb2011,
		lhcb2012,
		lhcb2013,
		lhcb2013KK,
		lhcb2013pipi,
		lhcb2013preliminary,
		lhcb2014,
		lhcb2018KK_extrap,
		lhcb_upgrade_extrap,
		lhcbphis,
		lhcbcomb,
    lhcb_old,
    lhcb_pipi,
    lhcb_kk,
    lhcb_3fb,
    lhcb_10fb,
    lhcb_22fb,
    lhcb_3fb_WA,
    lhcb_10fb_WA,
    lhcb_22fb_WA,
    lumi40pb,
		lumi1fb,
		lumi1fbConfcFit,
		lumi1fbConfsFit,
		lumi1fbNoAfav,
		lumi1fbPapercFit,
		lumi1fbPapercFitExpected,
		lumi1fbPapersFit,
		lumi1fbSystCor,
		lumi1fbprompt,
		lumi1fbsl,
		lumi2fb,
		lumi3fb,
		lumi3fbCPVA,
		lumi3fbFix,
		lumi3fbFullDP,
		lumi3fbPaper,
		lumi3fbDKstz,
    lumi3fb_estimate,
		lumi5ab,
		lumi5fb,
		lumi9fb,
		lumi50ab,
		lumi50fb,
		lambda1_3fb,
		lambdafree_3fb,
		manual,
    milc_update,
		none,
		nophicorr,
    onlyGsDGs,
    pdg,
		sneha,
    SM_prediction,
    statonly,
    test,
    test1,
    test2,
    test3,
    test4,
    test5,
    test6,
    test7,
    test8,
    test9,
		toy,
		truth,
		useBicubic,
		useCartCoords,
		useGaussian,
		useHistogram,
		useParametric,
		usePolarCoords,
		useTradObs,
    ut2014,
    world_average,
		zero
	};

	config  TStringToConfig(TString s);
	TString ConfigToTString(config c);
}

#endif
