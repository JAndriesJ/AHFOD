load('E:\GENEVA\Processed Data\Geneva\PAT_2012\255619\Data\HFOSummary\HFOSummaryMat.mat')
hfosummary1 = HFOSummaryMat;

load('E:\GENEVA\Processed Data\Geneva\PAT_2012\257707\Data\HFOSummary\HFOSummaryMat.mat')
hfosummary2 = HFOSummaryMat;


TEMPhfosummary.Ripples              = [hfosummary1.Ripples    , hfosummary2.Ripples];
TEMPhfosummary.FastRipples          = [hfosummary1.FastRipples, hfosummary2.FastRipples];  
TEMPhfosummary.RandFR               = [hfosummary1.RandFR     , hfosummary2.RandFR];

TEMPhfosummary.Rates.RHFO           = [hfosummary1.Rates.RHFO     ;  hfosummary2.Rates.RHFO];
TEMPhfosummary.Rates.FRHFO          = [hfosummary1.Rates.FRHFO    ;  hfosummary2.Rates.FRHFO];  
TEMPhfosummary.Rates.RandFRHFO      = [hfosummary1.Rates.RandFRHFO;  hfosummary2.Rates.RandFRHFO];

TEMPhfosummary.Noise.Ripples        = [hfosummary1.Noise.Ripples     ; hfosummary2.Noise.Ripples];
TEMPhfosummary.Noise.FastRipples    = [hfosummary1.Noise.FastRipples ; hfosummary2.Noise.FastRipples];

TEMPhfosummary.Baseline.Ripples     = [hfosummary1.Baseline.Ripples      ; hfosummary2.Baseline.Ripples ];
TEMPhfosummary.Baseline.FastRipples = [hfosummary1.Baseline.FastRipples  ; hfosummary2.Baseline.FastRipples ];


HFOSummaryMat = TEMPhfosummary;