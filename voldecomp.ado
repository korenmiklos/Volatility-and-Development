/* check GSECT */
program define voldecomp
version 9.0
syntax , shock(varname) share(varname) sector(varname) country(varlist) time(varname) [global(varname)] [gdp(varname)] [smooth(integer 1)]
/*  This program calculates the volatility components used in
    Koren and Tenreyro, "Volatility and Development", Quarterly Journal of Economics, 122(1):243-287, February 2007.
    Please cite this paper if you use it.

    It takes as input a long-shaped data of countries, sectors and years (or other time variable). The panel should be
    balanced for sectors: if "shock" or "share" data (see below) is missing for a sector in a country, that country is
    dropped from the analysis.

    Arguments:
    shock:  the variable containing the shock for which we will calculate the volatility
            for example, the growth of labor productivity in the sector
    share:  the variable containing the share with which we can aggregate up the shocks
            for example, employment in the sector
    sector: the variable indexing the sector
    country:the variable indexing the country
    time:   the variable indexing time
    global: (optional) the variable containing global sectoral shocks.
            If not specified, this will be calculated as the simple average of the "shock" for all countries.
    gdp:    (optional) the variable containing real GDP per capita so that it can be saved alongside the volatility measures
    smooth: (optional) an integer specifying the amount of smoothing.
            If "smooth"=1 (default), then "shock" is not smoothed, if "smooth"=5, the it takes
            5-year moving average of "shock"

    Output:
    A long-shaped dataset of countries and years with each of the risk measured as defined in "Volatility and Development".
    Because this program creates a new dataset, make sure to save your old one!

    Example:
    voldecomp, shock(lpgrowth) share(employment) sector(isic) country(uncode) time(year)
    */


/* these are temporary variables for the program */
tempvar cnt idio herf sumshare newglobal group tag mX mY XY N mean mean1 mean2 mean3 NNsect
tempname Sigma Nsect sigma2
tempfile widefile

/* if global shocks are not saved, create them first. */
if ("`global'"=="") {
    qui egen `newglobal' = mean(cond(`share'>0,`shock',.)), by(`sector' `time')
    local global `newglobal'
}

/* do the TSSET */
qui egen `group' = group(`country' `sector')
qui drop if missing(`group')
qui tsset `group' `time'

/* if required, perform some MA smoothing */
if ("`smooth'"!="1") {
    tssmooth ma `shock' = `shock', window(`smooth') replace
    tssmooth ma `global' = `global', window(`smooth') replace
}

/* save # of sectors */
drop `group'
egen `group' = group(`sector')
qui su `group'
local `Nsect'  `r(max)'

/* drop country*years with unbalanced sectors */
egen `NNsect' = sum(!missing(`shock')), by(`country' `time')
su `NNsect', d
replace `shock'=. if `NNsect'<``Nsect''
drop `NNsect'
egen `NNsect' = sum(!missing(`share')), by(`country' `time')
su `NNsect', d
replace `share'=. if `NNsect'<``Nsect''

/* make sure shares sum to 1 */
egen `sumshare' = sum(`share'), by(`country' `time')
qui su `sumshare'
if (r(max)>1)|(r(min)<1) {
    di in red "Something wrong with the share variable!"
    di in red "Rescaling by sum of shares."
    replace `share'=`share'/`sumshare'
}

/* reshape the data wide */
preserve
    keep `country' `group' `time' `global' `share'

    qui reshape wide `share' `global', i(`country' `time') j(`group')
    egen `tag' = tag(`time')

    /* do the global covariance matrix*/
    matrix accum Sigma = `global'* if `tag', nocons dev
    matrix Sigma = Sigma/r(N)
    matname Sigma `share'*, col(.)

    /* this are the columns of alpha'Sigma alpha */
    forval i=1/``Nsect'' {
        matrix nicelittlevector=Sigma["`global'`i'","`share'1".."`share'``Nsect''"]
        matrix score Sigmaalpha`i'=nicelittlevector
        qui replace Sigmaalpha`i'=Sigmaalpha`i'*`share'`i'
    }

    /* save the matrix for future work */

    /* now back to long format */
    keep `country' `time' Sigmaalpha*
    qui reshape long Sigmaalpha, i(`country' `time') j(`group')
    sort `country' `group' `time'
    save `widefile'
restore

/* now merge with the saved global risk */
capture drop _merge
sort `country' `group' `time'
merge `country' `group' `time' using `widefile', nokeep
capture drop _merge

/* create country and idio shocks */
egen `cnt' = mean(`shock'-`global'), by(`country' `time')
gen `idio' = `shock'-`global'-`cnt'

qui su `idio'
local `sigma2' `r(Var)'

/* calculate the variances. FORGET ABOUT rolling */
drop `group'
egen `group' = group(`country' `sector')
qui tsset `group' `time'

egen `N' = sum(!missing(`idio',`share')), by(`country' `sector')
egen `cnt'_var=sd(`cnt'), by(`country' `sector')
egen `idio'_var=sd(`idio'), by(`country' `sector')
replace `cnt'_var=`cnt'_var^2
replace `idio'_var=`idio'_var^2

egen `mX' = mean(cond(!missing(`cnt',`global'),`cnt',.)), by(`country' `sector')
egen `mY' = mean(cond(!missing(`cnt',`global'),`global',.)), by(`country' `sector')
egen `global'_cov = mean(cond(!missing(`cnt',`global'),(`cnt'-`mX')*(`global'-`mY'),.)), by(`country' `sector')

/* now aggregate up to risk components */
qui replace `idio'_var=`share'^2*`idio'_var
qui replace `global'_cov = 2*`share'*`global'_cov
gen `herf' = `share'^2


/* collapse to 1 obs per country*year. also save N */
if ("`gdp'"=="") {
    collapse (mean) idio = `idio'_var herf=`herf' cov = `global'_cov gsect = Sigmaalpha (mean) cnt = `cnt'_var (min) N=`N', by(`country' `time')
}
else {
    collapse (mean) idio = `idio'_var herf=`herf' cov = `global'_cov gsect = Sigmaalpha (mean) gdp =`gdp' cnt = `cnt'_var (min) N=`N', by(`country' `time')
}

foreach XX of var idio herf cov gsect {
    replace `XX'=`XX'*``Nsect''
}

gen isect = idio/herf

/* label variables appropriately */
label var idio "Idiosyncratic risk"
label var herf "Herfindahl index"
label var cov "Covariance"
label var cnt "Country risk"
label var gsect "Global sectoral risk"
label var isect "Average idiosyncratic variance"
label var N "Number of years in sample"
if ("`gdp'"!="") {
    label var gdp "GDP (log)"
}

tsset, clear

end
