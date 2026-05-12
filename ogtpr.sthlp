{smcl}
{* Mayo 2026}{...}
{hline}
help for {hi:ogtpr }{right:Dialog box:  {bf:{dialog ogtpr}}}
{hline}

{title: Optimal Targeting Groups for Poverty Reduction} 

{p 8 10}{cmd:ogtpr}  {it:varlist}  {cmd:,} [ 
{cmd:HSize(}{it:varname}{cmd:)} {cmd:HGroup(}{it:varname}{cmd:)} {cmd:PLine(}{it:real}{cmd:)}  {cmd:TRANS(}{it:real}{cmd:)  {cmd:PART(}{it:int}{cmd:)  {cmd:ERED(}{it:int}{cmd:)} 
{cmd:ALpha(}{it:real}{cmd:)} ]
 

{p}where {p_end}
{p 8 8} {cmd:varlist} is the variable of wellbeing (income). {p_end}


{title:Version} 15.0 and higher.

{title:Description}
 {p}{cmd:Poverty: Optimal Targeting by Population Groups}  {p_end}
 {p} Users should set their surveys' sampling design before using this module 
 (and save their data files). If the sampling design is not set, simple-random sampling (SRS) will be automatically assigned 
by default. {cmd:ogtpr} estimates the group-lumpsum transfers  to reduce optimally  the aggregate poverty for a given predefined budget of transfers.

{title:Options}

{p 0 4} {cmdab:hsize}        Variable that captures the size of the household. This variable is used to weight observations by household size (in addition to sampling weights, best set in survey design). {p_end}

{p 0 4} {cmdab:hgroup}       Variable that captures the socio-demographic group to be used in the decomposition. For example, for an urban-rural decomposition of poverty, this variable could equal 1 for rural households and 2 for urban ones. The associated varlist should contain only one variable. {p_end}

{p 0 4} {cmdab:alpha}        To set the FGT parameter (alpha). By default, alpha=0.   {p_end}

{p 0 4} {cmdab:pline}        To set the poverty line. {p_end}

{p 0 4} {cmdab:trans}        To set the fixed per capita lump sum transfer. {p_end}

{p 0 4} {cmdab:ered}         Set the value to one ( ered(1) ) to estimate the reduction in aggregate poverty and the quality of the group indicator. {p_end}

{p 0 4} {cmdab:dec}          To set the number of decimals used in the display of results. {p_end}


{title:Examples}
use bkf98I.dta, clear
ogtpr exppc, hgroup(gse) hsize(size) alpha(0) pline(80000) trans(4000) part(1000) ered(1)


{title:Reference}
{p 4 8}Araar, A. (2026). "Optimal Group Targeting for Poverty Reduction: Bridging Theory and Algorithm." Working Paper, Universite Laval & PEP. {p_end}


{title:Author(s)}
{p 4 8}Abdelkrim Araar, Universite Laval & PEP. {p_end}
{p 4 8}Email: {browse "mailto:aabd@ecn.ulaval.ca":aabd@ecn.ulaval.ca}{p_end}