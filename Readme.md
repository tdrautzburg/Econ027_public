\# Econometrics \& Time Series Analysis: Replication Repository



This repository contains the code, data, and functions required to replicate key econometric models, including ADL models, HAC standard error visualizations, and monetary policy shock analysis.



---------------------------------------------------------

FOLDER STRUCTURE

---------------------------------------------------------



\* Code/scripts: Primary execution files (e.g., RR\_Monetary.R, SW\_OJ.R). 

&nbsp; These perform the main regressions and hypothesis testing.



\* Code/functions: Helper scripts containing custom functions, such as 

&nbsp; lag operators (L) and data cleaning routines.



\* Data: Raw data files used in the analysis (e.g., Stock \& Watson, 

&nbsp; Romer \& Romer series).



\* Results: All newly generated output from current runs, including 

&nbsp; log files and raw estimates.



\* Paper/GraphsTables: Finalized replication materials, including 

&nbsp; high-quality figures and formatted tables.



---------------------------------------------------------

USAGE

---------------------------------------------------------



1\. Environment: Ensure you have R or Stata installed.



2\. Working Directory: This repository uses the 'here' package to manage 

&nbsp;  paths. To ensure scripts run correctly:

&nbsp;  - Always open the script file directly (this usually sets the 

&nbsp;    initial directory).

&nbsp;  - The scripts will use here::i\_am("Code/scripts/filename.R") to 

&nbsp;    automatically find the project root.

&nbsp;  - Do not manually use setwd(), as the scripts handle pathing 

&nbsp;    relative to the main folder.



3\. Execution: Scripts in Code/scripts are intended to be run in sequence. 

&nbsp;  They are designed to source necessary functions from Code/functions 

&nbsp;  automatically.



---------------------------------------------------------

LICENSE (MIT) - NO WARRANTY

---------------------------------------------------------



Copyright (c) 2026 \[Thorsten Drautzburg/Swarthmore College]



Permission is hereby granted, free of charge, to any person obtaining a copy 

of this software and associated documentation files (the "Software"), to deal 

in the Software without restriction, including without limitation the rights 

to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 

copies of the Software, and to permit persons to whom the Software is 

furnished to do so, subject to the following conditions:



The above copyright notice and this permission notice shall be included in 

all copies or substantial portions of the Software.



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 

IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 

FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 

AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 

LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 

OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 

THE SOFTWARE.

