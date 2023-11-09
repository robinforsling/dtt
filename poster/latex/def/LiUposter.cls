\ProvidesClass{LiUposter}[2009/10/21 LiU Poster Class v. 1.0]
%% Modification of the LINK-SIC-poster style [2015/09/08]
%% Modification of the rt-poster style [2006/09/26 RT Poster Class v. 1.0]
\NeedsTeXFormat{LaTeX2e}

%\newcommand*\logo@partner{partnerlogo.tex}

%% Pass options to mother class
\RequirePackage{xkeyval}
\newlength{\headh@ight}
\DeclareOptionX{headheight}{\setlength{\headh@ight}{#1}}
\DeclareOptionX{twocols}{\def\no@cols{2}}
\DeclareOptionX{threecols}{\def\no@cols{3}}
\DeclareOptionX{fourcols}{\def\no@cols{4}}
\DeclareOptionX{nocols}{\def\no@cols{0}}
\DeclareOptionX{a0paper}{\def\p@per{a0paper}}
\DeclareOptionX{a1paper}{\def\p@per{a1paper}}
\DeclareOptionX{landscape}{%
  \def\p@perorientation{landscape}%
  \setlength{\headh@ight}{103mm}}
\DeclareOptionX{portrait}{%
  \def\p@perorientation{portrait}%
  \setlength{\headh@ight}{140mm}}

\ExecuteOptionsX{fourcols,a0paper,landscape}

\ProcessOptionsX\relax

\LoadClass{article}

\RequirePackage{graphicx}

\RequirePackage{xcolor}
\RequirePackage{pgf}
\RequirePackage{calc}
\RequirePackage{geometry}
\RequirePackage{float}
\RequirePackage{ifthen}
\RequirePackage{tcolorbox}

\geometry{
  paper=\p@per,
  \p@perorientation,
  includeheadfoot,
  top=0pt, bottom=0pt,
  headheight=\headh@ight,
  headsep=10mm,
  footskip=80mm,
  right=40mm, left=40mm, % right=30mm, left=115mm
  includemp=false,
  marginpar=0pt
}

\def\oitem{\item[{\color{LiUorange}$\bullet$}]}
\def\gitem{\item[{\color{LiUgray}$\bullet$}]}
\def\pitem{\item[{\color{LiUpurple}$\bullet$}]}
\let\itemold\item
\renewcommand{\item}{\leftmargin3mm\labelsep5mm\parsep2mm\itemold}

%\usefonttheme{serif}
%\usepackage[no-math]{fontspec}
%\setsansfont{KorolevLiU}
%\setmainfont{Georgia}

%\usepackage{kpfonts}



%% Slide colors
\definecolor{LiUblue}{RGB}{0,185,231}
\definecolor{LiUturquoise}{RGB}{23,199,210}
\definecolor{LiUgreen}{RGB}{0,207,181}

\definecolor{LiUblue2}{RGB}{91,202,236}
\definecolor{LiUturquoise2}{RGB}{101,207,216}
\definecolor{LiUgreen2}{RGB}{86,214,192}

\definecolor{LiUblue3}{RGB}{143,215,241}
\definecolor{LiUturquoise3}{RGB}{150,218,225}
\definecolor{LiUgreen3}{RGB}{143,223,207}

\definecolor{LiUblue4}{RGB}{189,230,246}
\definecolor{LiUturquoise4}{RGB}{194,231,236}
\definecolor{LiUgreen4}{RGB}{191,235,225}

\definecolor{LiUblue5}{RGB}{210,238,249}
\definecolor{LiUturquoise5}{RGB}{215,239,242}
\definecolor{LiUgreen5}{RGB}{213,241,235}

\definecolor{LiUorange}{RGB}{255,100,66}
\definecolor{LiUpurple}{RGB}{137,129,211}
\definecolor{LiUyellow}{RGB}{253,239,93}
\definecolor{LiUgray}{RGB}{106,126,145}

\definecolor{LiUorange3}{RGB}{255,164,154}
\definecolor{LiUpurple3}{RGB}{183,180,225}
\definecolor{LiUyellow3}{RGB}{254,244,165}
\definecolor{LiUgray3}{RGB}{170,178,187}

\definecolor{LiUorange5}{RGB}{255,220,216}
\definecolor{LiUpurple5}{RGB}{226,224,242}
\definecolor{LiUyellow5}{RGB}{254,250,220}
\definecolor{LiUgray5}{RGB}{221,224,227}

\definecolor{myblue1}{rgb}{0,0.4470,0.7410}
\definecolor{myred1}{rgb}{0.6350,0.0780,0.1840}
\definecolor{mygreen1}{rgb}{0.4660,0.6740,0.1880}
\definecolor{myyellow1}{rgb}{0.9290,0.6940,0.1250}

\definecolor{myblue1}{rgb}{0,0.4470,0.7410}
\definecolor{myred1}{rgb}{0.6350,0.0780,0.1840}
\definecolor{mygreen1}{rgb}{0.4660,0.6740,0.1880}
\definecolor{myyellow1}{rgb}{0.9290,0.6940,0.1250}

%% Multicol stuff
\RequirePackage{multicol}
\setlength{\columnsep}{20mm}
\setlength{\columnseprule}{0mm}
\setlength{\multicolbaselineskip}{5pt}
\setlength{\multicolsep}{1ex}

\newlength\titlewidth
\setlength\titlewidth{\paperwidth-405mm} % 505
\define@key{LiUposter}{columns}{\def\no@cols{#1}}
%\define@key{LiUposter}{partnerlogo}{\def\logo@partner{#1}}
%\define@key{LiUposter}{nopartnerlogo}[true]{%
%  \renewcommand*\makepartnerlogo{\hspace*{80mm}\hfill}}
%\define@key{LiUposter}{rmpartnerlogo}[true]{%
%  \renewcommand*\makepartnerlogobox[1]{}%
%  \addtolength\titlewidth{202.5mm}% 252.5
%}
%\define@key{LiUposter}{rightlogo}{\def\logo@right{#1}}
%\define@key{LiUposter}{norightlogo}[true]{%
%  \renewcommand*\makerightlogo{\hfill\hspace*{80mm}}}
%\define@key{LiUposter}{rmrightlogo}[true]{%
%  \renewcommand*\makerightlogobox[1]{}%
%  \addtolength\titlewidth{202.5mm}% 252.5
%}
\newcommand\setupPoster[1]{\setkeys{LiUposter}{#1}}


%% Define some sections and stuff to be of a nice size
\renewcommand{\section}{%
  \@startsection {section}{1}{\z@}%
  {-3.5ex \@plus -1ex \@minus -.2ex}%
  {0.05ex \@plus.1ex}%
  {\sffamily\Large\color{black}}}
\renewcommand{\subsection}{%
  \@startsection{subsection}{2}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {0.9ex \@plus .2ex}%
  {\sffamily\Large\bfseries\color{black}}}
\renewcommand{\subsubsection}{%
  \@startsection{subsubsection}{3}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {0.9ex \@plus .2ex}%
  {\sffamily\large\bfseries\color{black}}}
\renewcommand{\paragraph}{%
  \@startsection{paragraph}{4}{\z@}%
  {3.25ex \@plus1ex \@minus.2ex}%
  {-1em}%
  {\normalfont\normalsize\bfseries\color{black}}}
\renewcommand{\subparagraph}{%
  \@startsection{subparagraph}{5}{\parindent}%
  {3.25ex \@plus1ex \@minus .2ex}%
  {-1em}%
  {\normalfont\normalsize\bfseries\color{black}}}
  
%  \newcommand*\makepartnerlogo{\hspace*{30mm}\includegraphics[height=50mm]{\logo@partner}\hfill}%
%\newcommand\makepartnerlogobox[1]{\begin{minipage}[c]{200mm}#1\end{minipage}}%

%% Define the head and foot of the poster
\newcommand*{\ps@poster}{%
  \renewcommand*{\@oddhead}{%
    \unskip\ignorespaces\noindent\hspace*{-2mm}%
    \colorbox{white}{%
    \begin{minipage}[b][\headh@ight]{\paperwidth-140mm}%100mm
      \unskip%
      \colorbox{white}{%
      {\begin{minipage}[t][\headh@ight-15mm]{\paperwidth-80mm} % -145mm
        \color{black}%
        %\centering%
        \vspace*{5mm}%
        \VeryHuge%
        {\hspace*{-2mm}\sffamily\mbox{}\@title\mbox{}}\\[4mm]%
        \Large%
        {\sffamily\mbox{}\@author\mbox{}}\\[-6mm]%
        \rule{\columnwidth}{1mm}
      \end{minipage}}}%
      \hfill\\%
    \end{minipage}\unskip}}%
\renewcommand*{\@evenhead}{\@evenhead}%
\renewcommand*{\@evenfoot}{%
\hspace*{-30mm}% \vspace*{-50mm}\\ % -110mm
 %\colorbox{red}{%
    \begin{minipage}[t][60mm]{\paperwidth-35mm}
    %\colorbox{blue}{%
    \hspace*{0mm}\begin{minipage}[b][80mm]{0.4\paperwidth}
    {\includegraphics[height=80mm]{LiU_primary_black.pdf}}\input{fig/partnerlogo.tex}
    \end{minipage}%}%
    \hfill
    %\colorbox{green}{%
    \begin{minipage}[b][80mm]{0.5\paperwidth} % 80mm
    \begin{flushright}
    {\sffamily\LARGE LINK\"OPING UNIVERSITY}\\[5mm] % 5mm
    {\sffamily\bfseries\LARGE \@institute}\\[18mm] % -1mm
    \ 
    \end{flushright}
    \end{minipage}%}
    \end{minipage}%}
 }
\renewcommand*{\@oddfoot}{\@evenfoot} }

\renewcommand*{\title}[1]{\gdef\@title{#1}}
\renewcommand*{\author}[1]{\gdef\@author{#1}}
\newcommand*{\institute}[1]{\gdef\@institute{#1}}

\input{def/LiUposter.fnt}

\renewcommand{\tiny}{\fontsize{17.28}{22}\selectfont}
\renewcommand{\scriptsize}{\fontsize{20.74}{25}\selectfont}
\renewcommand{\footnotesize}{\fontsize{24.88}{30}\selectfont}
\renewcommand{\small}{\fontsize{29.86}{37}\selectfont}
\renewcommand{\normalsize}{\fontsize{35.83}{45}\selectfont}
\renewcommand{\large}{\fontsize{43}{54}\selectfont}
\renewcommand{\Large}{\fontsize{51.6}{64}\selectfont}
\renewcommand{\LARGE}{\fontsize{61.92}{77}\selectfont}
\renewcommand{\huge}{\fontsize{74.3}{93}\selectfont}
\renewcommand{\Huge}{\fontsize{89.16}{112}\selectfont}
\newcommand{\HUGE}{\fontsize{110}{112}\selectfont}
\newcommand{\VeryHuge}{\fontsize{140}{134}\selectfont}
\normalsize

\setlength{\parindent}{0em}
\setlength{\smallskipamount}{6pt plus 2pt minus 2pt}
\setlength{\medskipamount}{12pt plus 4pt minus 4pt}
\setlength{\bigskipamount}{24pt plus 8pt minus 8pt}

\setlength{\abovecaptionskip}{25pt}
\setlength{\belowcaptionskip}{0pt}
\setlength{\abovedisplayskip}{25pt plus 6pt minus 15 pt}
\setlength{\abovedisplayshortskip}{0pt plus 6pt}
\setlength{\belowdisplayshortskip}{13pt plus 7pt minus 6pt}
\setlength{\belowdisplayskip}{\abovedisplayskip}

\setcounter{secnumdepth}{-2}
\setcounter{tocdepth}{-2}

\DeclareRobustCommand{\em}{%
  \@nomath\em \ifdim \fontdimen\@ne\font >\z@
  \upshape\else \slshape \fi\color{black}}

\ifnum\no@cols>1
  \AtBeginDocument{\pagestyle{poster}\begin{multicols}{\no@cols}}
  \AtEndDocument{\end{multicols}}
\else
  \AtBeginDocument{\pagestyle{poster}}
\fi

\newenvironment{bluebox}[1][]{%
\begin{tcolorbox}[colback=LiUblue5,coltitle=black,colbacktitle=LiUblue3,colframe=LiUblue,fonttitle=\sffamily\Large,title=#1,boxrule=0pt,arc=1mm,boxsep=3mm,toptitle=2mm,bottomtitle=2mm,oversize=-3mm]%width=1.05\columnwidth,
}{\end{tcolorbox}}

\newenvironment{greenbox}[1][]{%
\begin{tcolorbox}[colback=LiUgreen5,coltitle=black,colbacktitle=LiUgreen3,colframe=LiUgreen5,fonttitle=\sffamily\Large,title=#1,boxrule=0pt,arc=1mm,boxsep=3mm,toptitle=2mm,bottomtitle=2mm,oversize=-3mm]%width=1.05\columnwidth,
}{\end{tcolorbox}}

\newenvironment{turquoisebox}[1][]{%
\begin{tcolorbox}[colback=LiUturquoise5,coltitle=black,colbacktitle=LiUturquoise3,colframe=LiUturquoise5,fonttitle=\sffamily\Large,title=#1,boxrule=0pt,arc=1mm,boxsep=3mm,toptitle=2mm,bottomtitle=2mm,oversize=-3mm]%width=1.05\columnwidth,
}{\end{tcolorbox}}

