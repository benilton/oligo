\name{NEWS}
\title{News for Package 'oligo'}

\section{Changes in version 1.32}{
        \subsection{USER VISIBLE CHANGES}{
                \itemize{
                       \item Fixed vignette
                       \item Added initial support to Generic Arrays
                       \item Using BiocGenerics definition of normalize
                }
        }
}

\section{Changes in version 1.24}{
	\subsection{BUG FIXES}{
		\itemize{
			\item Removed dependency on RConverters.h
		}
	}
}

\section{Changes in version 1.22}{
	\subsection{USER VISIBLE CHANGES}{
		\itemize{
			\item Improved results by getProbeInfo();
			\item fitPLM, coefs and resids are now Defunct. Use fitProbeLevelModel,
			coef and residuals respectively. 'coef' and 'residuals' follow the
			standards used elsewhere in R;
		}
	}
	\subsection{BUG FIXES}{
		\itemize{
			\item Fixed problem caused by the fact that oligoClasses had its own
			annotation() method when BiocGenerics added a new one;
			\item Several fixes to probe selector, allowing 'target' to be used;
			\item PA Calls didn't know about target;
			}
	}
}

\section{Changes in version 1.20}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item New getProbeInfo() function added, to simplify probe
  selection without using SQL.
      \item New fitProbeLevelModel() function added. It allows
      probe-level models ('plm' and 'medianpolish'), which can be used for
  QC.
      \item fitPLM, coefs and resids are now Deprecated. Use
  fitProbeLevelModel, coef and residuals respectively. 'coef' and
  'residuals' follow the standards used elsewhere in R.
      \item Now using foreach for parallelization.
      }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Addressed issue in which probes without chr info would be
  removed from pmChr, leading to results whose dimensions did not match
  PM matrix (on TilingFeatureSet).
    }
  }
}

\section{Changes in version 1.18}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item pmStrand() now available for Affy-Tiling arrays (by Kristof
      De Beuf);
      }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed message method when normalizing with large dataset
      support;
    }
  }
}

\section{Changes in version 1.16}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Using NEWS.Rd;
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New summarization methods are available: PLM, PLMr, PLMrr
      and PLMrc. These are implemented in the basicPLM function and also
      available through the summarize() interface (method=c('plm', 'plmr',
      'plmrr', 'plmrc'));
      \item fitPLM(), NUSE() and RLE() are now available;
      \item Presence/Absence calls are available through paCalls();
      \item MAplot method with advanced options;
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item 'pm<-', 'mm<-', 'bg<-' now work with ff_matrix replacements;
      \item boxplot() now uses the correct (default) ylim if transfo !=
      log2;
      \item Fixed documentation for hist(), boxplot(), fitPLM();
      \item Fixed image for PLMset on NimbleGen arrays;
    }
  }
}

\section{Changes in version 1.14}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Added support to read gzipped XYS files directly.
      \item Added information for citation.
      \item Added getNetAffx method to retrieve NetAffx annotation for
        Exon ST and Gene ST arrays.
      \item Reading in uncompressed XYS files is much faster than in 1.12.2.
      \item Added vignette for Exon/Gene ST arrays (rma and NetAffx annot)
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed bug in boxplot / hist / MAplot
      \item Fixed issue that would appear with GeneFeatureSet objects,
        in which it selections like obj[c(2, 2, 1),] would generate
        invalid GeneFeatureSet objects
      \item Fixed bugs in image / bg / bgindex / bgSequence
        for Exon and Gene arrays (due to changes in annotation)
      \item Added date/time of scan (for XYS/CEL files) to the
        protocolData slot. For 1 channel applications, the field
        is called 'dates'; for 2 channel apps, there are 2 fields:
        'dates1' and 'dates2'
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item New method, runDate(), to extract scan dates from
        FeatureSet objects
    }
  }
}
\section{Changes in version 1.11}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item New vignette structure:
        Vignette 1: Overview (with multiple examples)
        Vignette 2: Expression (NimbleGen data)
        Vignette 3: SNP Genotyping (Affymetrix)
      \item Uses Imports rather than Depends. Users must
        explicitly load Biobase, affyio, preprocessCore
        if they need functions from these packages.
      \item Revised man page for rma.
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Support to large datasets via ff package.
        The user should load the ff package in order to use
        this feature. Files used by ff are saved to getwd(),
        unless the user changes the standard path.
      \item Support to parallel computing via snow package.
        The user should load *both* ff and snow in order to use
        the parallel features implemented in the package
        (backgroundCorrect, normalize, summarize, rma,
        normalizeToTarget). Cluster variable should be set as
        options(cluster=makeCluster(...)).
      \item seqColors and darkColors, respectively, generate a
        vector of sequential colors (in blue) and (qualitative)
        dark colors (based onn Dark2 from RColorBrewer).
      \item plotM method for SnpSuperSet (plot genotype clusters).
      \item getAffinitySplineCoefficients estimates the affinity
        spline coefficients.
      \item getBaseProfile estimates and optionally plots sequence
        effect like gcrma.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Fixed definition for probeNames() for Exon/Gene arrays.
        (Mikhail Pachkov)
    }
  }
}
