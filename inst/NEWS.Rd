\name{NEWS}
\title{News for Package 'oligo'}

\section{Changes in series 1.15.x}{
  \subsection{USER VISIBLE CHANGES}{
    \itemize{
      \item Using NEWS.Rd
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
