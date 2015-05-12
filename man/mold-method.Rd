\name{mold}
\alias{mold}
\alias{mold,IRanges-method}
\alias{mold,GRanges-method}
\alias{mold,GRangesList-method}
\alias{mold,Seqinfo-method}
\alias{mold,matrix-method}
\alias{mold,eSet-method}
\alias{mold,ExpressionSet-method}
\alias{mold,RangedSummarizedExperiment-method}
\alias{mold,SummarizedExperiment-method}
\alias{mold,Views-method}
\alias{mold,Rle-method}
\alias{mold,RleList-method}
\alias{mold,VRanges-method}
\title{mold data into data.frame}
\description{
  mold data into data.frame usued for visualization, so it's sometimes
  not as simple as coersion, for example, we add midpoint varialbe for
  mapping.
}
\usage{
\S4method{mold}{IRanges}(data)
\S4method{mold}{GRanges}(data)
\S4method{mold}{GRangesList}(data,indName = "grl_name")
\S4method{mold}{Seqinfo}(data)
\S4method{mold}{matrix}(data)
\S4method{mold}{eSet}(data)
\S4method{mold}{ExpressionSet}(data)
\S4method{mold}{RangedSummarizedExperiment}(data, assay.id = 1)
\S4method{mold}{Views}(data)
\S4method{mold}{Rle}(data)
\S4method{mold}{RleList}(data)
\S4method{mold}{VRanges}(data)
}
\arguments{
  \item{data}{
    data to be mold into data.frame with additional column that helped mapping.For example we add 'midpoint' variable to a ranged converted \link{data.frame}.
  }
  \item{indName}{
    when collapsing a \link{GRangesList}, list names are aggregated into extra column named by this parameter. 
  }
  \item{assay.id}{
    define the assay id you want to convert into a \link{data.frame}.
  }
}
\value{
  a data.frame object.
}
\author{Tengfei Yin}
