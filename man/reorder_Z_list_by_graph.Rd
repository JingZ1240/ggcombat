\name{reorder_Z_list_by_graph}
\alias{reorder_Z_list_by_graph}
\title{Reorder standardized site matrices by an estimated graph}
\usage{
reorder_Z_list_by_graph(Z_list, graph_fit)
}
\arguments{
\item{Z_list}{List of site-specific subject by feature matrices in original feature order.}
\item{graph_fit}{Graph object returned by \code{estimate_ggcombat_L_from_Z()}.}
}
\value{A list of matrices reordered by \code{graph_fit$results$Clist}.}
\description{
Convenience helper to align standardized data columns with the row order of the
membership matrix returned by \code{estimate_ggcombat_L_from_Z()}.
}
