% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/path_finders.R
\name{status_message_path}
\alias{status_message_path}
\title{Find the status message CSV file path}
\usage{
status_message_path(paths)
}
\arguments{
\item{paths}{The path to one or more folders that contain pathogensurveillance output.}
}
\value{
character vector of length 1
}
\description{
Return the file path to the CSV with the status reports, warnings, and errors
for a given pathogensurveillance output folder.
}
\seealso{
Other path finders: 
\code{\link{busco_ref_path}()},
\code{\link{busco_tree_path}()},
\code{\link{considered_ref_meta_path}()},
\code{\link{core_ref_path}()},
\code{\link{core_tree_path}()},
\code{\link{estimated_ani_matrix_path}()},
\code{\link{pocp_matrix_path}()},
\code{\link{ref_meta_path}()},
\code{\link{report_group_path}()},
\code{\link{sample_meta_path}()},
\code{\link{selected_ref_meta_path}()},
\code{\link{sendsketch_path}()},
\code{\link{software_version_path}()},
\code{\link{variant_align_path}()},
\code{\link{variant_ref_path}()},
\code{\link{variant_tree_path}()}
}
\concept{path finders}
